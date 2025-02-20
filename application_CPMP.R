######################################################################################
#
# Purpose: Program code to compute the estimators in Application to CPMP study result
#
# Author: Akinori Nimura
# 
######################################################################################

library(tidyverse)
library(gt)
library(cmdstanr)

# Read dataset of aggregated data in CPMP which was disclosed in clinicaltrials.gov
dat <- read_csv("CPMP_study_result_asof_16Jun2024.csv")

dat %>% 
  arrange(disease, study) %>% 
  group_by(disease) %>% 
  slice(n()) %>%
  dplyr::select(study, disease) %>% 
  rename(interest = study) -> study.interest

# Calculated the standard deviation in each treatment arm of each sub-study
dat %>% 
  filter(Week == 8) %>% 
  full_join(study.interest) %>% 
  filter((study == interest & trtgrp != "diff") | (study != interest & trtgrp == "Placebo")) %>% 
  pivot_wider(names_from = paramcd, values_from = aval) %>% 
  mutate(
    se = (upper - lower) / 2 / qnorm(0.975), 
    sd = se * sqrt(n)
  ) %>% 
  arrange(disease, study) -> dat2

# Function to calculate the W related to the ridge tuning parameters
calc.w <- function(
    s1.vec, 
    s1.mat, 
    s1.mat.inv, 
    s2.vec, 
    s2.mat, 
    s2.mat.inv, 
    d, 
    est.alpha, 
    type = c("min.var", "prec.measure")
){
  if(type == "min.var"){
    w <- as.numeric(s1.mat.inv %*% solve(est.alpha %*% t(est.alpha) + solve(d)) %*% (solve(d) %*% s1.vec - s2.vec))
  }else{
    w <- as.numeric(s1.mat.inv %*% solve(est.alpha %*% t(est.alpha) + 2 * solve(d)) %*% (solve(d) %*% s1.vec - s2.vec))
  }
  
  while(any(w < 0 | w > 1)){
    id.zero <- which(w < 0)
    id.one <- which(w > 1)
    
    w[id.zero] <- 0
    w[id.one] <- 1
    
    d.any <- d[-c(id.one, id.zero), -c(id.one, id.zero)]
    est.alpha.any <- est.alpha[-c(id.one, id.zero)]
    s1.vec.any <- s1.vec[-c(id.one, id.zero)]
    one.any.vec <- rep(1, length(est.alpha.any))
    s1.mat.any <- diag(s1.vec.any, nrow = length(s1.vec.any))
    s1.mat.inv.any <- diag(1 / s1.vec.any, nrow = length(s1.vec.any))
    s2.vec.any <- s2.vec[-c(id.one, id.zero)]
    s2.mat.any <- diag(s2.vec.any, nrow = length(s2.vec.any))
    s2.mat.inv.any <- diag(1 / s2.vec.any, nrow = length(s2.vec.any))
    
    d.one <- d[id.one, id.one]
    est.alpha.one <- est.alpha[id.one]
    s1.vec.one <- s1.vec[id.one]
    one.one.vec <- rep(1, length(est.alpha.one))
    s1.mat.one <- diag(s1.vec.one, nrow = length(s1.vec.one))
    s2.vec.one <- s2.vec[id.one]
    s2.mat.one <- diag(s2.vec.one, nrow = length(s2.vec.one))
    
    if(length(w[-c(id.zero, id.one)]) > 0){
      if(type == "min.var"){
        w[-c(id.one, id.zero)] <- s1.mat.inv.any %*% solve(est.alpha.any %*% t(est.alpha.any) + solve(d.any)) %*% (solve(d.any) %*% s1.mat.any %*% one.any.vec - s2.mat.any %*% one.any.vec - est.alpha.any %*% t(est.alpha.one) %*% s1.mat.one %*% one.one.vec)
      }else{
        w[-c(id.one, id.zero)] <- s1.mat.inv.any %*% solve(est.alpha.any %*% t(est.alpha.any) + 2 * solve(d.any)) %*% (solve(d.any) %*% s1.mat.any %*% one.any.vec - s2.mat.any %*% one.any.vec - est.alpha.any %*% t(est.alpha.one) %*% s1.mat.one %*% one.one.vec)
      }
    }
  }
  w
}

# Function to calculate the derivative of W
deriv.w <- function(s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, d, est.alpha, type = "prec.measure", h = 1e-4){
  ret <- 1:length(est.alpha) %>% 
    lapply(function(id){
      est.alpha.h1 <- est.alpha.h2 <- est.alpha
      est.alpha.h1[id] <- est.alpha.h1[id] + h
      est.alpha.h2[id] <- est.alpha.h2[id] - h
      
      w1 <- calc.w(s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, d, est.alpha.h1, type)
      w2 <- calc.w(s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, d, est.alpha.h2, type)
      
      (w1 - w2) / (2 * h)
    }) %>% 
    do.call(what = cbind)
  
  ret
}

# Function to calculate the variance estimate in proposed ridge estimator
calc.ss.beta <- function(txx, tzx, tx.cov.x, h1, h2, tzm.cov.x, s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, p, d, est.alpha, type = "prec.measure", h = 1e-4){
  n.dim <- length(est.alpha)
  
  w <- calc.w(s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, d, est.alpha, type)
  w.min <- calc.w(s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, d, numeric(length(est.alpha)), "min.var")
  w.deriv <- deriv.w(s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, d, est.alpha, type, h = 1e-4)
  
  A1 <- solve(txx) %*% (tx.cov.x + t(tzx) %*% solve(h1) %*% h2 %*% solve(h1) %*% tzx - t(tzx) %*% solve(h1) %*% tzm.cov.x - t(tzm.cov.x) %*% solve(h1) %*% tzx) %*% solve(txx)
  A2 <- solve(txx) %*% t(tzx) %*% p %*% (diag(as.numeric(est.alpha), nrow = length(as.numeric(est.alpha))) %*% w.deriv + diag(w, nrow = length(w))) %*% t(p) %*% solve(h1) %*% (tzm.cov.x - h2 %*% tzx) %*% solve(txx)
  A3 <- solve(txx) %*% t(tzx) %*% p %*% diag(w, nrow = length(w)) %*% solve(d) %*% diag(w, nrow = length(w)) %*% t(p) %*% tzx %*% solve(txx)
  ret <- A1 + A2 + t(A2) + A3
  
  B1 <- t(tzx) %*% p %*% diag(1 - w.min, nrow = length(w)) %*% t(p) %*% solve(h1) %*% h2 %*% solve(h1) %*% p %*% diag(1 - w.min, nrow = length(w)) %*% t(p) %*% tzx
  B2 <- t(tzx) %*% p %*% diag(1 - w.min, nrow = length(w)) %*% t(p) %*% solve(h1) %*% tzm.cov.x
  
  ret.min <- solve(txx) %*% (tx.cov.x + B1 - B2 - t(B2)) %*% solve(txx)
  
  max(as.numeric(t(c(0, 1)) %*% ret %*% c(0, 1)), as.numeric(t(c(0, 1)) %*% ret.min %*% c(0, 1)))
}

# Function to calculate the restricted likelihood for variance estimate
restricted.likelihood <- function(p, mua.est, muc.est, mub.est, sa2.est, sc2.est, sb2.est, na, nc, nb){
  if(length(nb) == 1){
    tx.v.x <- (1 / p[1]) * matrix(
      c(
        na + nc, na, 0, 
        na, na, 0, 
        0, 0, 0
      ), 
      ncol = 3, 
      byrow = T
    ) + (1 / p[2]) * matrix(
      c(
        nb, 0, nb, 
        0, 0, 0, 
        nb, 0, nb
      ), 
      ncol = 3, 
      byrow = T
    )
    tx.v.y <- c(
      (nc * muc.est + na * mua.est) / p[1] + nb * mub.est / p[2], 
      na * mua.est / p[1], 
      nb * mub.est / p[2]
    )
    ty.v.y <- (nc * muc.est^2 + sc2.est + na * mua.est^2 + sa2.est) / p[1] + (nb * mub.est^2 + sb2.est) / p[2]
    logl <- -0.5 * (nc + na) * log(p[1]) - 0.5 * nb * log(p[2]) - 0.5 * log(det(tx.v.x)) - 0.5 * ty.v.y + 0.5 * as.numeric(t(tx.v.y) %*% solve(tx.v.x) %*% tx.v.y)
    return(logl)
  }else{
    tx.v.x <- (1 / p[1]) * matrix(
      c(
        na + nc, na, 0, 0, 
        na, na, 0, 0, 
        0, 0, 0, 0, 
        0, 0, 0, 0
      ), 
      ncol = 4, 
      byrow = T
    ) + (1 / p[2]) * matrix(
      c(
        nb[1], 0, nb[1], 0, 
        0, 0, 0, 0, 
        nb[1], 0, nb[1], 0, 
        0, 0, 0, 0
      ), 
      ncol = 4, 
      byrow = T
    ) + (1 / p[3]) * matrix(
      c(
        nb[2], 0, 0, nb[2], 
        0, 0, 0, 0,
        0, 0, 0, 0, 
        nb[2], 0, 0, nb[2]
      ), 
      ncol = 4, 
      byrow = T
    )
    tx.v.y <- c(
      (nc * muc.est + na * mua.est) / p[1] + nb[1] * mub.est[1] / p[2] + nb[2] * mub.est[2] / p[3], 
      na * mua.est / p[1], 
      nb[1] * mub.est[1] / p[2], 
      nb[2] * mub.est[2] / p[3]
    )
    ty.v.y <- (nc * muc.est^2 + sc2.est + na * mua.est^2 + sa2.est) / p[1] + (nb[1] * mub.est[1]^2 + sb2.est[1]) / p[2] + (nb[2] * mub.est[2]^2 + sb2.est[2]) / p[3]
    logl <- -0.5 * (nc + na) * log(p[1]) - 0.5 * nb[1] * log(p[2]) - 0.5 * nb[2] * log(p[3]) - 0.5 * log(det(tx.v.x)) - 0.5 * ty.v.y + 0.5 * as.numeric(t(tx.v.y) %*% solve(tx.v.x) %*% tx.v.y)
    return(logl)
  }
  
}

# Function to compute the proposed ridge estimator
ridge.estimator <- function(
    mua.est, 
    muc.est, 
    mub.est, 
    sea.est, 
    sec.est, 
    seb.est, 
    na, 
    nc, 
    nb, 
    cov.est = "reml", 
    type = "prec.measure"
){
  sa2.est <- sea.est^2 * na * (na - 1)
  sc2.est <- sec.est^2 * nc * (nc - 1)
  sb2.est <- seb.est^2 * nb * (nb - 1)
  
  if(length(nb) == 2){
    # z^t m z
    h1 <- matrix(
      c(
        nb[1] * (nc + nb[2]), - nb[1] * nb[2], 
        -nb[2] * nb[1], nb[2] * (nc + nb[1])
      ), 
      ncol = 2
    ) / (nc + nb[1] + nb[2])
    
    if(cov.est == "san.se"){
      tzm.cov.mz <- (nc + sum(nb))^(-2) * matrix(
        c(
          nb[1] ^ 2 * sc2.est + (nc + nb[2])^2 * sb2.est[1] + nb[1] ^ 2 * sb2.est[2], nb[1] * nb[2] * sc2.est - nb[2] * (nc + nb[2]) * sb2.est[1] - nb[1] * (nc + nb[1]) * sb2.est[2], 
          nb[1] * nb[2] * sc2.est - nb[2] * (nc + nb[2]) * sb2.est[1] - nb[1] * (nc + nb[1]) * sb2.est[2], nb[2] ^ 2 * sc2.est + nb[2] ^ 2 * sb2.est[1] + (nc + nb[1]) ^ 2 * sb2.est[2]
        ), 
        ncol = 2, 
        byrow = T
      )
      # z^t m cov(y) x
      tzm.cov.x <- matrix(
        c(
          -nb[1] / (nc + nb[1] + nb[2]) * sc2.est + (nc + nb[2]) / (nc + nb[1] + nb[2]) * sb2.est[1] - nb[1] / (nc + nb[1] + nb[2]) * sb2.est[2], 0, 
          -nb[1] / (nc + nb[1] + nb[2]) * sc2.est - nb[2] / (nc + nb[1] + nb[2]) * sb2.est[1] + (nc + nb[1]) / (nc + nb[1] + nb[2]) * sb2.est[2], 0
        ), 
        ncol = 2, 
        byrow = T
      )
      # x^t cov(y) x
      tx.cov.x <- matrix(
        c(
          sa2.est + sc2.est + sb2.est[1] + sb2.est[2], sa2.est, 
          sa2.est, sa2.est
        ), 
        ncol = 2, 
        byrow = T
      )
    }else{
      ret.sigma <- optim(
        c((sa2.est + sc2.est) / (na + nc - 2), sb2.est[1] / (nb[1] - 1), sb2.est[2] / (nb[2] - 1)), 
        restricted.likelihood, 
        method = "L-BFGS-B", 
        lower = c(0, 0, 0), 
        upper = c(Inf, Inf, Inf), 
        mua.est = mua.est, 
        muc.est = muc.est, 
        mub.est = mub.est,
        sa2.est = sa2.est, 
        sc2.est = sc2.est, 
        sb2.est = sb2.est, 
        na = na, 
        nc = nc, 
        nb = nb, 
        control = list(fnscale = -1)
      )
      tzm.cov.mz <- (nc + sum(nb))^(-2) * matrix(
        c(
          nb[1] ^ 2 * ret.sigma$par[1] * nc + (nc + nb[2])^2 * ret.sigma$par[2] * nb[1] + nb[1] ^ 2 * ret.sigma$par[3] * nb[2], nb[1] * nb[2] * ret.sigma$par[1] * nc - nb[2] * (nc + nb[2]) * ret.sigma$par[2] * nb[1] - nb[1] * (nc + nb[1]) * ret.sigma$par[3] * nb[2], 
          nb[1] * nb[2] * ret.sigma$par[1] * nc - nb[2] * (nc + nb[2]) * ret.sigma$par[2] * nb[1] - nb[1] * (nc + nb[1]) * ret.sigma$par[3] * nb[2], nb[2] ^ 2 * ret.sigma$par[1] * nc + nb[2] ^ 2 * ret.sigma$par[2] * nb[1] + (nc + nb[1]) ^ 2 * ret.sigma$par[3] * nb[2]
        ), 
        ncol = 2, 
        byrow = T
      )
      # z^t m cov(y) x
      tzm.cov.x <- matrix(
        c(
          -nb[1] / (nc + nb[1] + nb[2]) * ret.sigma$par[1] * nc + (nc + nb[2]) / (nc + nb[1] + nb[2]) * ret.sigma$par[2] * nb[1] - nb[1] / (nc + nb[1] + nb[2]) * ret.sigma$par[3] * nb[2], 0, 
          -nb[1] / (nc + nb[1] + nb[2]) * ret.sigma$par[1] * nc - nb[2] / (nc + nb[1] + nb[2]) * ret.sigma$par[2] * nb[1] + (nc + nb[1]) / (nc + nb[1] + nb[2]) * ret.sigma$par[3] * nb[2], 0
        ), 
        ncol = 2, 
        byrow = T
      )
      # x^t cov(y) x
      tx.cov.x <- matrix(
        c(
          ret.sigma$par[1] * (na + nc) + ret.sigma$par[2] * nb[1] + ret.sigma$par[3] * nb[2], ret.sigma$par[1] * na, 
          ret.sigma$par[1] * na, ret.sigma$par[1] * na
        ), 
        ncol = 2, 
        byrow = T
      )
    }
    # z^t x
    tzx <- matrix(
      c(
        nb[1], 0, 
        nb[2], 0
      ), 
      ncol = 2, 
      byrow = T
    )
    # x^t x
    txx <- matrix(
      c(
        na + nc + nb[1] + nb[2], na, 
        na, na
      ), 
      ncol = 2, 
      byrow = T
    )
    # x^t y
    txy <- c(
      mua.est * na + muc.est * nc + mub.est[1] * nb[1] + mub.est[2] * nb[2], 
      mua.est * na
    )
  }else if(length(nb) == 1){
    # z^t m z
    h1 <- matrix(
      c(nb * nc), 
      ncol = 1
    ) / (nc + nb)
    
    # z^t m cov(y) m z
    if(cov.est == "san.se"){
      tzm.cov.mz <- (nc + sum(nb))^(-2) * matrix(
        c(
          nb ^ 2 * sc2.est + nc^2 * sb2.est
        ), 
        ncol = 1, 
        byrow = T
      )
      # z^t m cov(y) x
      tzm.cov.x <- matrix(
        c(
          -nb / (nc + nb) * sc2.est + nc / (nc + nb) * sb2.est, 0
        ), 
        ncol = 2, 
        byrow = T
      )
      # x^t cov(y) x
      tx.cov.x <- matrix(
        c(
          sc2.est + sa2.est + sb2.est, sa2.est, 
          sa2.est, sa2.est
        ), 
        ncol = 2, 
        byrow = T
      )
    }else{
      ret.sigma <- optim(
        c((sa2.est + sc2.est) / (na + nc - 2), sb2.est / (nb - 1)), 
        restricted.likelihood, 
        method = "L-BFGS-B", 
        lower = c(1e-5, 1e-5), 
        upper = c(Inf, Inf), 
        mua.est = mua.est, 
        muc.est = muc.est, 
        mub.est = mub.est,
        sa2.est = sa2.est, 
        sc2.est = sc2.est, 
        sb2.est = sb2.est, 
        na = na, 
        nc = nc, 
        nb = nb, 
        control = list(fnscale = -1)
      )
      
      tzm.cov.mz <- (nc + sum(nb))^(-2) * matrix(
        c(
          nb ^ 2 * ret.sigma$par[1] * nc + nc ^ 2 * ret.sigma$par[2] * nb
        ), 
        ncol = 1, 
        byrow = T
      )
      # z^t m cov(y) x
      tzm.cov.x <- (nc + nb)^(-1) * matrix(
        c(
          -nb * ret.sigma$par[1] * nc + nc * ret.sigma$par[2] * nb, 0
        ), 
        ncol = 2, 
        byrow = T
      )
      # x^t cov(y) x
      tx.cov.x <- matrix(
        c(
          ret.sigma$par[1] * (na + nc) + ret.sigma$par[2] * nb, ret.sigma$par[1] * na, 
          ret.sigma$par[1] * na, ret.sigma$par[1] * na
        ), 
        ncol = 2, 
        byrow = T
      )
    }
    
    # z^t x
    tzx <- matrix(
      c(
        nb, 0
      ), 
      ncol = 2, 
      byrow = T
    )
    # x^t x
    txx <- matrix(
      c(
        na + nc + nb, na, 
        na, na
      ), 
      ncol = 2, 
      byrow = T
    )
    # x^t y
    txy <- c(
      mua.est * na + muc.est * nc + mub.est * nb, 
      mua.est * na
    )
  }
  
  ret.eigen <- eigen(h1 %*% solve(tzm.cov.mz) %*% h1)
  p <- ret.eigen$vectors
  d <- diag(ret.eigen$values, ncol = length(ret.eigen$values))
  
  est.beta <- c(muc.est, mua.est - muc.est)
  est.gamma <- mub.est - muc.est
  
  est.alpha <- t(p) %*% est.gamma
  
  s1.vec <- as.numeric(t(p) %*% tzx %*% solve(txx) %*% c(0, 1))
  s1.mat <- diag(s1.vec, nrow = length(s1.vec))
  s1.mat.inv <- diag(1 / s1.vec, nrow = length(s1.vec))
  s2.vec <- as.numeric(t(p) %*% solve(h1) %*% tzm.cov.x %*% solve(txx) %*% c(0, 1))
  s2.mat <- diag(s2.vec, nrow = length(s2.vec))
  s2.mat.inv <- diag(1 / s2.vec, nrow = length(s2.vec))
  
  w <- calc.w(s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, d, est.alpha, type)
  
  est.alpha.ast <- diag(1 - w, ncol = length(w)) %*% est.alpha
  est.gamma.ast <- p %*% est.alpha.ast
  est.beta.ast <- solve(txx) %*% txy - solve(txx) %*% t(tzx) %*% est.gamma.ast
  
  ss.beta.ast <- calc.ss.beta(txx, tzx, tx.cov.x, h1, tzm.cov.mz, tzm.cov.x, s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, p, d, est.alpha, type)
  
  tibble(
    parm = c("trtdif", paste0("dif_control[", seq(1, length(nb)), "]")), 
    est = c(est.beta.ast[2], as.numeric(est.gamma.ast)), 
    se = sqrt(c(ss.beta.ast, rep(NA, length(nb)))), 
    alpha = 0.05, 
    conf.level = 1 - alpha, 
    lower = est - se * qt(1 - alpha / 2, df = na + nc + sum(nb) - (2 + length(nb))), 
    upper = est + se * qt(1 - alpha / 2, df = na + nc + sum(nb) - (2 + length(nb))), 
    estimator = paste0("Ridge (", ifelse(cov.est == "san.se", "Sandwich, ", "REML, "), type)
  )
}

# Function to compute the full borrowing estimator
full.estimator <- function(
    mua.est, 
    muc.est, 
    mub.est, 
    sea.est, 
    sec.est, 
    seb.est, 
    na, 
    nc, 
    nb
){
  sa2.est <- sea.est^2 * na * (na - 1)
  sc2.est <- sec.est^2 * nc * (nc - 1)
  sb2.est <- seb.est^2 * nb * (nb - 1)
  
  if(length(nb) == 2){
    h1 <- matrix(
      c(
        nb[1] * (nc + nb[2]), - nb[1] * nb[2], 
        -nb[2] * nb[1], nb[2] * (nc + nb[1])
      ), 
      ncol = 2
    ) / (nc + nb[1] + nb[2])
    
    # z^t m
    tzm <- (nc + nb[1] + nb[2])^(-1) * matrix(
      c(
        nb[1], nc + nb[2], -nb[1], 
        nb[2], -nb[2], nc + nb[1]
      ), 
      ncol = 3, 
      byrow = T
    )
    
    h2 <- tzm %*% diag(c(sc2.est, sb2.est)) %*% t(tzm)
    
    # z^t x
    tzx <- matrix(
      c(
        nb[1], 0, 
        nb[2], 0
      ), 
      ncol = 2, 
      byrow = T
    )
    # x^t x
    txx <- matrix(
      c(
        na + nc + nb[1] + nb[2], na, 
        na, na
      ), 
      ncol = 2, 
      byrow = T
    )
    # x^t y
    txy <- c(
      mua.est * na + muc.est * nc + mub.est[1] * nb[1] + mub.est[2] * nb[2], 
      mua.est * na
    )
    
    tx.cov.x <- matrix(
      c(
        sa2.est + sc2.est + sb2.est[1] + sb2.est[2], sa2.est, 
        sa2.est, sa2.est
      ), 
      ncol = 2, 
      byrow = T
    )
  }else if(length(nb) == 1){
    h1 <- matrix(
      c(nb * nc), 
      ncol = 1
    ) / (nc + nb)
    
    # z^t m
    tzm <- (nc + nb)^(-1) * matrix(
      c(
        nb, nc
      ), 
      ncol = 2, 
      byrow = T
    )
    
    h2 <- tzm %*% diag(c(sc2.est, sb2.est)) %*% t(tzm)
    
    # z^t x
    tzx <- matrix(
      c(
        nb, 0
      ), 
      ncol = 2, 
      byrow = T
    )
    # x^t x
    txx <- matrix(
      c(
        na + nc + nb, na, 
        na, na
      ), 
      ncol = 2, 
      byrow = T
    )
    # x^t y
    txy <- c(
      mua.est * na + muc.est * nc + mub.est * nb, 
      mua.est * na
    )
    
    tx.cov.x <- matrix(
      c(
        sa2.est + sc2.est + sb2.est, sa2.est, 
        sa2.est, sa2.est
      ), 
      ncol = 2, 
      byrow = T
    )
  }
  est.beta <- solve(txx) %*% txy
  ss.beta1 <- as.numeric(t(c(0, 1)) %*% solve(txx) %*% tx.cov.x %*% solve(txx) %*% c(0, 1))
  
  tibble(
    parm = c("trtdif"), 
    est = est.beta[2], 
    se = sqrt(ss.beta1), 
    alpha = 0.05, 
    conf.level = 1 - alpha, 
    lower = est - se * qt(1 - alpha / 2, df = na + nc + sum(nb) - 2), 
    upper = est + se * qt(1 - alpha / 2, df = na + nc + sum(nb) - 2), 
    estimator = "Full borrow"
  )
}

# Function to compute the no borrowing estimator
noborrow.estimator <- function(
    mua.est, 
    muc.est, 
    mub.est, 
    sea.est, 
    sec.est, 
    seb.est, 
    na, 
    nc, 
    nb
){
  sa2.est <- sea.est^2 * na * (na - 1)
  sc2.est <- sec.est^2 * nc * (nc - 1)
  
  # x^t x
  txx <- matrix(
    c(
      na + nc, na, 
      na, na
    ), 
    ncol = 2, 
    byrow = T
  )
  # x^t y
  txy <- c(
    mua.est * na + muc.est * nc, 
    mua.est * na
  )
  
  tx.cov.x <- matrix(
    c(
      sa2.est + sc2.est, sa2.est, 
      sa2.est, sa2.est
    ), 
    ncol = 2, 
    byrow = T
  )
  
  est.beta <- solve(txx) %*% txy
  ss.beta1 <- as.numeric(t(c(0, 1)) %*% solve(txx) %*% tx.cov.x %*% solve(txx) %*% c(0, 1))
  
  tibble(
    parm = c("trtdif", paste0("dif_control[", seq(1, length(nb)), "]")), 
    est = c(est.beta[2], mub.est - muc.est), 
    se = c(sqrt(ss.beta1), rep(NA, length(nb))), 
    alpha = 0.05, 
    conf.level = 1 - alpha, 
    lower = est - se * qt(1 - alpha / 2, df = na + nc - 2), 
    upper = est + se * qt(1 - alpha / 2, df = na + nc - 2), 
    estimator = "No borrow"
  )
}

# Function to compute the power prior estimator
power.prior.eb.estimator <- function(
    mua.est, 
    muc.est, 
    mub.est, 
    sea.est, 
    sec.est, 
    seb.est, 
    na, 
    nc, 
    nb
){
  sa2.est <- sea.est^2 * na * (na - 1)
  sc2.est <- sec.est^2 * nc * (nc - 1)
  sb2.est <- seb.est^2 * nb * (nb - 1)
  
  txx <- nc
  txx0 <- nb
  txy <- muc.est * nc
  txy0 <- mub.est * nb
  ss.est <- (sea.est^2 * na * (na - 1) + sec.est^2 * nc * (nc - 1) + sum(seb.est^2 * nb * (nb - 1))) / (na + nc + sum(nb) - 2 - length(nb))
  
  alpha.opt <- optim(rep(0.5, length(nb)), marginal.a, method = "L-BFGS-B", lower = 1e-5, upper = 1, ss.est = ss.est, txx = txx, txx0 = txx0, txy = txy, txy0 = txy0)
  
  if(length(txx0) == 1){
    muc.est.updated <- solve(txx + alpha.opt$par * txx0) %*% (txy + alpha.opt$par * txy0)
    sec.est.updated <- sqrt(ss.est * solve(txx + alpha.opt$par * txx0))
  }else{
    muc.est.updated <- solve(txx + alpha.opt$par[1] * txx0[1] + alpha.opt$par[2] * txx0[2]) %*% (txy + alpha.opt$par[1] * txy0[1] + alpha.opt$par[2] * txy0[2])
    sec.est.updated <- sqrt(ss.est * solve(txx + alpha.opt$par[1] * txx0[1] + alpha.opt$par[2] * txx0[2]))
  }
  
  trtdif <- mua.est - muc.est.updated
  ss.beta1 <- sea.est ^ 2 + sec.est.updated^2
  
  tibble(
    parm = c("trtdif"), 
    est = as.numeric(trtdif), 
    se = as.numeric(sqrt(ss.beta1)), 
    alpha = 0.05, 
    conf.level = 1 - alpha, 
    lower = est - se * qnorm(1 - alpha / 2), 
    upper = est + se * qnorm(1 - alpha / 2), 
    estimator = "Power prior with EB"
  )
}

# Function to compute the marginal likelihood in power prior
marginal.a <- function(alpha, ss.est, txx, txx0, txy, txy0){
  if(length(txx0) == 1){
    log(txx + alpha * txx0) - log(txx0) - log(alpha) + (1 / ss.est) * (alpha * t(txy0) %*% solve(txx0) %*% txy0 - t(txy + alpha * txy0) %*% solve(txx + alpha * txx0) %*% (txy + alpha * txy0))
  }else{
    log(txx + alpha[1] * txx0[1] + alpha[2] * txx0[2]) - log(alpha[1] * txx0[1] + alpha[2] * txx0[2]) + (1 / ss.est) * (t(alpha[1] * txy0[1] + alpha[2] * txy0[2]) %*% solve(alpha[1] * txx0[1] + alpha[2] * txx0[2]) %*% (alpha[1] * txy0[1] + alpha[2] * txy0[2]) - t(txy + alpha[1] * txy0[1] + alpha[2] * txy0[2]) %*% solve(txx + alpha[1] * txx0[1] + alpha[2] * txx0[2]) %*% (txy + alpha[1] * txy0[1] + alpha[2] * txy0[2]))
  }
}

# Read stan code for hierarchical bayesian model
mod <- cmdstan_model("hbm_sum_stats_case_study.stan")

set.seed(12345)

# Compute each estimator
ret <- dat2 %>% 
  group_nest(disease) %>% 
  mutate(tmp = lapply(data, function(x){
    mua.est <- x$mean[x$study == x$interest & x$trtgrp != "Placebo"]
    muc.est <- x$mean[x$study == x$interest & x$trtgrp == "Placebo"]
    mub.est <- x$mean[x$study != x$interest]
    sea.est <- x$se[x$study == x$interest & x$trtgrp != "Placebo"]
    sec.est <- x$se[x$study == x$interest & x$trtgrp == "Placebo"]
    seb.est <- x$se[x$study != x$interest]
    
    na <- x$n[x$study == x$interest & x$trtgrp != "Placebo"]
    nc <- x$n[x$study == x$interest & x$trtgrp == "Placebo"]
    nb <- x$n[x$study != x$interest]
    
    sa2.est <- x$sd[x$study == x$interest & x$trtgrp != "Placebo"]^2 * (na - 1)
    sc2.est <- x$sd[x$study == x$interest & x$trtgrp == "Placebo"]^2 * (nc - 1)
    sb2.est <- x$sd[x$study != x$interest]^2 * (nb - 1)
    
    study.borrow <- x$study[x$study != x$interest]
    
    ret1.ridge.san.se.post <- ridge.estimator(
      mua.est, 
      muc.est, 
      mub.est, 
      sea.est, 
      sec.est, 
      seb.est, 
      na, 
      nc, 
      nb, 
      cov.est = "san.se", 
      type = "prec.measure"
    )
    ret1.ridge.reml.post <- ridge.estimator(
      mua.est, 
      muc.est, 
      mub.est, 
      sea.est, 
      sec.est, 
      seb.est, 
      na, 
      nc, 
      nb, 
      cov.est = "reml", 
      type = "prec.measure"
    )
    
    ret1.full <- full.estimator(
      mua.est, 
      muc.est, 
      mub.est, 
      sea.est, 
      sec.est, 
      seb.est, 
      na, 
      nc, 
      nb
    )
    
    ret1.power <- power.prior.eb.estimator(
      mua.est, 
      muc.est, 
      mub.est, 
      sea.est, 
      sec.est, 
      seb.est, 
      na, 
      nc, 
      nb
    )
    print(ret1.power)
    
    ret2.ridge.san.se.post <- ridge.estimator(
      mua.est, 
      muc.est, 
      mub.est[1], 
      sea.est, 
      sec.est, 
      seb.est[1], 
      na, 
      nc, 
      nb[1], 
      cov.est = "san.se", 
      type = "prec.measure"
    )
    ret2.ridge.reml.post <- ridge.estimator(
      mua.est, 
      muc.est, 
      mub.est[1], 
      sea.est, 
      sec.est, 
      seb.est[1], 
      na, 
      nc, 
      nb[1], 
      cov.est = "reml", 
      type = "prec.measure"
    )
    
    ret2.full <- full.estimator(
      mua.est, 
      muc.est, 
      mub.est[1], 
      sea.est, 
      sec.est, 
      seb.est[1], 
      na, 
      nc, 
      nb[1]
    )
    
    ret2.power <- power.prior.eb.estimator(
      mua.est, 
      muc.est, 
      mub.est[1], 
      sea.est, 
      sec.est, 
      seb.est[1], 
      na, 
      nc, 
      nb[1]
    )
    
    ret3.ridge.san.se.post <- ridge.estimator(
      mua.est, 
      muc.est, 
      mub.est[2], 
      sea.est, 
      sec.est, 
      seb.est[2], 
      na, 
      nc, 
      nb[2], 
      cov.est = "san.se", 
      type = "prec.measure"
    )
    ret3.ridge.reml.post <- ridge.estimator(
      mua.est, 
      muc.est, 
      mub.est[2], 
      sea.est, 
      sec.est, 
      seb.est[2], 
      na, 
      nc, 
      nb[2], 
      cov.est = "reml", 
      type = "prec.measure"
    )
    
    ret3.full <- full.estimator(
      mua.est, 
      muc.est, 
      mub.est[2], 
      sea.est, 
      sec.est, 
      seb.est[2], 
      na, 
      nc, 
      nb[2]
    )
    
    ret3.power <- power.prior.eb.estimator(
      mua.est, 
      muc.est, 
      mub.est[2], 
      sea.est, 
      sec.est, 
      seb.est[2], 
      na, 
      nc, 
      nb[2]
    )
    
    ret1.no <- noborrow.estimator(
      mua.est, 
      muc.est, 
      mub.est, 
      sea.est, 
      sec.est, 
      seb.est, 
      na, 
      nc, 
      nb
    )
    
    ret2.no <- noborrow.estimator(
      mua.est, 
      muc.est, 
      mub.est[1], 
      sea.est, 
      sec.est, 
      seb.est[1], 
      na, 
      nc, 
      nb[1]
    )
    
    ret3.no <- noborrow.estimator(
      mua.est, 
      muc.est, 
      mub.est[2], 
      sea.est, 
      sec.est, 
      seb.est[2], 
      na, 
      nc, 
      nb[2]
    )
    
    bind_rows(
      bind_rows(
        ret1.no, 
        ret1.ridge.san.se.post, 
        ret1.ridge.reml.post, 
        ret1.full, 
        ret1.power
      ) %>% 
        mutate(
          study.borrow = str_c(study.borrow[1], study.borrow[2], sep = " & ")
        ), 
      bind_rows(
        ret2.no, 
        ret2.ridge.san.se.post, 
        ret2.ridge.reml.post, 
        ret2.full, 
        ret2.power
      ) %>% 
        mutate(
          study.borrow = study.borrow[1]
        ), 
      bind_rows(
        ret3.no, 
        ret3.ridge.san.se.post, 
        ret3.ridge.reml.post, 
        ret3.full, 
        ret3.power
      ) %>% 
        mutate(
          study.borrow = study.borrow[2]
        )
    ) -> result
    
    set.seed(12345)
    
    fit.mod1_05 <- mod$sample(
      data = list(
        mua_est = mua.est, 
        muc_est = muc.est, 
        mub_est = mub.est, 
        ssa_est = sea.est, 
        ssc_est = sec.est, 
        ssb_est = seb.est, 
        n_sub = 2, 
        s_tau = 0.5
      ), 
      iter_sampling = 1000, 
      init = function(){
        list(
          mua = rnorm(1, mua.est, sea.est), 
          mu_sub = c(rnorm(1, muc.est, sec.est), rnorm(1, mub.est[1], seb.est[1]), rnorm(1, mub.est[2], seb.est[2])), 
          mu = rnorm(1, (muc.est + mub.est[1] + mub.est[2]) / 3, 0.5), 
          tau = runif(1, 0, 0.5)
        )
      }
    )
    fit.mod1_1 <- mod$sample(
      data = list(
        mua_est = mua.est, 
        muc_est = muc.est, 
        mub_est = mub.est, 
        ssa_est = sea.est, 
        ssc_est = sec.est, 
        ssb_est = seb.est, 
        n_sub = 2, 
        s_tau = 1.0
      ), 
      iter_sampling = 1000, 
      init = function(){
        list(
          mua = rnorm(1, mua.est, sea.est), 
          mu_sub = c(rnorm(1, muc.est, sec.est), rnorm(1, mub.est[1], seb.est[1]), rnorm(1, mub.est[2], seb.est[2])), 
          mu = rnorm(1, (muc.est + mub.est[1] + mub.est[2]) / 3, 1.0), 
          tau = runif(1, 0, 1.0)
        )
      }
    )
    fit.mod2_05 <- mod$sample(
      data = list(
        mua_est = mua.est, 
        muc_est = muc.est, 
        mub_est = mub.est[1], 
        ssa_est = sea.est, 
        ssc_est = sec.est, 
        ssb_est = seb.est[1], 
        n_sub = 1, 
        s_tau = 0.5
      ), 
      iter_sampling = 1000, 
      init = function(){
        list(
          mua = rnorm(1, mua.est, sea.est), 
          mu_sub = c(rnorm(1, muc.est, sec.est), rnorm(1, mub.est[1], seb.est[1])), 
          mu = rnorm(1, (muc.est + mub.est[1]) / 2, 0.5), 
          tau = runif(1, 0, 0.5)
        )
      }
    )
    fit.mod2_1 <- mod$sample(
      data = list(
        mua_est = mua.est, 
        muc_est = muc.est, 
        mub_est = mub.est[1], 
        ssa_est = sea.est, 
        ssc_est = sec.est, 
        ssb_est = seb.est[1], 
        n_sub = 1, 
        s_tau = 1.0
      ), 
      iter_sampling = 1000, 
      init = function(){
        list(
          mua = rnorm(1, mua.est, sea.est), 
          mu_sub = c(rnorm(1, muc.est, sec.est), rnorm(1, mub.est[1], seb.est[1])), 
          mu = rnorm(1, (muc.est + mub.est[1]) / 2, 1.0), 
          tau = runif(1, 0, 1.0)
        )
      }
    )
    fit.mod3_05 <- mod$sample(
      data = list(
        mua_est = mua.est, 
        muc_est = muc.est, 
        mub_est = mub.est[2], 
        ssa_est = sea.est, 
        ssc_est = sec.est, 
        ssb_est = seb.est[2], 
        n_sub = 1, 
        s_tau = 0.5
      ), 
      iter_sampling = 1000, 
      init = function(){
        list(
          mua = rnorm(1, mua.est, sea.est), 
          mu_sub = c(rnorm(1, muc.est, sec.est), rnorm(1, mub.est[2], seb.est[2])), 
          mu = rnorm(1, (muc.est + mub.est[2]) / 2, 0.5), 
          tau = runif(1, 0, 0.5)
        )
      }
    )
    fit.mod3_1 <- mod$sample(
      data = list(
        mua_est = mua.est, 
        muc_est = muc.est, 
        mub_est = mub.est[2], 
        ssa_est = sea.est, 
        ssc_est = sec.est, 
        ssb_est = seb.est[2], 
        n_sub = 1, 
        s_tau = 1.0
      ), 
      iter_sampling = 1000, 
      init = function(){
        list(
          mua = rnorm(1, mua.est, sea.est), 
          mu_sub = c(rnorm(1, muc.est, sec.est), rnorm(1, mub.est[2], seb.est[2])), 
          mu = rnorm(1, (muc.est + mub.est[2]) / 2, 1.0), 
          tau = runif(1, 0, 1.0)
        )
      }
    )
    result %>% 
      bind_rows(
        bind_rows(
          fit.mod1_05$summary(c("dif", "dif_control"), "mean", "sd", ~quantile(., probs = c(0.025, 0.975))) %>% 
            rename(
              est = mean, 
              se = sd, 
              lower = "2.5%", 
              upper = "97.5%"
            ) %>% 
            mutate(
              alpha = 0.05, 
              estimator = "HBM", 
              s_tau = 0.5, 
              parm = case_when(
                variable == "dif" ~ "trtdif", 
                variable == "dif_control" ~ "dif_control", 
                T ~ variable
              )
            ) %>% 
            dplyr::select(-variable), 
          fit.mod1_1$summary(c("dif", "dif_control"), "mean", "sd", ~quantile(., probs = c(0.025, 0.975))) %>% 
            rename(
              est = mean, 
              se = sd, 
              lower = "2.5%", 
              upper = "97.5%"
            ) %>% 
            mutate(
              alpha = 0.05, 
              estimator = "HBM", 
              s_tau = 1.0, 
              parm = case_when(
                variable == "dif" ~ "trtdif", 
                variable == "dif_control" ~ "dif_control", 
                T ~ variable
              )
            ) %>% 
            dplyr::select(-variable)
        ) %>% 
          mutate(
            study.borrow = str_c(study.borrow[1], study.borrow[2], sep = " & ")
          ), 
        bind_rows(
          fit.mod2_05$summary(c("dif", "dif_control"), "mean", "sd", ~quantile(., probs = c(0.025, 0.975))) %>% 
            rename(
              est = mean, 
              se = sd, 
              lower = "2.5%", 
              upper = "97.5%"
            ) %>% 
            mutate(
              alpha = 0.05, 
              estimator = "HBM", 
              s_tau = 0.5, 
              parm = case_when(
                variable == "dif" ~ "trtdif", 
                variable == "dif_control" ~ "dif_control", 
                T ~ variable
              )
            ) %>% 
            dplyr::select(-variable), 
          fit.mod2_1$summary(c("dif", "dif_control"), "mean", "sd", ~quantile(., probs = c(0.025, 0.975))) %>% 
            rename(
              est = mean, 
              se = sd, 
              lower = "2.5%", 
              upper = "97.5%"
            ) %>% 
            mutate(
              alpha = 0.05, 
              estimator = "HBM", 
              s_tau = 1.0, 
              parm = case_when(
                variable == "dif" ~ "trtdif", 
                variable == "dif_control" ~ "dif_control", 
                T ~ variable
              )
            ) %>% 
            dplyr::select(-variable)
        ) %>% 
          mutate(
            study.borrow = study.borrow[1]
          ), 
        bind_rows(
          fit.mod3_05$summary(c("dif", "dif_control"), "mean", "sd", ~quantile(., probs = c(0.025, 0.975))) %>% 
            rename(
              est = mean, 
              se = sd, 
              lower = "2.5%", 
              upper = "97.5%"
            ) %>% 
            mutate(
              alpha = 0.05, 
              estimator = "HBM", 
              s_tau = 0.5, 
              parm = case_when(
                variable == "dif" ~ "trtdif", 
                variable == "dif_control" ~ "dif_control", 
                T ~ variable
              )
            ) %>% 
            dplyr::select(-variable), 
          fit.mod3_1$summary(c("dif", "dif_control"), "mean", "sd", ~quantile(., probs = c(0.025, 0.975))) %>% 
            rename(
              est = mean, 
              se = sd, 
              lower = "2.5%", 
              upper = "97.5%"
            ) %>% 
            mutate(
              alpha = 0.05, 
              estimator = "HBM", 
              s_tau = 1.0, 
              parm = case_when(
                variable == "dif" ~ "trtdif", 
                variable == "dif_control" ~ "dif_control", 
                T ~ variable
              )
            ) %>% 
            dplyr::select(-variable)
        ) %>% 
          mutate(
            study.borrow = study.borrow[2]
          )
      )
  }))

