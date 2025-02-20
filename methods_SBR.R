#################################################################################
#
# Purpose: Program code to compute the estimators used for the simulation study
#
# Author: Akinori Nimura
# 
#################################################################################


library(tidyverse)
library(MASS)
library(cmdstanr)

# Calculate W related to ridge tuning parameters
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

# Calculate the derivative W
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

# Calculate the variance estimate of proposed ridge estimator
calc.ss.beta <- function(x, z, m.x, cov.y, s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, p, d, est.alpha, type = "prec.measure", h = 1e-4){
  n.dim <- length(est.alpha)
  
  w <- calc.w(s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, d, est.alpha, type)
  w.min <- calc.w(s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, d, numeric(length(est.alpha)), "min.var")
  w.deriv <- deriv.w(s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, d, est.alpha, type, h = 1e-4)
  
  A <- t(diag(nrow(x)) - m.x %*% z %*% solve(t(z) %*% m.x %*% z) %*% t(z)) %*% cov.y %*% (diag(nrow(x)) - m.x %*% z %*% solve(t(z) %*% m.x %*% z) %*% t(z))
  B <- z %*% p %*% (diag(as.numeric(est.alpha), nrow = length(as.numeric(est.alpha))) %*% w.deriv + diag(w, nrow = length(w))) %*% t(p) %*% solve(t(z) %*% m.x %*% z) %*% t(z) %*% m.x %*% cov.y %*% (diag(nrow(x)) - m.x %*% z %*% solve(t(z) %*% m.x %*% z) %*% t(z))
  AB.min <- (diag(nrow(x)) - m.x %*% z %*% solve(t(z) %*% m.x %*% z) %*% p %*% diag(1 - w.min, nrow = length(w.min)) %*% t(p) %*% t(z))
  
  #ret <- solve(t(x) %*% x) %*% t(x) %*% (A + B + t(B)) %*% x %*% solve(t(x) %*% x)
  ret <- solve(t(x) %*% x) %*% t(x) %*% (A + B + t(B) + z %*% p %*% diag(w, nrow = length(w)) %*% solve(d) %*% diag(w, nrow = length(w)) %*% t(p) %*% t(z)) %*% x %*% solve(t(x) %*% x)
  ret.min <- solve(t(x) %*% x) %*% t(x) %*% t(AB.min) %*% cov.y %*% AB.min %*% x %*% solve(t(x) %*% x)
  
  max(as.numeric(t(c(0, 1)) %*% ret %*% c(0, 1)), as.numeric(t(c(0, 1)) %*% ret.min %*% c(0, 1)))
}

# Calculate REML estimate of residual variance
calc.cov.y <- function(m.z, y, z){
  tmp <- eigen(m.z)
  a <- tmp$vectors[, abs(tmp$values) > 1e-5]
  resid <- as.numeric(m.z %*% y)
  
  ss.init <- c(sum(resid[z[, 1] == 0]^2) / (sum(1-z[, 1]) - 2), sum(resid[z[, 1] == 1]^2) / (sum(z[, 1]) - 1))
  ay <- t(a) %*% y
  a1 <- a[seq(1, sum(1-z[,1])), ]
  a2 <- a[seq(sum(1-z[,1]) + 1, nrow(z)), ]
  
  g <- function(ss){
    log(det(ss[1] * t(a1) %*% a1 + ss[2] * t(a2) %*% a2)) + t(ay) %*% solve(ss[1] * t(a1) %*% a1 + ss[2] * t(a2) %*% a2) %*% ay
  }
  
  ss <- nlm(g, ss.init)
  
  diag(
    c(
      rep(ss$estimate[1], sum(1-z[,1])), 
      rep(ss$estimate[2], sum(z[,1]))
    )
  )
}

# Compute the proposed ridge estimator
ridge.estimator <- function(dat, cov.est = "reml", type = "prec.measure"){
  m.x <- diag(nrow(dat)) - dat$x %*% solve(t(dat$x) %*% dat$x) %*% t(dat$x)
  m.z <- m.x - m.x %*% dat$z %*% solve(t(dat$z) %*% m.x %*% dat$z) %*% t(dat$z) %*% m.x
  resid <- as.numeric(m.z %*% dat$y)
  
  if(cov.est == "san.se"){
    cov.y <- diag(resid^2)
  }else{
    cov.y <- calc.cov.y(m.z, dat$y, dat$z)
  }
  
  h1 <- t(dat$z) %*% m.x %*% dat$z
  h2 <- t(dat$z) %*% m.x %*% cov.y %*% m.x %*% dat$z
  est.gamma <- solve(t(dat$z) %*% m.x %*% dat$z) %*% t(dat$z) %*% m.x %*% dat$y
  eigen.ret <- eigen(h1 %*% solve(h2) %*% h1)
  p <- eigen.ret$vectors
  d <- diag(eigen.ret$values, nrow = length(eigen.ret$values))
  est.alpha <- t(p) %*% est.gamma
  s1.vec <- as.numeric(t(p) %*% t(dat$z) %*% dat$x %*% solve(t(dat$x) %*% dat$x) %*% c(0, 1))
  s1.mat <- diag(s1.vec, nrow = length(s1.vec))
  s1.mat.inv <- diag(1 / s1.vec, nrow = length(s1.vec))
  s2.vec <- as.numeric(t(p) %*% solve(h1) %*% t(dat$z) %*% m.x %*% cov.y %*% dat$x %*% solve(t(dat$x) %*% dat$x) %*% c(0, 1))
  s2.mat <- diag(s2.vec, nrow = length(s2.vec))
  s2.mat.inv <- diag(1 / s2.vec, nrow = length(s2.vec))
  
  w <- calc.w(s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, d, est.alpha, type = "prec.measure")
  est.alpha.ast <- diag(1 - w, nrow = length(w)) %*% est.alpha
  
  est.beta.ast <- solve(t(dat$x) %*% dat$x) %*% t(dat$x) %*% (dat$y - dat$z %*% p %*% est.alpha.ast)
  ss.beta.ast <- calc.ss.beta(dat$x, dat$z, m.x, cov.y, s1.vec, s1.mat, s1.mat.inv, s2.vec, s2.mat, s2.mat.inv, p, d, est.alpha, type)
  
  tibble(
    method = paste0("Ridge (", cov.est, ")"), 
    hyp.parm = NA, 
    trtdif = est.beta.ast[2], 
    se = as.numeric(sqrt(ss.beta.ast)), 
    alpha = 0.05, 
    lower = trtdif - qt(1 - alpha / 2, df = sum(diag(m.z))) * se, 
    upper = trtdif + qt(1 - alpha / 2, df = sum(diag(m.z))) * se, 
    w = w, 
    dif.study = as.numeric(est.gamma)
  )
}

# Compute the full borrowing estimator
full.borrow <- function(dat){
  est.beta <- solve(t(dat$x) %*% dat$x) %*% t(dat$x) %*% dat$y
  resid <- as.numeric(dat$y - dat$x %*% est.beta)
  ss.beta <- solve(t(dat$x) %*% dat$x) %*% t(dat$x) %*% diag(resid ^ 2) %*% dat$x %*% solve(t(dat$x) %*% dat$x)
  
  tibble(
    method = "Full borrowing", 
    hyp.parm = NA, 
    trtdif = est.beta[2], 
    se = sqrt(ss.beta[2, 2]), 
    alpha = 0.05, 
    lower = trtdif - qt(1 - alpha / 2, df = nrow(dat) - 2) * se, 
    upper = trtdif + qt(1 - alpha / 2, df = nrow(dat) - 2) * se
  )
}

# Compute the no borrowing estimator
no.borrow <- function(dat){
  ret <- lm(y ~ x - 1, data = dat %>% filter(z[, 1] == 0))
  tibble(
    method = "No borrow", 
    hyp.parm = NA, 
    trtdif = coef(ret)[2], 
    se = sqrt(vcov(ret)[2, 2]), 
    alpha = 0.05, 
    lower = trtdif - qt(1 - alpha / 2, df = ret$df.residual) * se, 
    upper = trtdif + qt(1 - alpha / 2, df = ret$df.residual) * se
  )
}

# Compute the power prior estimator
power.prior.eb <- function(dat){
  xc <- 1 - (dat %>% filter(x[, 2] == 0 & z[, 1] == 0) %>% .$x %>% .[, 2])
  yc <- dat %>% filter(x[, 2] == 0 & z[, 1] == 0) %>% .$y
  xb <- dat %>% filter(x[, 2] == 0 & z[, 1] == 1) %>% .$z
  yb <- dat %>% filter(x[, 2] == 0 & z[, 1] == 1) %>% .$y
  ya <- dat %>% filter(x[, 2] == 1 & z[, 1] == 0) %>% .$y
  na <- length(ya)
  nc <- length(yc)
  nb <- length(yb)
  
  txx <- as.numeric(t(xc) %*% xc)
  txx0 <- as.numeric(t(xb) %*% xb)
  txy <- as.numeric(t(xc) %*% yc)
  txy0 <- as.numeric(t(xb) %*% yb)
  ss.est <- (var(ya) * (na - 1) + var(yc) * (nc - 1) + var(yb) * (nb - 1)) / (na + nc + nb - 3)
  
  alpha.opt <- optim(
    rep(0.5, length(nb)), 
    marginal.a, 
    method = "L-BFGS-B", 
    lower = 1e-5, 
    upper = 1, 
    ss.est = ss.est, 
    txx = txx,
    txx0 = txx0,
    txy = txy, 
    txy0 = txy0
  )
  
  muc.est.updated <- solve(txx + alpha.opt$par * txx0) %*% (txy + alpha.opt$par * txy0)
  sec.est.updated <- sqrt(ss.est * solve(txx + alpha.opt$par * txx0))
  
  trtdif <- mean(ya) - muc.est.updated
  ss.beta1 <- ss.est / na + sec.est.updated^2
  
  tibble(
    trtdif = as.numeric(trtdif), 
    se = as.numeric(sqrt(ss.beta1)), 
    alpha = 0.05, 
    conf.level = 1 - alpha, 
    lower = trtdif - se * qnorm(1 - alpha / 2), 
    upper = trtdif + se * qnorm(1 - alpha / 2), 
    method = "Power prior with EB"
  )
}

# Calculate the marginal likelihood in the power prior
marginal.a <- function(alpha, ss.est, txx, txx0, txy, txy0){
  log(txx + alpha * txx0) - log(txx0) - log(alpha) + (1 / ss.est) * (alpha * t(txy0) %*% solve(txx0) %*% txy0 - t(txy + alpha * txy0) %*% solve(txx + alpha * txx0) %*% (txy + alpha * txy0))
}

# Read stan code for hierarchical bayesian model
mod <- cmdstan_model("hbm_sum_stats.stan")

# Compute the hierarchical bayesian model estimator
hbm.sum.stats <- function(dat, s_tau = 1.0, n_sample = 1000){
  yc <- dat %>% filter(x[, 2] == 0 & z[, 1] == 0) %>% .$y
  yb <- dat %>% filter(x[, 2] == 0 & z[, 1] == 1) %>% .$y
  ya <- dat %>% filter(x[, 2] == 1 & z[, 1] == 0) %>% .$y
  mua_est <- mean(ya)
  muc_est <- mean(yc)
  mub_est <- mean(yb)
  ssa_est <- sqrt(var(ya) / length(ya))
  ssc_est <- sqrt(var(yc) / length(yc))
  ssb_est <- sqrt(var(yb) / length(yb))
  
  fit.mod <- mod$sample(
    data = list(
      mua_est = mua_est, 
      ssa_est = ssa_est, 
      muc_est = muc_est, 
      ssc_est = ssc_est, 
      n_sub = 1, 
      mub_est = mub_est, 
      ssb_est = ssb_est, 
      s_tau = s_tau
    ), 
    init = function(){
      list(
        mua = rnorm(1, mua_est, ssa_est), 
        mu_sub = c(rnorm(1, muc_est, ssc_est), rnorm(1, mub_est, ssb_est)), 
        mu = rnorm(1, (muc_est + mub_est) / 2, s_tau), 
        tau = runif(1, 0, s_tau)
      )
    }, 
    iter_sampling = n_sample
  )
  fit.mod$summary(
    "dif", 
    mean, 
    sd, 
    ~quantile(., probs = c(seq(0.01, 0.1, 0.001), 0.975))
  ) %>% 
    rename(
      trtdif = mean, 
      se = sd, 
      lower = "2.5%", 
      upper = "97.5%"
    ) %>% 
    mutate(
      alpha = 0.05, 
      method = paste0("HBM (s_tau=", s_tau, ")"), 
      hyp.parm = NA
    ) %>% 
    dplyr::select(-variable)
}
