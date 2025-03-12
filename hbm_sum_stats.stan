data {
  real mua_est;
  real ssa_est;
  real muc_est;
  real ssc_est;
  int n_sub;
  array[n_sub] real mub_est;
  array[n_sub] real ssb_est;
  real s_tau;
}

parameters {
  real mu; 
  array[n_sub+1] real mu_sub; 
  real mua; 
  real<lower=0> tau; 
}

model {
  tau ~ uniform(0, s_tau);
  // tau ~ normal(0, s_tau); 
  // tau ~ inv_gamma(0.4, 0.002);
  mu ~ normal(0, 100);
  mu_sub ~ normal(mu, tau); 
  mua ~ normal(0, 100); 
  
  mua_est ~ normal(mua, ssa_est);
  muc_est ~ normal(mu_sub[1], ssc_est);
  mub_est ~ normal(mu_sub[2], ssb_est);
}

generated quantities {
  real dif;
  dif = mua - mu_sub[1];
}
