data {
  int<lower=1> N; // number of treatments
  int<lower=1> I; // number of strains
  
  int<lower=1,upper=I> sid[N]; // strain indices
  
  int<lower=0> x[N];     // dose in ob
  int<lower=0> y[N];     // dependent variable, number of virus-killed
  int<lower=0> total[N]; // total caterpillars in treatment
  
  real prior_mu_alpha;
  real prior_mu_beta;
  
  real<lower=0> prior_sigma_alpha;
  real<lower=0> prior_sigma_beta;
}

parameters{
  real raw_alpha[I];
  real raw_beta[I];
  
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
}

transformed parameters {
  real alpha[I];
  real beta[I];
  
  vector[N] theta; //predicted proportion virus-killed on logit scale
  vector[N] inv_logit_theta;
  
  for(i in 1:I) {
    alpha[i] = prior_mu_alpha + raw_alpha[i] * sigma_alpha;
    beta[i] = prior_mu_beta + raw_beta[i] * sigma_beta;
  }

  for(n in 1:N) {
    theta[n] = alpha[sid[n]] + beta[sid[n]] * x[n];
  }
  
  inv_logit_theta = inv_logit(theta);
}

model {
  //priors
  sigma_alpha ~ normal(prior_sigma_alpha, prior_sigma_alpha*3);
  sigma_beta ~ normal(prior_sigma_beta, prior_sigma_beta*3);
  
  raw_alpha ~ normal(0,1);
  raw_beta ~ normal(0,1);
  

  //likelihood
  y ~ binomial_logit(total, theta);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_logit_lpmf(y[n] | total[n], theta[n]);
  }
}
