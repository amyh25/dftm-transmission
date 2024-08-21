data {
  int<lower=1> N; // number of treatments
  int<lower=1> H; // number of capsids
  int<lower=1> I; // number of strains
  
  int<lower=1,upper=H> cid[I]; // capsid indices
  int<lower=1,upper=I> sid[N]; // strain indices
  
  int<lower=0> x[N];     // dose in ob
  int<lower=0> y[N];     // dependent variable, number of virus-killed
  int<lower=0> total[N]; // total caterpillars in treatment
  
  vector[I] prior_mu_alpha;
  vector[I] prior_mu_beta;
  
  vector<lower=0>[H] prior_sigma_alpha;
  vector<lower=0>[H] prior_sigma_beta;
}

parameters{
  vector[I] raw_alpha;
  vector[I] raw_beta;
  
  vector<lower=0>[H] sigma_alpha;
  vector<lower=0>[H] sigma_beta;
}

transformed parameters {
  vector[I] alpha;
  vector[I] beta;
  
  vector[N] theta; //predicted proportion virus-killed on logit scale
  vector[N] inv_logit_theta;
  
  for(i in 1:I) {
    alpha[i] = prior_mu_alpha[i] + raw_alpha[i] * sigma_alpha[cid[i]];
    beta[i] = prior_mu_beta[i] + raw_beta[i] * sigma_beta[cid[i]];
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
  
  for (i in 1:I) {
    raw_alpha[i] ~ normal(0,1);
    raw_beta[i] ~ normal(0,1);
  }
  

  //likelihood
  y ~ binomial_logit(total, theta);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_logit_lpmf(y[n] | total[n], theta[n]);
  }
}
