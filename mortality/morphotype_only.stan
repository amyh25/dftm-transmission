data {
  int<lower=1> N; // number of treatments
  int<lower=1> H; // number of capsids
  int<lower=1> I; // number of strains
  
  int<lower=1,upper=H> cid[I]; // capsid indices
  int<lower=1,upper=I> sid[N]; // strain indices
  
  int<lower=0> x[N];     // dose in ob
  int<lower=0> y[N];     // dependent variable, number of virus-killed
  int<lower=0> total[N]; // total caterpillars in treatment
}

parameters{
  vector[I] raw_alpha;
  vector[I] raw_beta;
  
  vector<lower=-10,upper=10>[H] mu_alpha;
  vector<lower=0,upper=.01>[H] mu_beta;
  
  vector<lower=0,upper=10>[H] sigma_alpha;
  vector<lower=0,upper=.01>[H] sigma_beta;
}

transformed parameters {
  vector[I] alpha;
  vector[I] beta;
  
  vector[N] theta; // predicted proportion virus-killed on logit scale
  vector[N] inv_logit_theta;
  
  for(i in 1:I) {
    alpha[i] = mu_alpha[cid[i]] + raw_alpha[i] * sigma_alpha[cid[i]];
    beta[i] = fmax(mu_beta[cid[i]] + raw_beta[i] * sigma_beta[cid[i]],0);
  }

  for(n in 1:N) {
    theta[n] = alpha[sid[n]] + beta[sid[n]] * x[n];
    
    if (theta[n] > 700) { // to avoid floating point errors with logit evaluating to inf or -inf
      theta[n] = 700;
    } else if (theta[n] < -700) {
      theta[n] = -700;
    }
    
    inv_logit_theta[n] = inv_logit(theta[n]);
  }
}

model {
  // priors
  for (i in 1:I) {
    raw_alpha[i] ~ normal(0,1);
    raw_beta[i] ~ normal(0,1);
  }
  
  for (h in 1:H) {
    sigma_alpha[h] ~ normal(0,1);
    sigma_beta[h] ~ normal(0,.001);
  
    mu_alpha[h] ~ normal(0,1);
    mu_beta[h] ~ normal(0,1);
  }
  

  // likelihood
  for (n in 1:N) {
    y[n] ~ binomial_logit(total[n], theta[n]);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_logit_lpmf(y[n] | total[n], theta[n]);
  }
}
