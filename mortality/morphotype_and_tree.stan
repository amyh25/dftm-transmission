data {
  int<lower=1> N; // number of treatments
  int<lower=1> H; // number of capsids
  int<lower=1> I; // number of strains
  int<lower=1> J; // number of tree species
  
  int<lower=1,upper=H> cid[I]; // capsid indices
  int<lower=1,upper=I> sid[N]; // strain indices
  int<lower=1,upper=J> tid[N]; // tree indices
  
  int<lower=0> x[N];     // dose in ob
  int<lower=0> y[N];     // dependent variable, number of virus-killed
  int<lower=0> total[N]; // total caterpillars in treatment
}

parameters{
  matrix[I,J] raw_alpha;
  matrix[I,J] raw_beta;
  
  matrix<lower=-10,upper=10>[H,J] mu_alpha;
  matrix<lower=0,upper=.01>[H,J] mu_beta;
  
  vector<lower=0,upper=10>[H] sigma_alpha;
  vector<lower=0,upper=.01>[H] sigma_beta;
}

transformed parameters {
  matrix[I,J] alpha;
  matrix[I,J] beta;
  
  vector[N] theta; // predicted proportion virus-killed on logit scale
  vector[N] inv_logit_theta;
  
  for (i in 1:I) {
    for (j in 1:J) {
      alpha[i,j] = mu_alpha[cid[i],j] + raw_alpha[i,j] * sigma_alpha[cid[i]];
      beta[i,j] = fmax(mu_beta[cid[i],j] + raw_beta[i,j] * sigma_beta[cid[i]],0);
    }
  }

  for(n in 1:N) {
    theta[n] = alpha[sid[n],tid[n]] + beta[sid[n],tid[n]] * x[n];
    
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
    for (j in 1:J) {
      raw_alpha[i,j] ~ normal(0,1);
      raw_beta[i,j] ~ normal(0,1);
    }
  }
  
  for (h in 1:H) {
    sigma_alpha[h] ~ normal(0,1);
    sigma_beta[h] ~ normal(0,.001);
    
    for (j in 1:J) {
      mu_alpha[h,j] ~ normal(0,1);
      mu_beta[h,j] ~ normal(0,1);
    }
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
