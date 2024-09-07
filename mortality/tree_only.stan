data {
  int<lower=1> N; // number of treatments
  int<lower=1> I; // number of strains
  int<lower=1> J; // number of tree species
  
  int<lower=1,upper=I> sid[N]; // strain indices
  int<lower=1,upper=J> tid[N]; // tree indices
  
  int<lower=0> x[N];     // dose in ob
  int<lower=0> y[N];     // dependent variable, number of virus-killed
  int<lower=0> total[N]; // total caterpillars in treatment
}

parameters{
  matrix[I,J] raw_alpha;
  matrix[I,J] raw_beta;
  
  row_vector<lower=-10,upper=10>[J] mu_alpha;
  row_vector<lower=0,upper=.01>[J] mu_beta;
  
  real<lower=0,upper=10> sigma_alpha;
  real<lower=0,upper=.01> sigma_beta;
}

transformed parameters {
  matrix[I,J] alpha;
  matrix[I,J] beta;
  
  vector[N] theta; // predicted proportion virus-killed on logit scale
  vector[N] inv_logit_theta;
  
  for(i in 1:I) {
    alpha[i] = mu_alpha + raw_alpha[i] * sigma_alpha;
    beta[i] = fmax(mu_beta + raw_beta[i] * sigma_beta,0);
  }

  for(n in 1:N) {
    theta[n] = alpha[sid[n],tid[n]] + beta[sid[n],tid[n]] * x[n];
    
    if (theta[n] > 700) { // to avoid floating point errors with logit evaluating to inf or -inf
      theta[n] = 700;
    } else if (theta[n] < -700) {
      theta[n] = -700;
    }
  }
  
  inv_logit_theta = inv_logit(theta);
}

model {
  //priors
  sigma_alpha ~ normal(0,1);
  sigma_beta ~ normal(0,.001);
  
  for (i in 1:I) {
    raw_alpha[i] ~ normal(0,1);
    raw_beta[i] ~ normal(0,1);
  }
  
  mu_alpha ~ normal(0,1);
  mu_beta ~ normal(0,1);
  

  // likelihood
  y ~ binomial_logit(total, theta);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_logit_lpmf(y[n] | total[n], theta[n]);
  }
}
