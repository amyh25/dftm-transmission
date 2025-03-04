data {
  int<lower=0> N; // number of samples
  int<lower=0> K; // n_tree + n_strains + ctrl
  matrix[N, K] X; // data MUST be ordered {intercept, CTRL, tree_sp, MNPV, SNPV}
  vector[N] y; // outcome
  int map[K];
  int hier; // hier == 2 when hierarchical model is applied; else == 0
}
parameters {
  vector[K] beta;
  real<lower=0> sigma; 
  vector[K] gamma; 
  vector<lower=0>[K-hier] tau; 
}
model {
  vector[N] mu; 

  // priors
  // fill the data matrices
  for (k in 1:(K-hier)) {
    tau[k] ~ gamma(2,0.1); //see section 6.9 in STAN user guide
  }
  for (k in 1:K) {
    gamma[k] ~ normal(0,5);
    beta[k] ~ normal(gamma[k], tau[map[k]]);
  }
  sigma ~ gamma(2,0.1);

  // compute the linear predictor
  mu = X * beta;  
  for (n in 1:N) {
    y[n] ~ normal(mu[n], sigma); 
  }
}
generated quantities {
  vector[N] log_lik; 
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | X[n,] * beta, sigma); 
  }
}

