data {
  int<lower=1> N; // number of treatments
  int<lower=1> J; // number of tree species
  
  int<lower=1,upper=J> tid[N]; // tree indices
  
  int<lower=0> y[N]; // speed of kill data
}

parameters{
  vector<lower=0>[J] alpha;
  vector<lower=0>[J] beta;
}

model {
  //priors
  for (j in 1:J) {
    alpha[j] ~ normal(0,10);
    beta[j] ~ normal(0,10);
  }
  

  //likelihood
  for (n in 1:N) {
    y[n] ~ gamma(alpha[tid[n]], 1/beta[tid[n]]);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = gamma_lpdf(y[n] | alpha[tid[n]], 1/beta[tid[n]]);
  }
}
