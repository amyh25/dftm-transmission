data {
  int<lower=1> N; // number of treatments
  
  int<lower=0> y[N]; // speed of kill data
}

parameters{
  real<lower=0> alpha;
  real<lower=0> beta;
}

model {
  //priors
  alpha ~ normal(0,10);
  beta ~ normal(0,10);
  

  //likelihood
  for (n in 1:N) {
    y[n] ~ gamma(alpha, 1/beta);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = gamma_lpdf(y[n] | alpha, 1/beta);
  }
}
