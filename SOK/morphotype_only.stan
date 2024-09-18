data {
  int<lower=1> N; // number of treatments
  int<lower=1> H; // number of capsids
  
  int<lower=1,upper=H> cid[N]; // capsid indices
  
  int<lower=0> y[N]; // speed of kill data
}

parameters{
  vector<lower=0>[H] alpha;
  vector<lower=0>[H] beta;
}

model {
  //priors
  for (h in 1:H) {
    alpha[h] ~ normal(0,10);
    beta[h] ~ normal(0,10);
  }  

  //likelihood
  for (n in 1:N) {
    y[n] ~ gamma(alpha[cid[n]], 1/beta[cid[n]]);
  }
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = gamma_lpdf(y[n] | alpha[cid[n]], 1/beta[cid[n]]);
  }
}
