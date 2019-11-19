// Stan code for the BHM for subtask learning
// Varying Intercept and Varying slopes

// Data Block

data {
  int<lower=1> N; //number of data items
  real time[N];
  real<lower=1, upper=4> day[N];
  int<lower=1> J; // number of subjects
  int<lower=1> K; // number of subtasks
  int<lower=1, upper=J> pID[N];
  int<lower=1, upper=K> subtask[N];
  //real day[N]; // predictor (covariate)
}

// parameters block
parameters {
  vector[2] beta; // intercept and slope
  real<lower=0> sigma_e;
  matrix[2,J] u;
  vector<lower=0>[2] sigma_u;
  matrix[2,K] w;
  vector<lower=0>[2] sigma_w;
  }

// model block  
model {
  real mu; // local variable declaration that is the mean of time
  // priors
  u[1] ~ normal(0, sigma_u[1]);
  u[2] ~ normal(0, sigma_u[2]);
  w[1] ~ normal(0, sigma_w[1]);
  w[2] ~ normal(0, sigma_w[2]);
  
  //likelihood
  for (i in 1:N){
    mu = beta[1] +u[1, pID[i]] + w[1, subtask[i]]
        + (beta[2] + u[2, pID[i]] + w[2, subtask[i]]) * day[i];
    time[i] ~ lognormal(mu,sigma_e);
  }
}

