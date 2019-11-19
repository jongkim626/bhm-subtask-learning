// Author: Jong Kim, jongkim626@gmail.com
// Last Updated: 19 Nov 2019

// Stan code for the BHM for subtask learning
// Varying Intercept

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
  vector[J] u;  // pID intercepts
  vector[K] w;  // subtask intercepts
  real<lower=0> sigma_e;
  real<lower=0> sigma_u;
  real<lower=0> sigma_w;
  }

// model block  
model {
  real mu; // local variable declaration that is the mean of time
  // priors
  u ~ normal(0, sigma_u);
  w ~ normal(0, sigma_w);
  
  //likelihood
  for (i in 1:N){
    mu = beta[1] +u[pID[i]] + w[subtask[i]] + beta[2] * day[i];
    time[i] ~ lognormal(mu,sigma_e);
  }
}

