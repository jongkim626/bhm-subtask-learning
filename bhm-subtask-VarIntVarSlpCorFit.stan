// Author: Jong Kim, jongkim626@gmail.com
// Last Updated: 19 Nov 2019

// Stan code for the BHM for subtask learning
// Varying Intercept and Varying slopes with Consideration of Correlation

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

// Parameters Block
parameters {
  vector[2] beta; // intercept and slope
  real<lower=0> sigma_e;
  
  matrix[2,J] z_u;
  vector<lower=0>[2] sigma_u;
  cholesky_factor_corr[2] L_u;
  
  matrix[2,K] z_w;
  vector<lower=0>[2] sigma_w;
  cholesky_factor_corr[2] L_w;
  }
  
transformed parameters{
  matrix[2,J] u;
  matrix[2,K] w;
  
  u = diag_pre_multiply(sigma_u, L_u) * z_u;    // Subject varying effects
  w = diag_pre_multiply(sigma_w, L_w) * z_w;    // Subtask varying effects
  
}



// model block  
model {
  real mu; // local variable declaration that is the mean of time

    // priors
  L_u ~ lkj_corr_cholesky(2.0);
  L_w ~ lkj_corr_cholesky(2.0);
  to_vector(z_u) ~ normal(0,1);
  to_vector(z_w) ~ normal(0,1);

  
  //likelihood
  for (i in 1:N){
    mu = beta[1] +u[1, pID[i]] + w[1, subtask[i]]
        + (beta[2] + u[2, pID[i]] + w[2, subtask[i]]) * day[i];
    time[i] ~ lognormal(mu,sigma_e);
  }
}
