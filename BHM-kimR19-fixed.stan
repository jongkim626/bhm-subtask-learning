// Stan code for the BHM for subtask learning
// Jong Kim

// The fixed effects model

// data block
data {
  int<lower=1> N; //number of data items
//  int<lower=1, upper=4>J; // number of predictors
  real day[N]; // predictor (covariate)
  real time[N];  // task completion time
//  int<lower=1> J;
//  real time[I];                    //task completion time
//  int<lower=1,upper=4> day[J];  //predictor
//  int<lower=1, upper=14> subtask;
//  int<lower=1, upper=30> pID;
}

// parameters block
parameters {
  vector[2] beta; // intercept and slope
  real<lower=0> sigma_e;
  }

// model block  
model {
  real mu; // local variable declaration that is the mean of time
  for (i in 1:N){
    mu = beta[1] + beta[2] * day[i];
    time[i] ~ lognormal(mu,sigma_e);
  }
}

