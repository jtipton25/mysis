
// regression
data {
  int<lower=1> N;             // number of observations
  int<lower=1> p;             // number of predictors
  int<lower=1> q;             // number of levels
  matrix[N, p] X;             // predictor
  vector[N] y;                // response
  int idx[N];                 // indicator variable of treatment
} parameters {  
  vector<lower=0>[q] s;       // heterogeneous regression standard deviation
  vector[p] beta;             // mean parameters
} model {
  // priors:
  //s ~ cauchy(0, 20);
  for(i in 2:p){
    beta[i] ~ normal(0, 100);
  }

  // data model:
  for(i in 1:N){
    y[i] ~ normal(X[i] * beta, s[idx[i]]);
  }
}

