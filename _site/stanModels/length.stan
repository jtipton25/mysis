
// regression
data {
  int<lower=1> N;             // number of observations
  int<lower=1> p;             // number of predictors
  matrix[N, p] X;             // predictor
  vector[N] y;                // response
} parameters {  
  real<lower=0> s;            // regression standard deviation
  vector[p] beta;             // mean parameters
  vector<lower=0>[p] gamma2;  // LASSO mixing parameters
  real<lower=0> lambda2;      // LASSO shrinkage parameters
} model {
  // priors:
  //s ~ cauchy(0, 20);
  lambda2 ~ gamma(1, 1);      // Prior for weak LASSO shrinkage 
  gamma2 ~ exponential(lambda2 / 2.0);
  for(i in 2:p){
    beta[i] ~ normal(0, s * sqrt(gamma2[i]));
  }

  // data model:
  for(i in 1:N){
    y[i] ~ normal(X[i] * beta, s[idx[i]]);
  }
}

