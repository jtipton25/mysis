
// regression
data {
  int<lower=1> N;       // number of observations
  int<lower=1> p;       // number of predictors
  matrix[N, p] X;       // predictor
  vector[N] y;          // response
} parameters {
  real<lower=0> s;      // regression standard deviation
  vector[p] beta;       // mean parameters
} model {
  // priors:
  //s ~ cauchy(0, 20);
  // beta ~ multi_normal(0, 100);
  // data model:
  y ~ normal(X * beta, s);
}

