
// Negative Binomial GLM
// see Gelman and Hill (2007)
data {
  int<lower=1> N;       // number of observations
  int<lower=0> p;
  matrix[N, p] X;       // predictor
  int<lower=0> y[N];    // response
} transformed data {
  vector[p] J;
  matrix[p, p] I_p;
  // calculate placeholder values
  J <- rep_vector(0, p);
  I_p <- diag_matrix(rep_vector(1, p));
}
parameters {
  real<lower=0> phi;   // Negative Binomial dispersion parameter
  vector[p] beta;      // Negative Binomial mean parameters
} model {
  // priors:
  phi ~ cauchy(0, 20);
  beta ~ multi_normal(J, 100 * I_p);
  // data model:
  y ~ neg_binomial_2_log(X * beta, phi);
}

