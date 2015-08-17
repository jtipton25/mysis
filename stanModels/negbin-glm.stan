
// negative binomial parameterized as eta (log(mu)) and dispersion (phi)
// see p286 in stan-reference-2.4.0.pdf
// a basic GLM example
data {
  int<lower=1> N;    // rows of data
  int<lower=0> p;
  matrix[N, p] X;       // predictor
  int<lower=0> y[N]; // response
} transformed data {
  vector[p] J;
  matrix[p, p] I_p;
  // calculate placeholder values
  J <- rep_vector(0, p);
  I_p <- diag_matrix(rep_vector(1, p));
}
parameters {
  real<lower=0> phi; // neg. binomial dispersion parameter
  vector[p] beta;  // neg. binomial mean parameters
}
model {
  // priors:
  phi ~ cauchy(0, 20);
  beta ~ multi_normal(J, 100 * I_p);
  // data model:
  y ~ neg_binomial_2_log(X * beta, phi);
}

