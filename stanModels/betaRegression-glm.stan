
// beta Regression
data {
  int<lower=1> N;                   // sample size
  int<lower=1> p;                   // p predictors
  vector<lower=0,upper=1>[N] y;     // response 
  matrix[N, p] X;                   // predictor matrix
}
parameters {
  vector[p] theta;                  // reg coefficients
  real<lower=0> phi;                // dispersion parameter
}
transformed parameters{
  vector[p] beta;
  beta <- theta * 10;               // matt trick if desired
}
model {
  // model calculations
  vector[N] Xbeta;                  // linear predictor
//  vector<lower=0,upper=1>[N] mu;    // transformed linear predictor
  vector[N] mu;    // transformed linear predictor
  vector[N] A;             // parameter for beta distn
  vector[N] B;             // parameter for beta distn
  Xbeta <- X * beta;
  for (i in 1:N) { 
    mu[i] <- inv_logit(Xbeta[i]);   
  }
  A <- mu * phi;
  B <- (1.0 - mu) * phi;
  // priors
  theta ~ normal(0, 1);   
  phi ~ cauchy(0, 5);               // different options for phi  
  //phi ~ inv_gamma(.001, .001);
  //phi ~ uniform(0, 500);          // put upper on phi if using this
  // likelihood
  y ~ beta(A, B);
}

