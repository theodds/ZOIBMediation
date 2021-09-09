data{
  int<lower=1> n;
  int<lower=1> np_y; // number of columns for the design matrix of outcome
  int<lower=1> np_m; // number of columns for design matrix of mediator
  vector<lower=0, upper=1>[n] y;
  vector<lower=0, upper=1>[n] m;
  matrix[n, np_y] X_y;
  matrix[n, np_m] X_m;
}

parameters {
  matrix[np_m, 4] beta_mediator;
  matrix[np_y, 4] beta_outcome;
}

transformed parameters {
  matrix[n,4] eta_y;
  matrix[n,4] eta_m;

  eta_y = X_y * beta_outcome;
  eta_m = X_m * beta_mediator;

  eta_y[1:n,3] = inv_logit(eta_y[1:n,3]);
  eta_y[1:n,4] = exp(eta_y[1:n,4]);
  eta_m[1:n,3] = inv_logit(eta_m[1:n,3]);
  eta_m[1:n,4] = exp(eta_m[1:n,4]);

}

model {

  // Coefficient Model
  for(j in 1:4) {
    for(i in 1:np_y) {
      beta_outcome[i,j] ~ normal(0.0,10.0);
    }
    for(i in 1:np_m) {
      beta_mediator[i,j] ~ normal(0.0,10.0);
    }
  }

  // zero one inflated beta likelihood
  for (i in 1:n) {
    if (y[i] == 0.0) {
      target += bernoulli_logit_lpmf(1 | eta_y[i,1]) ;
    } else if (y[i] == 1) {
      target += bernoulli_logit_lpmf(1 | eta_y[i,2]) + bernoulli_logit_lpmf(0 | eta_y[i,1]);
    } else {
      target += bernoulli_logit_lpmf(0 | eta_y[i,2]) + bernoulli_logit_lpmf(0 | eta_y[i,1]) +
                  beta_lpdf(y[i] | eta_y[i,3] * eta_y[i,4], (1 - eta_y[i,3]) * eta_y[i,4]);
    }

    if(m[i] == 0.0) {
      target +=bernoulli_logit_lpmf(1 | eta_m[i,1]);
    }
    else if(m[i] == 1.0) {
      target += bernoulli_logit_lpmf(1 | eta_m[i,2]) + bernoulli_logit_lpmf(0 | eta_m[i,1]);
    }
    else {
      target += bernoulli_logit_lpmf(0 | eta_m[i,2]) + bernoulli_logit_lpmf(0 | eta_m[i,1]) +
                  beta_lpdf(m[i] | eta_m[i,3] * eta_m[i,4], (1 - eta_m[i,3]) * eta_m[i,4]);
    }
  }

}



