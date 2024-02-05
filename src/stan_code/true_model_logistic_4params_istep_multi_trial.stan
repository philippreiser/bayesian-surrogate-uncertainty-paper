data {
  int<lower=1> N_exp;          // Number of exp (measured) values
  int<lower=1> N_measures;      // Number of measurements per w_exp
  int<lower=1> N_clusters;
  array[N_exp] vector[N_measures] y_exp; // Output exp (measured) variables

  real c_0;                                      // surrogate coefficient
  vector[3] c;                                   // surrogate coefficients
  int prior_only;  // should the likelihood be ignored?
}

parameters {
  real<lower=0.005, upper=0.05> sigma_exp;
  matrix<lower=-1, upper=1>[N_exp, 1] w_exp;
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += uniform_lpdf(sigma_exp | 0.005, 0.05); // prior on sigma
  for (n in 1:N_exp){
      lprior += normal_lpdf(w_exp[n, 1] | 0, 0.5); // prior on input variables
  }
}

model {
  // Observational model
  if (!prior_only) {
    vector[N_exp] mu_y_exp =  c_0 / (1+exp(-c[1] * (w_exp[:, 1] - c[2]))) + c[3];
    for (n in 1:N_exp) {
      target += normal_lpdf(to_vector(y_exp[n, :]) | mu_y_exp[n], sigma_exp);
    }
  }
  // Prior model
  target += lprior;
}

generated quantities {
  matrix<lower=-1, upper=1>[N_exp, 1] w_prior;
  w_prior = w_exp;
}
