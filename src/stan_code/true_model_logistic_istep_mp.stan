data {
  int<lower=1> N_exp;          // Number of exp (measured) values
  int<lower=1> N_measures;      // Number of measurements per w_exp
  int<lower=1> N_clusters;
  array[N_exp] vector[N_measures] y_exp; // Output exp (measured) variables

  real c_0;                                      // coefficient for constant "polynomial"
  vector[2] c;                                   // coefficients of non-constant polynomials
  int prior_only;  // should the likelihood be ignored?
}

parameters {
  // real<lower=0> sigma_exp;
  matrix<lower=-1, upper=1>[N_exp, 1] w_exp;
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  // lprior += exponential_lpdf(sigma_exp | 1);
  for (n in 1:N_exp){
      lprior += normal_lpdf(w_exp[n, 1] | 0, 0.5);
  }
}

model {
  // Observational model
  if (!prior_only) {
    vector[N_exp] mu_y_exp =  c[1] / (1+exp(-c[2] * (w_exp[:, 1] - c_0))) - 1;
    for (n in 1:N_exp) {
      for (m in 1:N_measures) {
        target += normal_lpdf(y_exp[n, m] | mu_y_exp[n], 0.01);
      }
    }
  }
  // Prior model
  target += lprior;
}

generated quantities {
  matrix<lower=-1, upper=1>[N_exp, 1] w_prior;
  w_prior = w_exp;
}
