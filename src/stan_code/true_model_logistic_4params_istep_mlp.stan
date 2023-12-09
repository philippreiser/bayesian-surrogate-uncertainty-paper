data {
  int<lower=1> N_exp;          // Number of exp (measured) values
  int<lower=1> N_measures;      // Number of measurements per w_exp
  int<lower=1> N_clusters;
  array[N_exp] vector[N_measures] y_exp; // Output exp (measured) variables

  vector[N_clusters] c_0;                                      // coefficient for constant "polynomial"
  matrix[N_clusters, 3] c;                                   // coefficients of non-constant polynomials
  vector[N_clusters] cluster_weights;
  int prior_only;  // should the likelihood be ignored?
}

parameters {
  real<lower=0.005, upper=0.05> sigma_exp;
  matrix<lower=-1, upper=1>[N_exp, 1] w_exp;
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += uniform_lpdf(sigma_exp | 0.005, 0.05); //exponential_lpdf(sigma_exp | 30);
  for (n in 1:N_exp){
      lprior += normal_lpdf(w_exp[n, 1] | 0, 0.5);
  }
}

model {
  // Observational model
  if (!prior_only) {
    for (k in 1:N_clusters) {
      vector[N_exp] mu_y_exp =  c_0[k] / (1+exp(-c[k, 1] * (w_exp[:, 1] - c[k, 2]))) + c[k, 3];
      for (n in 1:N_exp) {
        target += cluster_weights[k] * normal_lpdf(to_vector(y_exp[n, :]) | mu_y_exp[n], sigma_exp);
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
