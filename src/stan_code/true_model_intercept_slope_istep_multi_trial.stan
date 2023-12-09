data {
  int<lower=1> N_i;          // Number of exp (measured) values
  int<lower=1> N_measures;      // Number of measurements per omega_i
  array[N_i] vector[N_measures] y_i; // Output exp (measured) variables
  
  real sigma_i;
  vector[2] c;                                   // slope
  int prior_only;  // should the likelihood be ignored?
}

parameters {
  matrix[N_i, 1] omega_i;
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  for (n in 1:N_i){
      lprior += normal_lpdf(omega_i[n, 1] | 0, 1);
  }
}

model {
  // Observational model
  if (!prior_only) {
    vector[N_i] mu_y_i =  c[1] + c[2]*omega_i[:, 1];
    for (n in 1:N_i) {
      for (m in 1:N_measures) {
        target += normal_lpdf(y_i[n, m] | mu_y_i[n], sigma_i); // check standard deviation
      }
    }
  }
  // Prior model
  target += lprior;
}

generated quantities {
  matrix[N_i, 1] w_prior;
  w_prior = omega_i;
}
