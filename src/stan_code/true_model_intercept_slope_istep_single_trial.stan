data {
  int<lower=1> N_i;          // Number of exp (measured) values
  int<lower=1> N_measures;      // Number of measurements per omega_i
  int<lower=1> N_clusters;
  array[N_i] vector[N_measures] y_i; // Output exp (measured) variables

  real sigma_i;
  matrix[N_clusters, 2] c;                                   // coefficients of non-constant polynomials
  vector[N_clusters] cluster_weights;
  int prior_only;  // should the likelihood be ignored?
}

parameters {
  // real<lower=0> sigma_i;
  matrix[N_i, 1] omega_i;
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  // lprior += exponential_lpdf(sigma_i | 1);
  for (n in 1:N_i){
      lprior += normal_lpdf(omega_i[n, 1] | 0, 1);
  }
}

model {
  // Observational model
  if (!prior_only) {
    real llik = 0;
    llik += -log(N_clusters);
    vector[N_clusters] mu_y_i = c[:, 1] + c[:, 2]*omega_i[1, 1];
    vector[N_clusters] cluster_lliks;
    for (k in 1:N_clusters) {
      cluster_lliks[k] = normal_lpdf(y_i[1, 1] | mu_y_i[k], sigma_i);
    }
    
    llik += log_sum_exp(cluster_lliks);
    target += llik;
  }
  // Prior model
  target += lprior;
}

generated quantities {
  matrix[N_i, 1] w_prior;
  w_prior = omega_i;
}
