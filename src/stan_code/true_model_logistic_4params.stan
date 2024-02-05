data {
  int<lower=1> N_sim;           // Number of input/output simulation pairs
  vector[N_sim] y_sim;          // Output simulation variables
  matrix[N_sim, 1] w_sim;       // Input simulation variables

  real sigma_sim_lower;         // Lower bound for sigma_sim
  int prior_only;
}

parameters {
  // surrogate coefficients
  vector[3] c;
  real<lower=0> c_0;
  real<lower=sigma_sim_lower> sigma_sim;         // sigma for simulation values
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  // priors for surrogate coefficients
  lprior += normal_lpdf(c_0 | 2, 1);
  lprior += normal_lpdf(c[1] | 10, 10);
  lprior += normal_lpdf(c[2] | 0, 1);
  lprior += normal_lpdf(c[3] | -1, 1);
  lprior += exponential_lpdf(sigma_sim | 1);
}

model {
  // likelihood including constants
  if (!prior_only) {
    // surrogate likelihood
    target += normal_lpdf(y_sim | c_0 / (1+exp(-c[1] * (w_sim[:, 1] - c[2]))) + c[3], sigma_sim);
  }
  // priors including constants
  target += lprior;
}

generated quantities {
  vector[N_sim] mu_pred_sim;
  mu_pred_sim = c_0 / (1+exp(-c[1] * (w_sim[:, 1] - c[2]))) + c[3];

  array[N_sim] real y_rep_sim;
  y_rep_sim = normal_rng(mu_pred_sim, sigma_sim);

  vector[N_sim] log_lik_sim;
  for (n in 1:N_sim) log_lik_sim[n] = normal_lpdf(y_sim[n] | mu_pred_sim[n], sigma_sim);
}

