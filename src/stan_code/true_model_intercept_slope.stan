data {
  int<lower=1> N_sim;           // Number of input/output simulation pairs
  vector[N_sim] y_t;          // Output simulation variables
  matrix[N_sim, 1] omega_t;       // Input simulation variables

  real sigma_a;
  int prior_only;
}

parameters {
  vector[2] c;                                   // intercept & slope
  // real<lower=sigma_a_lower> sigma_a;         // sigma for simulation values
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += normal_lpdf(c[1] | 0, 10);
  lprior += normal_lpdf(c[2] | 0, 10);
  // lprior += exponential_lpdf(sigma_a | 1);
}

model {
  // likelihood including constants
  if (!prior_only) {
    // initialize slope predictor term
    target += normal_lpdf(y_t | c[1] + c[2]*omega_t[:, 1], sigma_a);
  }
  // priors including constants
  target += lprior;
}

generated quantities {
  vector[N_sim] mu_pred_sim;
  mu_pred_sim = c[1] + c[2]*omega_t[:, 1];

  array[N_sim] real y_rep_sim;
  y_rep_sim = normal_rng(mu_pred_sim, sigma_a);

  vector[N_sim] log_lik_sim;
  for (n in 1:N_sim) log_lik_sim[n] = normal_lpdf(y_t[n] | mu_pred_sim[n], sigma_a);
}

