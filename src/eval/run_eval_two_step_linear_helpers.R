get_i_eval_method_df <- function(omega_t, y_t, sigma_a, c_gt, prior_c, omega_i_gt, 
                                 y_i, omega_i_range, sigma_i, posterior_omega_i, method)
{
  eval_i_step_df <- data.frame("omega_t" = rep(omega_t, length(omega_i_range)), 
                               "y_t" = rep(y_t, length(omega_i_range)),
                               "sigma_a" = rep(sigma_a, length(omega_i_range)),
                               "c_gt" = rep(c_gt, length(omega_i_range)),
                               "prior_c" = rep(prior_c, length(omega_i_range)),
                               "omega_i_gt" = rep(omega_i_gt, length(omega_i_range)),
                               "y_i" = rep(y_i, length(omega_i_range)),
                               "omega_i_range" = omega_i_range,
                               "sigma_i" = rep(sigma_i, length(omega_i_range)),
                               "posterior_omega_i" = rep(posterior_omega_i, length(omega_i_range)),
                               "method" = rep(method, length(omega_i_range)))
}

get_i_eval_df <- function(omega_t, y_t, sigma_a, c_gt, prior_c, omega_i_gt,
                          y_i, omega_i_range, sigma_i, posterior_point, posterior_st, posterior_mt){
  eval_i_df_point <- get_i_eval_method_df(omega_t, y_t, sigma_a, c_gt, prior_c, omega_i_gt,
                                         y_i, omega_i_range, sigma_i, posterior_point, "point")
  eval_i_df_st <- get_i_eval_method_df(omega_t, y_t, sigma_a, c_gt, prior_c, omega_i_gt,
                                       y_i, omega_i_range, sigma_i, posterior_st, "mixture_likelihood")
  eval_i_df_mt <- get_i_eval_method_df(omega_t, y_t, sigma_a, c_gt, prior_c, omega_i_gt,
                                       y_i, omega_i_range, sigma_i, posterior_mt, "mixture_posterior")
  rbind(eval_i_df_point, eval_i_df_st, eval_i_df_mt)
}

get_i_eval_df_non_tidy <- function(omega_i_range, posterior_densities_point, posterior_densities_e_lik, posterior_densities_e_post, posterior_densities_e_log_lik, sigma_a, omega_i_gt){
  data.frame("omega_i" = omega_i_range,
             "point" = posterior_densities_point,
             "e_lik" = posterior_densities_e_lik,
             "e_post" = posterior_densities_e_post,
             "e_log_lik" = posterior_densities_e_log_lik,
             "sigma_a" = sigma_a,
             "omega_i_gt" = omega_i_gt
             )
}
