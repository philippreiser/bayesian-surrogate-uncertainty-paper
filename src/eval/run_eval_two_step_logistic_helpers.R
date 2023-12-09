library(here)
library(patchwork)
library(xtable)
library(latex2exp)
source(file.path(here(), "src/utils/utils.R"))
source(file.path(here(), "src/data_generation/get_data.R"))
theme_set(theme_bw())

run_mean_ml_mlp_mp <- function(idata, smodel_fit, number_clusters, n_exp, w_exp_gt,
                               poly_degree, M, fct, adapt_delta,
                               file_imodel_ml, file_imodel_mp, file_imodel_mlp,
                               iter_sampling, iter_warmup, chains,
                               use_clustering=TRUE){
  # create data for I-Step using the mean of the T-Posterior
  i_data_mean <- i_step_1d_generator_single(n_exp, poly_degree, M,
                                            smodel_fit, fct, method="mean",
                                            number_draws=1,
                                            upper=1, lower=-1,
                                            fixed_w_exp=w_exp_gt,
                                            fixed_y_exp=idata$y_exp,
                                            sigma_exp=idata$sigma_exp,
                                            prior_only=0,
                                            kmeans.nstart=25, kmeans.iter.max=100)
  # create data for I-Step using the median of the T-Posterior
  i_data_median <- i_step_1d_generator_single(n_exp, poly_degree, M,
                                              smodel_fit, fct, method="median",
                                              number_draws=1,
                                              upper=1, lower=-1,
                                              fixed_w_exp=w_exp_gt,
                                              fixed_y_exp=idata$y_exp,
                                              sigma_exp=idata$sigma_exp,
                                              prior_only=0,
                                              kmeans.nstart=25, kmeans.iter.max=100)
  # create data for I-Step using clusters or all draws of the T-Posterior
  if (use_clustering){
    s_posterior_processing = "single_trial" # use "single_trial" (for clustering)
    iter_sampling_mp <- iter_sampling
    chains_mp <- chains
  } else{
    s_posterior_processing = "single_trial_all"
    smodel_num_chains <- length(smodel_fit$metadata()$id)
    smodel_num_samples <- smodel_fit$metadata()$iter_sampling*smodel_num_chains
    iter_sampling_mp <- iter_sampling*chains/smodel_num_samples
    chains_mp <- 1
  }
  i_data_ml <- i_step_1d_generator_single(n_exp, poly_degree, M,
                                          smodel_fit, fct, method=s_posterior_processing,
                                          number_draws=number_clusters,
                                          upper=1, lower=-1,
                                          fixed_w_exp=w_exp_gt,
                                          fixed_y_exp=idata$y_exp,
                                          sigma_exp=idata$sigma_exp,
                                          prior_only=0,
                                          kmeans.nstart=25, kmeans.iter.max=2000)
  # Run different I-Steps
  # Mean
  i_fit_mean <- get_imodel_fit(i_data_mean$generated, file_imodel_mp, smodel_fit, i_data_mean, method = "mean", adapt_delta = adapt_delta, init = init_w_exp_gt,
                               iter_sampling = iter_sampling, iter_warmup = iter_warmup, chains = chains)
  # Median
  i_fit_median <- get_imodel_fit(i_data_median$generated, file_imodel_mp, smodel_fit, i_data_median, method = "mean", adapt_delta = adapt_delta, init = init_w_exp_gt,
                               iter_sampling = iter_sampling, iter_warmup = iter_warmup, chains = chains)
  # Mixture Likelihood
  i_fit_ml <- get_imodel_fit(i_data_ml$generated, file_imodel_ml, smodel_fit, i_data_ml, method = "single_trial", adapt_delta = adapt_delta, init = init_w_exp_gt,
                             iter_sampling = iter_sampling, iter_warmup = iter_warmup, chains = chains)
  # Mixture Posterior
  i_fit_mp <- get_imodel_fit(i_data_ml$generated, file_imodel_mp, smodel_fit, i_data_ml, method = "multi_draws", adapt_delta = adapt_delta, init = init_w_exp_gt,
                             iter_sampling = iter_sampling_mp, iter_warmup = iter_warmup, chains = chains_mp)
  # Mixture Log Posterior
  i_fit_mlp <- get_imodel_fit(i_data_ml$generated, file_imodel_mlp, smodel_fit, i_data_ml, method = "single_trial", adapt_delta = adapt_delta, init = init_w_exp_gt,
                              iter_sampling = iter_sampling, iter_warmup = iter_warmup, chains = chains)

  if (use_clustering){
    # resample posterior draws according to cluster weights
    posterior_draws_weighted_mp <- sample_draws_weighted_mp(i_fit_mp, i_data_ml$generated, chains, iter_sampling,
                                               variables=c("w_exp[1,1]", "sigma_exp"))
  } else{
    posterior_draws_weighted_mp <- merge_chains(i_fit_mp$post_warmup_draws)
  }

  # Bind w_exp and sigma_exp posterior draws of all methods in a df
  posterior_mean <- merge_chains(rename_variables(i_fit_mean$draws("w_exp"), "mean" = "w_exp[1,1]"))
  posterior_median <- merge_chains(rename_variables(i_fit_median$draws("w_exp"), "median" = "w_exp[1,1]"))
  posterior_ml <- merge_chains(rename_variables(i_fit_ml$draws("w_exp"), "E-lik" = "w_exp[1,1]"))
  posterior_mp <- rename_variables(subset_draws(posterior_draws_weighted_mp, variable = "w_exp[1,1]"), "E-post" = "w_exp[1,1]")
  posterior_mlp <- merge_chains(rename_variables(i_fit_mlp$draws("w_exp"), "E-log-lik" = "w_exp[1,1]"))
  posterior_mean_ml_mp <- bind_draws(posterior_mean, posterior_median, posterior_ml, 
                                     posterior_mlp, posterior_mp, along="variable")
  posterior_df_new <- as_draws_df(posterior_mean_ml_mp)
  posterior_df_new$w_exp_gt <- rep(w_exp_gt, nrow(posterior_df_new))
  posterior_df_new$n_sim <- rep(n_sim, nrow(posterior_df_new))
  posterior_df_new$sigma_sim <- rep(sigma_sim, nrow(posterior_df_new))
  posterior_df_new$n_exp <- rep(n_exp, nrow(posterior_df_new))
  posterior_df_new$sigma_exp <- rep(sigma_exp, nrow(posterior_df_new))
  posterior_df_new$n_clusters <- rep(number_clusters, nrow(posterior_df_new))
  posterior_df_new$poly_degree <- rep(poly_degree, nrow(posterior_df_new))

  sigma_exp_mean <- merge_chains(rename_variables(i_fit_mean$draws("sigma_exp"), "mean" = "sigma_exp"))
  sigma_exp_median <- merge_chains(rename_variables(i_fit_median$draws("sigma_exp"), "median" = "sigma_exp"))
  sigma_exp_ml <- merge_chains(rename_variables(i_fit_ml$draws("sigma_exp"), "E-lik" = "sigma_exp"))
  sigma_exp_mp <- rename_variables(subset_draws(posterior_draws_weighted_mp, variable = "sigma_exp"), "E-post" = "sigma_exp")
  sigma_exp_mlp <- merge_chains(rename_variables(i_fit_mlp$draws("sigma_exp"), "E-log-lik" = "sigma_exp"))
  sigma_exp_mean_ml_mp <- bind_draws(sigma_exp_mean, sigma_exp_median, sigma_exp_ml, sigma_exp_mp, sigma_exp_mlp, along="variable")
  sigma_exp_df_new <- as_draws_df(sigma_exp_mean_ml_mp)
  sigma_exp_df_new$w_exp_gt <- rep(w_exp_gt, nrow(sigma_exp_df_new))
  sigma_exp_df_new$n_sim <- rep(n_sim, nrow(sigma_exp_df_new))
  sigma_exp_df_new$sigma_sim <- rep(sigma_sim, nrow(sigma_exp_df_new))
  sigma_exp_df_new$n_exp <- rep(n_exp, nrow(sigma_exp_df_new))
  sigma_exp_df_new$sigma_exp <- rep(sigma_exp, nrow(sigma_exp_df_new))
  sigma_exp_df_new$n_clusters <- rep(number_clusters, nrow(sigma_exp_df_new))
  sigma_exp_df_new$poly_degree <- rep(poly_degree, nrow(sigma_exp_df_new))

  fit_diagnostics <- tibble(
    "sampling_and_warmup_time" = c(sum(i_fit_mean$metadata()$time$total),
                                   sum(i_fit_ml$metadata()$time$total),
                                   sum(i_fit_mlp$metadata()$time$total),
                                   sum(i_fit_mp$time$chains$total)),
    "sampling_time" = c(sum(i_fit_mean$metadata()$time$sampling),
                        sum(i_fit_ml$metadata()$time$sampling),
                        sum(i_fit_mlp$metadata()$time$sampling),
                        sum(i_fit_mp$time$chains$sampling)),
    "warmup_time" = c(sum(i_fit_mean$metadata()$time$warmup),
                        sum(i_fit_ml$metadata()$time$warmup),
                        sum(i_fit_mlp$metadata()$time$warmup),
                        sum(i_fit_mp$time$chains$warmup)),
    "num_divergent" = c(sum(i_fit_mean$diagnostic_summary()$num_divergent),
                        sum(i_fit_ml$diagnostic_summary()$num_divergent),
                        sum(i_fit_mlp$diagnostic_summary()$num_divergent),
                        sum(subset_draws(i_fit_mp$post_warmup_sampler_diagnostics, "divergent__"))),
    "method" = c("mean", "mix-lik", "mix-log-post", "mix-post")
  )
  fit_diagnostics$w_exp_gt <- w_exp_gt
  fit_diagnostics$n_sim <- n_sim
  fit_diagnostics$sigma_sim <- sigma_sim
  fit_diagnostics$n_exp <- n_exp
  fit_diagnostics$sigma_exp <- sigma_exp
  fit_diagnostics$n_clusters <- number_clusters
  fit_diagnostics$poly_degree <- poly_degree
  return(list("w_posterior" = posterior_df_new,
              "sigma_posterior" = sigma_exp_df_new,
              "fit_diagnostics" = fit_diagnostics))
}
