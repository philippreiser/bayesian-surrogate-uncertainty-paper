get_dirs <- function(){
  stan_code_dir <- file.path(here(),"src/stan_code")
  plot_dir <- file.path(here(), "plots")
  fitted_models_dir <- file.path(here(), "fitted_models")
  c(stan_code_dir, plot_dir, fitted_models_dir)
}

rmse_draws <- function(theta_gt, theta) {
  sqrt(mean((theta_gt - theta)^2))
}

mae_draws <- function(theta_gt, theta) {
  mean(abs(theta_gt - theta))
}

bias_draws <- function(theta_gt, theta) {
  mean(theta) - theta_gt
}

sharpness_draws <- function(x, q) {
  quantile(x, probs=c(1-((1-q)/2))) - quantile(x, probs=c((1-q)/2))
}

sharpness_s50s66s90s95s99 <- function(x) {
  list(
    s50 = sharpness_draws(x, 0.5),
    s66 = sharpness_draws(x, 0.66),
    s90 = sharpness_draws(x, 0.90),
    s95 = sharpness_draws(x, 0.95),
    s99 = sharpness_draws(x, 0.99)
  )
}

w_gt_in_ci_draws <- function(x, w_gt, q) {
  q_upper <- quantile(x, probs=c(1-((1-q)/2)))
  q_lower <- quantile(x, probs=c((1-q)/2))
  1*((q_lower <= w_gt) && (w_gt <= q_upper))
}

w_gt_in_ci50ci66ci90ci95ci99 <- function(x, w_gt) {
  list(
    ci50 = w_gt_in_ci_draws(x, w_gt, 0.5),
    ci66 = w_gt_in_ci_draws(x, w_gt, 0.66),
    ci90 = w_gt_in_ci_draws(x, w_gt, 0.90),
    ci95 = w_gt_in_ci_draws(x, w_gt, 0.95),
    ci99 = w_gt_in_ci_draws(x, w_gt, 0.99)
  )
}

# modified code from
# https://github.com/martinmodrak/sbc_test_quantities_paper
compute_log_gamma_history <- function(ranks, max_rank) {
  rank_t <- rep(0, max_rank + 1)
  log_gamma <- rep(NA_real_, length(ranks))
  dummy <- rep(NA_real_, length(ranks))
  
  K <- max_rank + 1
  z <- (1:(K - 1)) / K
  
  for(i in 1:length(ranks)) {
    rank_t[ranks[i] + 1] <- rank_t[ranks[i] + 1] + 1
    scaled_ecdf <- cumsum(rank_t[1:max_rank])
    
    log_gamma[i] <- log(2) + min(
      pbinom(scaled_ecdf, i, z, log = TRUE),
      pbinom(scaled_ecdf - 1, i, z, lower.tail = FALSE, log = TRUE)
    )
  }
  log_gamma
}

# modified code from
# https://github.com/martinmodrak/sbc_test_quantities_paper
compute_log_gamma_threshold_history <- function(results, min_sim_id = 0, max_sim_id = Inf, wrap_cols = 4, variables_regex = NULL, ylim = NULL) {
  unique_max_rank <- unique(results$stats$max_rank)
  if(length(unique_max_rank) > 1) {
    stop("Requires all max_rank to be equal")
  }
  stats <- results$stats
  if(!is.null(variables_regex)) {
    stats <- stats %>% filter(grepl(variables_regex, variable))
  }
  max_sim_id_to_show <- min(max_sim_id, max(stats$sim_id))
  log_gamma_threshold <- log(SBC:::adjust_gamma(N = max(stats$sim_id), L = 1, K = unique_max_rank + 1))
  log_gamma <- compute_log_gamma_history(stats$rank, unique(stats$max_rank))[nrow(stats)]
  log_gamma_gamma_th <- log_gamma - log_gamma_threshold

  log_gamma_gamma_th
}

get_cluster_draws <- function(c_0, c, sigma_sim=NULL, number_centers=25, nstart=25, iter.max=100) {
  smodel_coeffs <- cbind(c_0, c, sigma_sim)
  cl <- kmeans(smodel_coeffs, centers = number_centers, nstart = nstart, 
               iter.max = iter.max)
  kmeans_centroids <- cl$centers
  c_0 <-  kmeans_centroids[, 1]
  c <- kmeans_centroids[, 2:(ncol(c)+1)]
  sigma_sim <- kmeans_centroids[, ncol(c)+2]
  cluster_weights <- cl$size/(sum(cl$size))
  return(list("c_0" = c_0,
              "c" = c,
              "sigma_sim" = sigma_sim,
              "cluster_weights" = cluster_weights))
}

init_w_exp_gt <- function() list(w_exp = as.matrix(w_exp_gt))

get_experiment_dir_name <- function(head_path, experiment_name, model,
                                    experiment_details=NULL, date_time=NULL){
  if (is.null(date_time)){
    date_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  }
  if (is.null(experiment_details)){
    dir_name <- file.path(head_path, experiment_name, model, date_time)
  } else{
    dir_name <- file.path(head_path, experiment_name, model,
                          paste0(date_time, experiment_details))
  }
  return(dir_name)
}

get_plot_dir <- function(base_path, n_sim, propagate_sigma_t, date,
                         n_sims_sbc_s, n_sims_sbc_i, use_clustering){
  file_list <- list.files(base_path)
  head_path <- file_list[grepl(paste0("^", date, "_\\d{2}-\\d{2}-\\d{2}_nsim_", n_sim, 
                                      "+_sigmasim_0\\.01_ncluster_25_useclustering_", 
                                      use_clustering, "_polydegree_5sbc_s",
                                      n_sims_sbc_s,"sbc_i",n_sims_sbc_i,"prop_sigt", propagate_sigma_t, "$"),
                  file_list)]
  file.path(base_path, head_path)
}

get_plot_dir_old <- function(head_path, n_sim, sigma_sim, n_sims_sbc_s,
                             n_sims_sbc_i, propagate_sigma_t){
  paste0(head_path, "nsim", n_sim, "sigmasim_", sigma_sim, "sbc_s", n_sims_sbc_s, "sbc_i", n_sims_sbc_i,
         "prop_sigt", propagate_sigma_t)
}

get_stan_files <- function(stan_code_dir, model){
  i_model_file_point <- file.path(stan_code_dir,
                                  paste0("true_model_", model,
                                         "_istep_multi_trial.stan"))
  i_model_file_e_lik <- file.path(stan_code_dir,
                                  paste0("true_model_", model,
                                         "_istep_single_trial.stan"))
  i_model_file_e_log_lik <- file.path(stan_code_dir,
                                      paste0("true_model_", model,
                                             "_istep_llm.stan"))
  t_model_file <- file.path(stan_code_dir,
                            paste0("true_model_", model, ".stan"))
  return(data.frame("i_model_point"=i_model_file_point,
                    "i_model_e_lik"=i_model_file_e_lik,
                    "i_model_e_log_lik"=i_model_file_e_log_lik,
                    "t_model"=t_model_file))
}


#' Return the stan code file names in a directory
#' smodel = model for T-Step
#' imodel = model for I-Step
#' multi_trial = E-Post
#' single_trial = E-Lik
#' mlp = E-Log-Lik
#' _sigma = propagate sigma_A
#' @param stan_code_dir 
#' @param surrogate_model 
#' @param propagate_sigma_t 
#'
#' @return
#' @export
#'
#' @examples
get_stan_models_logistic <- function(stan_code_dir, surrogate_model,
                                     propagate_sigma_t){
  if (surrogate_model == "true_model"){
    smodel_file <- file.path(stan_code_dir, "true_model_logistic_4params.stan")
    file_imodel_mp <- file.path(stan_code_dir, "true_model_logistic_4params_istep_multi_trial.stan")
    file_imodel_ml <- file.path(stan_code_dir, "true_model_logistic_4params_istep_single_trial.stan")
    file_imodel_mlp <- file.path(stan_code_dir, "true_model_logistic_4params_istep_mlp.stan")
    # filenames of stan file that propagate sigma from t-step
  } else if (surrogate_model == "pce"){
    smodel_file <- file.path(stan_code_dir, "multi_legendre_pce.stan")
    if (propagate_sigma_t){
      file_imodel_mp <- file.path(stan_code_dir, "multi_legendre_pce_istep_mp_sigma.stan")
      file_imodel_ml <- file.path(stan_code_dir, "multi_legendre_pce_istep_ml_sigma.stan")
      file_imodel_mlp <- file.path(stan_code_dir, "multi_legendre_pce_istep_mlp_sigma.stan")
    } else{
      file_imodel_mp <- file.path(stan_code_dir, "multi_legendre_pce_istep_mp.stan")
      file_imodel_ml <- file.path(stan_code_dir, "multi_legendre_pce_istep_ml.stan")
      file_imodel_mlp <- file.path(stan_code_dir, "multi_legendre_pce_istep_mlp.stan")
    }
  }
  return(data.frame("imodel_mp"=file_imodel_mp,
                    "imodel_ml"=file_imodel_ml,
                    "imodel_mlp"=file_imodel_mlp,
                    "smodel"=smodel_file))
}
