get_sbc_results <- function(base_path, n_sims, sigma_sims, date,
                            surrogate_model, i_methods, n_sims_sbc_s,
                            n_sims_sbc_i, use_clustering
                            ){
  if (surrogate_model == "pce"){
    propagate_sigma_t_vec = c(TRUE, FALSE)
  } else{
    propagate_sigma_t_vec = c(FALSE)
  }
  metrics_list <- list()
  eval_sbcs_list <- list()
  results_stats_list <- list()
  for (i in seq_along(sigma_sims)){
    sigma_sim <- sigma_sims[[i]]
    for (j in seq_along(n_sims)){
      n_sim <- n_sims[[j]]
      for (k in seq_along(propagate_sigma_t_vec)){
        propagate_sigma_t <- propagate_sigma_t_vec[[k]]
        plot_dir <- get_plot_dir(base_path, n_sim, propagate_sigma_t, date,
                                 n_sims_sbc_s, n_sims_sbc_i, use_clustering)
        print(plot_dir)
        results_file_name_start <- paste0(plot_dir, "/true_logistic_nsim_"
                                          , n_sim, "_sigmasim_", sigma_sim)
        # results_mean <- readRDS(paste0(results_file_name_start, "_cmean.rds"))
        results_median <- readRDS(paste0(results_file_name_start, "_cmedian.rds"))
        results_ml <- readRDS(paste0(results_file_name_start, "_cml.rds"))
        results_mlp <- readRDS(paste0(results_file_name_start, "_cmlp.rds"))
        results_mp <- readRDS(paste0(results_file_name_start, "_cmp.rds"))
        sbc_results <- list(results_median$result, results_ml$result, results_mlp$result, results_mp$result)
        eval_sbcs <- eval_i_sbcs(sbc_results, i_methods, compute_log_gamma_threshold_history, variables_regex="w_exp")
        eval_sbcs$n_sim <- n_sim
        eval_sbcs$propagate_sigma_t <- propagate_sigma_t
        eval_sbcs$sigma_sim <- sigma_sim
        metrics <- summarize_eval_df(eval_sbcs)
        metrics$n_sim <- n_sim
        metrics$sigma_sim <- sigma_sim
        metrics$propagate_sigma_t <- propagate_sigma_t
        metrics_list[[length(metrics_list)+1]] <- metrics
        eval_sbcs_list[[length(eval_sbcs_list)+1]] <- eval_sbcs
      }
    }
  }
  metrics <- bind_rows(metrics_list)
  if ("mix-log-lik" %in% colnames(metrics)){
    metrics <- metrics %>% 
      rename(
        "mix-log-post" = "mix-log-lik"
      )
  }
  if ("inf_method" %in% colnames(metrics)){
    metrics <- metrics %>% 
      rename(
        "method" = "inf_method"
      )
  }
  metrics
}

eval_i_sbcs <- function(sbc_results, i_methods, compute_log_gamma_threshold_history, variables_regex = NULL){
  n <- length(sbc_results)
  eval_df_list <- vector("list", length = n)
  for (i in c(1:n)){
    # eval_df <- eval_i_sbc(sbc_results[[i]], i_methods[[i]], variables_regex = variables_regex)
    log_gamma_th <- compute_log_gamma_threshold_history(sbc_results[[i]], variables_regex = variables_regex)
    eval_df <- data.frame(inf_method=i_methods[[i]], gamma_sbc = log_gamma_th)
    # eval_df$gamma_sbc <- log_gamma_th
    eval_df_list[[i]] <- eval_df
    
  }
  eval_dfs <- do.call(rbind, eval_df_list)
  eval_dfs
}

summarize_eval_df <- function(eval_df){
  metrics <- eval_df %>%
    group_by(inf_method) %>%
    summarise(
      gamma_sbc = mean(gamma_sbc)
    )
  metrics
}
