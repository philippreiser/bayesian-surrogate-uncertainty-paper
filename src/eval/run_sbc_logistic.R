library(SBC)
library(patchwork)
library(here)
library(xtable)
source(file.path(here(), "src/utils/utils.R"))
source(file.path(here(), "src/eval/sbc_logistic_backend.R"))
source(file.path(here(), "src/utils/true_models_helpers.R"))
source(file.path(here(), "src/data_generation/get_data.R"))
source(file.path(here(), "src/models/fit_2-step.R"))

theme_set(theme_bw())

# Setup
cache_dir <- file.path(here(),"evaluation/SBC/_basic_usage_SBC_cache")
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}

args = commandArgs(trailingOnly=TRUE)
config_dir <- file.path(here(), "config/logistic_case_study")
if (length(args)==0){
  config_file_name <- "sbc_logistic_config.yml"
} else{
  config_file_name <- args[1]
}
config <- yaml::read_yaml(file.path(config_dir, config_file_name), eval.expr=TRUE)

n_sims_sbc_s <- config$n_sims_sbc_s # Number of different s-steps for SBC
n_sims_sbc_i <- config$n_sims_sbc_i # Number of i-steps per s-step
keep_fits <- config$keep_fits

surrogate_model <- config$surrogate_model
seed <- config$seed
fct <- logistic_fct
n_sims <- config$n_sims
sigma_sims <- config$sigma_sim
n_exp <- config$n_exp
lower_w_exp <- config$lower_w_exp
upper_w_exp <- config$upper_w_exp
use_clustering <- config$use_clustering
number_clusters <- config$number_clusters
adapt_delta <- config$adapt_delta
refresh_cmdstan <- config$refresh_cmdstan
chains <- config$chains
iter_sampling <- config$iter_sampling
iter_warmup <- config$iter_warmup
propagate_sigma_t <- config$propagate_sigma_a

set.seed(seed)
## PCE##
poly_degree <- config$poly_degree
M <- config$M

## Files ##
stan_code_dir <- file.path(here(),"src/stan_code")
stan_models_logistic <- get_stan_models_logistic(stan_code_dir, surrogate_model,
                                                 propagate_sigma_t)
imodel_mp <- cmdstan_model(stan_models_logistic$imodel_mp)
imodel_ml <- cmdstan_model(stan_models_logistic$imodel_ml)
imodel_mlp <- cmdstan_model(stan_models_logistic$imodel_mlp)
smodel_file <- stan_models_logistic$smodel

## Data ##
sdatasets <- list()
for (s_id in c(1:n_sims_sbc_s)){
  sdata <- get_1d_sdata(129, poly_degree, M, fct, sigma_sims[1], prior_only = 0,
                        fixed_w_sim = FALSE)
  sdatasets_id <- list()
  for (n_sim in n_sims){
    indices <- get_training_data_indices(129, n_sim)
    sdata_tmp <- sdata
    sdata_tmp$N_sim <- n_sim
    sdata_tmp$y_sim <- sdata$y_sim[indices]
    sdata_tmp$w_sim <- matrix(sdata$w_sim[indices], ncol=M)
    sdatasets_id[[length(sdatasets_id)+1]] <- sdata_tmp
  }
  sdatasets[[length(sdatasets)+1]] <- sdatasets_id
}

plot_dirs <- list()
for (n_sim_id in seq_along(n_sims)){
  n_sim <- n_sims[n_sim_id]
  experiment_details <- paste0("_nsim_", n_sim, "_sigmasim_", sigma_sims[1],
                               "_ncluster_", number_clusters,
                               "_useclustering_",use_clustering,
                               "_polydegree_", poly_degree,
                               "sbc_s", n_sims_sbc_s,
                               "sbc_i", n_sims_sbc_i,
                               "prop_sigt", propagate_sigma_t)
  plot_dir <- get_experiment_dir_name("plots", "logistic_case_sbc",
                                      surrogate_model, experiment_details)
  dir.create(plot_dir, showWarnings = FALSE)
  file.copy(file.path(config_dir, config_file_name), plot_dir)
  plot_dirs[[length(plot_dirs)+1]] <- plot_dir
}
smodel_fits <- list()
for (s_id in c(1:n_sims_sbc_s)){
  smodel_fits_id <- list()
  for (n_sim_id in seq_along(n_sims)){
    n_sim <- n_sims[n_sim_id]
    sdata <- sdatasets[[s_id]][[n_sim_id]]
    smodel_fit_list <- get_smodel_fit(sdata, smodel_file, fct, plot_pp=TRUE,
                                      model_title="True Model: \nlogistic (4params)",
                                      iter_sampling=iter_sampling, chains = chains,
                                      adapt_delta=0.999)
    smodel_fit <- smodel_fit_list[[1]]
    (p_pp_surr <- smodel_fit_list[[2]])
    smodel_fits_id[[length(smodel_fits_id)+1]] <- smodel_fit
    ggsave(file.path(plot_dirs[[n_sim_id]], paste0("smodel", s_id, "_", n_sim, ".png")), width=10, height=4)
  }
  smodel_fits[[length(smodel_fits)+1]] <- smodel_fits_id
}

idatasets <- list()
for (idx in c(1:(n_sims_sbc_s*n_sims_sbc_i))){
  w_exp_gt <- get_w_exp(1, mean=0, sd=0.5, lower=-1, upper=1)$w_exp_1
  sigma_exp <- runif(1, min=0.005, max=0.05)
  idata <- get_1d_idata(n_exp, w_exp_gt, logistic_fct, sigma_exp)
  idata$sigma_exp <- sigma_exp
  idatasets[[length(idatasets)+1]] <- idata
}

library(future)
plan(multisession, gc=TRUE, workers=10)
for (j in seq_along(sigma_sims)){
  sigma_sim <- sigma_sims[[j]]
  for (i in seq_along(n_sims)){
    plot_dir <- plot_dirs[[i]]
    n_sim <- n_sims[[i]]
    i_step_datasets <- list()
    for (s_step in c(1:n_sims_sbc_s)){
      # S-step
      smodel_fit <- smodel_fits[[s_step]][[i]]
      for (i_step in c(1:n_sims_sbc_i)){
        # I-step dataset
        w_exp_gt <- idatasets[[length(i_step_datasets) + 1]]$w_exp_1[1]
        y_exp <- idatasets[[length(i_step_datasets) + 1]]$y_exp
        sigma_exp <- idatasets[[length(i_step_datasets) + 1]]$sigma_exp
        if (use_clustering){
          idata_method <- "weighted_cluster_draws"
        } else{
          idata_method <- "thin"
        }
        i_step_generator <- SBC_generator_function(i_step_1d_generator_single, n_exp, poly_degree, M,
                                                   smodel_fit, fct,
                                                   method=idata_method,
                                                   number_draws=number_clusters,
                                                   upper=1, lower=-1,
                                                   fixed_w_exp=w_exp_gt,
                                                   fixed_y_exp=y_exp,
                                                   sigma_exp=sigma_exp,
                                                   prior_only=0,
                                                   kmeans.nstart=25, kmeans.iter.max=100)
        i_step_dataset <- generate_datasets(
          i_step_generator,
          1)
        i_step_datasets[[length(i_step_datasets) + 1]] <- i_step_dataset
      }
    }
    i_step_dataset <- do.call(bind_datasets, i_step_datasets)
    
    # i-step c mean
    i_step_backend <- SBC_backend_mean(
      imodel_mp, iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains, adapt_delta = adapt_delta)
    results_mean <- compute_SBC(i_step_dataset, i_step_backend,
                                globals = c("SBC_fit.SBC_backend_mean", "combine_args"),
                                cache_mode = "results", cache_location = file.path(plot_dir, paste0("true_logistic_nsim_", n_sim, "_sigmasim_", sigma_sim, "_cmean")),
                                keep_fits = keep_fits)
    
    # i-step c median
    i_step_backend <- SBC_backend_median(
      imodel_mp, iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains, adapt_delta = adapt_delta)
    results_median <- compute_SBC(i_step_dataset, i_step_backend,
                                globals = c("SBC_fit.SBC_backend_median", "combine_args"),
                                cache_mode = "results", cache_location = file.path(plot_dir, paste0("true_logistic_nsim_", n_sim, "_sigmasim_", sigma_sim, "_cmedian")),
                                keep_fits = keep_fits)

    # i-step c single trial
    i_step_backend <- SBC_backend_singletrial(
      imodel_ml, iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains, adapt_delta = adapt_delta)
    
    results_st <- compute_SBC(i_step_dataset, i_step_backend,
                              globals = c("SBC_fit.SBC_backend_singletrial", "combine_args", "get_cluster_draws"),
                              cache_mode = "results", cache_location = file.path(plot_dir, paste0("true_logistic_nsim_", n_sim, "_sigmasim_", sigma_sim, "_cml")),
                              keep_fits = keep_fits)
    
    # i-step c mix-log-post
    i_step_backend <- SBC_backend_singletrial(
      imodel_mlp, iter_warmup = iter_warmup, iter_sampling = iter_sampling, chains = chains, adapt_delta = adapt_delta)
    results_mlp <- compute_SBC(i_step_dataset, i_step_backend,
                               globals = c("SBC_fit.SBC_backend_singletrial", "combine_args", "get_cluster_draws"),
                               cache_mode = "results", cache_location = file.path(plot_dir, paste0("true_logistic_nsim_", n_sim, "_sigmasim_", sigma_sim, "_cmlp")),
                               keep_fits = keep_fits)
    
    # SBC multi trial
    if (use_clustering){
      iter_sampling_mp <- iter_sampling
      chains_mp <- chains
    } else{
      iter_sampling_mp <- chains*iter_sampling/length(i_step_dataset$generated[[1]]$c_0)
      chains_mp <- 1
    }
    i_step_backend <- SBC_backend_multitrial(
      imodel_mp, iter_warmup = iter_warmup, iter_sampling = iter_sampling_mp, chains = chains_mp, adapt_delta = adapt_delta)
    results_mt <- compute_SBC(i_step_dataset, i_step_backend,
                              globals = c("SBC_fit.SBC_backend_multitrial", "combine_args",
                                          "get_cluster_draws", "here",
                                          "SBC_fit_to_draws_matrix.list", "read_cmdstan_csv",
                                          "subset_draws", "resample_draws", "as_draws_matrix"),
                              cache_mode = "results", cache_location = file.path(plot_dir, paste0("true_logistic_nsim_", n_sim, "_sigmasim_", sigma_sim, "_cmp")),
                              keep_fits = keep_fits, thin_ranks = 1)
  }
}
