library(here)
library(patchwork)
library(xtable)
source(file.path(here(), "src/utils/utils.R"))
source(file.path(here(), "src/data_generation/get_data.R"))
source(file.path(here(), "src/utils/true_models_helpers.R"))
source(file.path(here(), "src/eval/run_eval_two_step_logistic_helpers.R"))
source(file.path(here(), "src/models/fit_2-step.R"))

library(ggh4x)
theme_set(theme_bw())

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  config_file_name <- "logistic_config.yml"
} else{
  config_file_name <- args[1]
}
config_dir <- file.path(here(), "config/logistic_case_study")
config <- yaml::read_yaml(file.path(config_dir, config_file_name), eval.expr=TRUE)

### Experiment parameters
surrogate_model <- config$surrogate_model
seed <- config$seed
fct <- logistic_fct
n_sims <- config$n_sims
n_sim_all <- n_sims[length(n_sims)]
sigma_sim <- config$sigma_sim
n_exp <- config$n_exp
w_exp_gts <- config$w_exp_gts
sigma_exp <- config$sigma_exp
lower_w_exp <- config$lower_w_exp
upper_w_exp <- config$upper_w_exp
poly_degree <- config$poly_degree
M <- config$M
use_clustering <- config$use_clustering
number_clusters <- config$number_clusters
adapt_delta <- config$adapt_delta
refresh_cmdstan <- config$refresh_cmdstan
chains <- config$chains
iter_sampling <- config$iter_sampling
iter_warmup <- config$iter_warmup
propagate_sigma_a <- config$propagate_sigma_a
experiment_details <- paste0("_nsim_", n_sim_all, "_nexp_", n_exp, "_sigmaexp_",
                             sigma_exp, "_sigmasim_", sigma_sim,"_ncluster_",
                             number_clusters, "_useclustering_",use_clustering,
                             "_polydegree_", poly_degree, "_propagatesigmaa_",
                             propagate_sigma_a)

### Create results directory and store config
plot_dir <- get_experiment_dir_name("plots", "logistic_case_densities",
                                    surrogate_model, experiment_details)
if (!dir.exists(plot_dir)){
  dir.create(plot_dir, recursive = TRUE)
}
file.copy(file.path(config_dir, config_file_name), plot_dir)

### Load Stan files
stan_code_dir <- file.path(here(),"src/stan_code")
stan_models_logistic <- get_stan_models_logistic(stan_code_dir, surrogate_model,
                                                 config$propagate_sigma_a)
file_imodel_mp <- stan_models_logistic$imodel_mp
file_imodel_ml <- stan_models_logistic$imodel_ml
file_imodel_mlp <- stan_models_logistic$imodel_mlp
smodel_file <- stan_models_logistic$smodel

### Generate Datasets
set.seed(seed)
sdatasets <- get_1d_sdatasets(n_sims, poly_degree, M, fct, sigma_sim)
idatasets <- lapply(w_exp_gts, function(w_exp_gt) get_1d_idata(n_exp, w_exp_gt, fct, sigma_exp))

### Perform two-step procedure with different I-Steps
posterior_df <- data.frame()
sigma_exp_df <- data.frame()
p_pp_surrs <- list()
p_n_sim <- list()
for (idx_sim in seq_along(n_sims)){
  # columns of fig 1 (vary n_sim)
  n_sim <- n_sims[[idx_sim]]
  sdata <- sdatasets[[idx_sim]]
  smodel_fit_list <- get_smodel_fit(sdata, smodel_file, fct, plot_pp=TRUE,
                                    model_title="True Model: \nlogistic (4params)",
                                    hide_legend=FALSE, title=TeX(sprintf(r'($N_{T} = %d$)', n_sim)),
                                    adapt_delta=adapt_delta,
                                    chains=4, iter_sampling=250,
                                    linewidth=1.5,
                                    propagate_sigma_a=propagate_sigma_a)
  smodel_fit <- smodel_fit_list[[1]]
  (p_pp_surr <- smodel_fit_list[[2]])
  p_pp_surrs[[length(p_pp_surrs) + 1]] <- p_pp_surr
  for (idx_exp in seq_along(w_exp_gts)){
    # rows of fig 1 (vary w_exp_gt)
    w_exp_gt <- w_exp_gts[[idx_exp]]
    idata <- idatasets[[idx_exp]]
    posterior_list <- run_mean_ml_mlp_mp(idata, smodel_fit, number_clusters, n_exp, w_exp_gt,
                                           poly_degree, M, fct, adapt_delta,
                                           file_imodel_ml, file_imodel_mp, file_imodel_mlp,
                                         iter_sampling, iter_warmup, chains,
                                         use_clustering=use_clustering)
    posterior_df_new <- posterior_list$w_posterior
    sigma_exp_df_new <- posterior_list$sigma_posterior
    posterior_df <- rbind(posterior_df, posterior_df_new)
    sigma_exp_df <- rbind(sigma_exp_df, sigma_exp_df_new)
  }
}

### Save results
if (!file.exists(file.path(plot_dir, paste0("posterior_df", seed, ".rds")))){
  saveRDS(posterior_df, file.path(plot_dir, paste0("posterior_df", seed, ".rds")))
  saveRDS(sigma_exp_df, file.path(plot_dir, paste0("sigma_exp_df", seed, ".rds")))
  saveRDS(p_pp_surrs, file.path(plot_dir, paste0("p_pp_surrs", seed, ".rds")))
}

### Load results
# uncomment to load pre-computed posterior draws
# posterior_df <- readRDS(file.path(plot_dir, paste0("posterior_df", seed, ".rds")))
# p_pp_surrs <- readRDS(file.path(plot_dir, paste0("p_pp_surrs", seed, ".rds")))

### Plot and save plots
p_iposterior_densities <- plot_iposterior_densities(posterior_df, p_pp_surrs, w_exp_gts, n_sims,
                                                    linewidth = 1.5)
p_iposterior_densities
ggsave(file.path(plot_dir, paste0("pp_prop_sigma_t_", surrogate_model, "_seed", seed, "_axis.pdf")), width=30, height=16)
