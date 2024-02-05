library(here)
library(patchwork)
source(file.path(here(), "src/utils/utils.R"))
source(file.path(here(), "src/data_generation/get_data.R"))
source(file.path(here(), "src/utils/true_models_helpers.R"))
source(file.path(here(), "src/eval/run_eval_two_step_logistic_helpers.R"))
source(file.path(here(), "src/models/fit_2-step.R"))

# library(ggh4x)
theme_set(theme_bw())

### Load config
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
  config_file_name <- "logistic_example_config.yml"
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


### plot T-/I-data
data.frame(
  w_sim = c(sdatasets[[1]]$w_sim),
  y_sim = c(sdatasets[[1]]$y_sim),
  w_exp = c(idatasets[[1]]$w_exp_1),
  y_exp = c(idatasets[[1]]$y_exp)
) %>%
  ggplot()+
  geom_point(aes(x = w_sim, y = y_sim))+
  geom_point(aes(x = w_exp, y = w_exp), color = "blue")

### T-Step
idx_sim <- 1
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
# Plot posterior predictive
(p_pp_surr <- smodel_fit_list[[2]])

### I-Step
idx_exp <- 1
w_exp_gt <- w_exp_gts[[idx_exp]]
idata <- idatasets[[idx_exp]]
posterior_list <- run_mean_ml_mlp_mplp_mp(idata, smodel_fit, number_clusters, n_exp, w_exp_gt,
                                     poly_degree, M, fct, adapt_delta,
                                     file_imodel_ml, file_imodel_mp, file_imodel_mlp,
                                     iter_sampling, iter_warmup, chains,
                                     use_clustering=TRUE)
posterior_df <- posterior_list$w_posterior
appender_w_exp_gt <- function(value){
  return(TeX(paste("$\\omega_I^* = $", value)))
}
appender_n_sim <- function(value){
  return(TeX(paste("$\\N_T = $", value)))
}
to_string_w_exp_gt <- as_labeller(appender_w_exp_gt, default = label_parsed)
to_string_n_sim <- as_labeller(appender_n_sim, default = label_parsed)
p_i <- posterior_df %>%
  rename("Point" = "median", "E-Post" = "E-post", 
         "E-Log-Lik" = "E-log-lik", "E-Lik" = "E-lik")%>%
  pivot_longer(c("Point", "E-Post"), names_to = "I-Step")%>%
  ggplot()+
  geom_density(aes(value, color=`I-Step`))+
  facet_grid2(vars(n_sim), vars(w_exp_gt), 
              labeller = labeller(w_exp_gt = to_string_w_exp_gt,
                                  n_sim = to_string_n_sim),
              scales="free", independent="y")+
  geom_vline(aes(xintercept = w_exp_gt)) +
  xlab(TeX(r'($\omega_{I}$)'))+
  ylab("I-posterior distribution")+
  scale_color_manual(values = colors)+
  theme(legend.position = "bottom")+
  xlim(-0.13, 0.0)

p_i
