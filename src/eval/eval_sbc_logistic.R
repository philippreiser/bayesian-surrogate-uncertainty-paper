library(SBC)
library(here)
library(xtable)
source(file.path(here(),"src/models/fit_2-step.R"))
source(file.path(here(), "src/utils/utils.R"))
source(file.path(here(), "src/eval/eval_sbc_logistic_helpers.R"))
source(file.path(here(), "src/eval/sbc_logistic_backend.R"))
theme_set(theme_bw())

# figure: sbc gamma
surrogate_model <- "true_model"
if (surrogate_model == "true_model"){
  n_sims <- c(5, 6, 7, 8, 9, 10)
} else{
  n_sims <- c(10, 15, 20, 30, 40)
}
sigma_sims <- c(0.01)
n_sims_sbc_s <- 10 # Number of different s-steps for SBC
n_sims_sbc_i <- 20 # Number of i-steps per s-step
propagate_sigma_t <- FALSE
use_clustering <- TRUE
date <- format(Sys.time(), "%Y-%m-%d")

base_path <- file.path(here(), "plots/logistic_case_sbc", surrogate_model)
head_path <- file.path(base_path, paste0(date, "wsimfix_", surrogate_model, "_"))

i_methods <- c("Point", "E-Lik", "E-Log-Lik", "E-Post")
metrics <- get_sbc_results(base_path, n_sims, sigma_sims, date,
                           surrogate_model, i_methods, n_sims_sbc_s,
                           n_sims_sbc_i, use_clustering)
prop_sigma_a_labels <- list(
  'FALSE'=TeX('$\\theta=c$', output="character"),
  'TRUE'=TeX('$\\theta=(c, \\sigma_A)$', output="character")
)
prop_sigma_a_labeller <- function(variable,value){
  return(prop_sigma_a_labels[value])
}
to_string_prop_sigma_a <- as_labeller(prop_sigma_a_labeller, default = label_parsed)
p_log_gamma <- metrics %>%
  ggplot(aes(x=n_sim, y=gamma_sbc, color=`method`))+
  geom_point(show.legend = FALSE)+
  geom_line(show.legend = FALSE)+
  scale_color_manual(values=colors)+
  geom_hline(yintercept = 0)+
  ylab(TeX(r'($\log(\gamma)$)'))+
  xlab(TeX(r'($N_T$)'))+
  theme(legend.position = "bottom")
if (surrogate_model == "pce"){
  p_log_gamma <- p_log_gamma+
    facet_wrap(~propagate_sigma_t,  scales = "free_y", 
               labeller = labeller(propagate_sigma_t = to_string_prop_sigma_a),
               nrow=1)
}
p_log_gamma

# end
plot_dir <- get_plot_dir(base_path, n_sims[1], propagate_sigma_t, date,
                         n_sims_sbc_s, n_sims_sbc_i, use_clustering)

# figure: sbc ecdf plots
get_results_df <- function(results, method, n_sim, propagate_sigma_t){
  stats_sigma_exp <- results$result$stats%>%
    filter(variable=="sigma_exp")
  stats_sigma_exp_sim_value <- stats_sigma_exp$simulated_value
  stats_w_exp <- results$result$stats%>%
    filter(variable=="w_exp[1,1]")%>%
    mutate(variable="w",
           n_sim=n_sim,
           method=method,
           propagate_sigma_t=propagate_sigma_t)
  stats_w_exp$simulated_value_sigma_exp <- stats_sigma_exp_sim_value
  return(stats_w_exp)
}
if (surrogate_model == "true_model"){
  n_sims <- c(5, 10)
} else{
  n_sims <- c(10)
}
if (surrogate_model == "pce"){
  propagate_sigma_t_vec = c(TRUE, FALSE)
} else{
  propagate_sigma_t_vec = c(FALSE)
}
results_stats_list <- list()
for (j in seq_along(n_sims)){
  n_sim <- n_sims[[j]]
  for (k in seq_along(propagate_sigma_t_vec)){
    propagate_sigma_t <- propagate_sigma_t_vec[[k]]
    plot_dir <- get_plot_dir(base_path, n_sim, propagate_sigma_t, date,
                             n_sims_sbc_s, n_sims_sbc_i, use_clustering)
    results_file_name_start <- paste0(plot_dir, "/true_logistic_nsim_"
                                      , n_sim, "_sigmasim_", sigma_sims[[1]])
    results_median <- readRDS(paste0(results_file_name_start, "_cmedian.rds"))
    results_median_stats <- get_results_df(results_median, "Point", n_sim,
                                           propagate_sigma_t)
    results_ml <- readRDS(paste0(results_file_name_start, "_cml.rds"))
    results_ml_stats <- get_results_df(results_ml, "E-Lik", n_sim,
                                       propagate_sigma_t)
    results_mlp <- readRDS(paste0(results_file_name_start, "_cmlp.rds"))
    results_mlp_stats <- get_results_df(results_mlp, "E-Log-Lik", n_sim,
                                        propagate_sigma_t)
    results_mp <- readRDS(paste0(results_file_name_start, "_cmp.rds"))
    results_mp_stats <- get_results_df(results_mp, "E-Post", n_sim,
                                       propagate_sigma_t)
    results <- list(results_median_stats, results_ml_stats, results_mlp_stats,
                    results_mp_stats)
    results_stats_list[[length(results_stats_list)+1]] <- results
  }
}
all_sbc_stats <- bind_rows(results_stats_list)
results_mp_stats$simulated_value[[1]]
plot_i_posteriors <- function(index){
  w_exp_gt <- results_mp$stats$simulated_value[[index]]
  ((mcmc_dens(results_mp$result$fits[[index]], "w_exp[1,1]")+ggtitle("E-Post")+geom_vline(xintercept = w_exp_gt))+
      (mcmc_dens(results_ml$result$fits[[index]]$draws(), "w_exp[1,1]")+ggtitle("E-Lik"))+geom_vline(xintercept = w_exp_gt))/
    ((mcmc_dens(results_median$result$fits[[index]]$draws(), "w_exp[1,1]")+ggtitle("Median")+geom_vline(xintercept = w_exp_gt))+
       (mcmc_dens(results_mlp$result$fits[[index]]$draws(), "w_exp[1,1]")+ggtitle("E-Log-Lik"))+geom_vline(xintercept = w_exp_gt))
}

p_all_ecdf_plot <- plot_ecdf_diff_facett_method_nsim(all_sbc_stats, i_methods, n_sims, surrogate_model)
p_all_ecdf_plot
ggsave(file.path(plot_dir, paste0("true_logistic_nsim_ecdf_", surrogate_model, ".png")), p_all_ecdf_plot, width = 9, height = 5)
ggsave(file.path(plot_dir, paste0("true_logistic_nsim_ecdf_", surrogate_model, ".pdf")), p_all_ecdf_plot, width = 9, height = 5)


appender <- function(value){
  print(value)
  if ((substring(value[1], 1, 1)=="E")||(substring(value[1], 1, 1)=="m")){
    return(value)
  } else if ((substring(value[1], 1, 1)=="T")||(substring(value[1], 1, 1)=="F")){
    return(prop_sigma_a_labels[value])
  } else{
    return(TeX(paste("$\\N_T = $", value)))
  }
}
p_sharpness <- all_sbc_stats %>%
  mutate(width_ci_90 = q95-q5)%>%
  group_by(method, n_sim, propagate_sigma_t)%>%
  summarise(mean_width_ci_90 = mean(width_ci_90),
            mean_width_ci_90_ci = 1.96 * sd(width_ci_90) / sqrt(n_sims_sbc_s*n_sims_sbc_i),
            q_width_ci_25 = quantile(width_ci_90, 0.025),
            q_width_ci_975 = quantile(width_ci_90, 0.975))%>%
  ggplot(aes(x=n_sim, y=mean_width_ci_90, color=method))+
  scale_color_manual(values=colors)+
  geom_point()+
  geom_line()+
  ylab("Sharpness (90% CI)")+
  xlab(TeX(r'($N_T$)'))+
  theme(legend.position = "bottom")
if (surrogate_model == "pce"){
  p_sharpness <- p_sharpness+
    facet_wrap(~propagate_sigma_t, nrow = 1, scales = "free_y",
               labeller = as_labeller(appender, 
                                      default = label_parsed))
}
p_sharpness

p_all_ecdf_gamma <- (p_all_ecdf_plot / p_log_gamma / p_sharpness)+
  plot_layout(heights = c(2, 1, 1)) &
  theme(legend.position = "bottom")
p_all_ecdf_gamma
ggsave(file.path(plot_dir, paste0("ecdf_gamma_", surrogate_model, ".png")), p_all_ecdf_gamma, width = 9, height = 9)
ggsave(file.path(plot_dir, paste0("ecdf_gamma_", surrogate_model, ".pdf")), p_all_ecdf_gamma, width = 9, height = 9)


