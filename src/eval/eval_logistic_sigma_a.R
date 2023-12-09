library(here)
library(patchwork)
library(xtable)
source(file.path(here(), "src/utils/utils.R"))
source(file.path(here(), "src/data_generation/get_data.R"))
source(file.path(here(), "src/eval/run_eval_two_step_logistic_helpers.R"))
library(ggh4x)
theme_set(theme_bw())

surrogate_model <- "pce"
date_prop_true <- "2023-09-18_13-49-56"
date_prop_false <- "2023-09-18_11-15-45"
experiment_details <- "_nsim_20_nexp_5_sigmaexp_0.01_sigmasim_0.01_ncluster_1000_useclustering_FALSE_polydegree_5"
plot_dir_prop_true <- get_experiment_dir_name("plots", 
                                              "logistic_case_densities",
                                              surrogate_model,
                                              paste0(experiment_details,
                                                     "_propagatesigmaa_TRUE"),
                                              date_prop_true)
plot_dir_prop_false <- get_experiment_dir_name("plots", 
                                              "logistic_case_densities",
                                              surrogate_model,
                                              paste0(experiment_details,
                                                     "_propagatesigmaa_FALSE"),
                                              date_prop_false)
config_file_name <- "pce_config.yml"
config_prop_true <- yaml::read_yaml(file.path(plot_dir_prop_true, 
                                              "pce_propsigt_TRUE_config.yml"),
                                    eval.expr=TRUE)
config_prop_false <- yaml::read_yaml(file.path(plot_dir_prop_false, 
                                               "pce_propsigt_FALSE_config.yml"),
                                    eval.expr=TRUE)
posterior_df_prop_true <- readRDS(file.path(plot_dir_prop_true, paste0("posterior_df", config_prop_true$seed, ".rds")))
posterior_df_prop_true$propagate_sigma_t <- TRUE
posterior_df_prop_false <- readRDS(file.path(plot_dir_prop_false, paste0("posterior_df", config_prop_false$seed, ".rds")))
p_pp_surrs <- readRDS(file.path(plot_dir_prop_true, paste0("p_pp_surrs",config_prop_true$seed, ".rds")))
p_pp_surrs_post <- post_process_t_step_pp_plot(p_pp_surrs)
posterior_df <- rbind(posterior_df_prop_true, posterior_df_prop_false)
linewidth <- 1.5
appender_w_exp_gt <- function(value){
  return(TeX(paste("$\\omega_I^* = $", value)))
}
prop_sigma_a_labels <- list(
  'FALSE'=TeX('$\\theta=c$', output="character"),
  'TRUE'=TeX('$\\theta=(c, \\sigma_A)$', output="character")
)
prop_sigma_a_labeller <- function(variable,value){
  return(prop_sigma_a_labels[value])
}

appender_n_sim <- function(value){
  return(TeX(paste("$\\N_T = $", value)))
}
to_string_n_sim <- as_labeller(appender_n_sim, default = label_parsed)
to_string_prop_sigma_a <- as_labeller(prop_sigma_a_labeller, default = label_parsed)
#levels(posterior_df$propagate_sigma_t) <- c("FALSE" = TeX('$\\theta$'),
#                                            "TRUE" = TeX('$\\theta$'))
p_i <- posterior_df %>%
  filter(w_exp_gt == -0.05)%>%
  filter(n_sim %in% c(10, 20))%>%
  rename("Point" = "median", "E-Post" = "E-post", 
         "E-Log-Lik" = "E-log-lik", "E-Lik" = "E-lik")%>%
  pivot_longer(c("Point", "E-Post", "E-Log-Lik", "E-Lik"), names_to = "I-Step")%>%
  ggplot()+
  scale_color_viridis_d()+
  geom_density(aes(value, color=`I-Step`), linewidth=linewidth)+
  #theme(axis.text.y=element_blank(),
  #      axis.ticks.y=element_blank())+
  facet_grid2(vars(n_sim), vars(propagate_sigma_t), 
              labeller = labeller(n_sim = to_string_n_sim,
                                  propagate_sigma_t = to_string_prop_sigma_a),
              scales="free", independent="y")+
  #facet_wrap(vars(n_sim, w_exp_gt), nrow=length(n_sims), scales="free", labeller = "label_both")+
  xlab(TeX(r'($\omega_{I}$)'))+
  ylab("I-posterior distribution")+
  scale_color_manual(values = colors)+
  theme(legend.position = "bottom")+
  theme(text = element_text(size = 41),
        axis.text = element_text(size = 20))+
  geom_vline(aes(xintercept = w_exp_gt), linewidth=linewidth) +
  xlim(c(-0.2, 0.05)) # w_exp_gt == -0.05
  # xlim(c(0.0, 0.2)) # w_exp_gt == 0.1
  # xlim(c(0.2, 0.4)) # w_exp_gt == 0.3
  # xlim(c(0.3, 0.5)) # w_exp_gt == 0.4
p_i
((p_pp_surrs_post[[2]]/p_pp_surrs_post[[3]]) | p_i)+
  plot_layout(guides = 'collect', widths = c(1, 3)) &
  theme(legend.position = "bottom")
ggsave(file.path(plot_dir_prop_true, paste0("i_posterior_densities_", surrogate_model, "_compare_sigma_a", config_prop_true$seed, "_axis.pdf")), width=30, height=11)
