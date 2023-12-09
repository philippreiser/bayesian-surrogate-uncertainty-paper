library(here)
library(tidyverse)
library(ggridges)
library(bayesplot)
source(file.path(here(), "src/models/fit_2-step.R"))
source(file.path(here(), "src/utils/utils.R"))
source(file.path(here(), "src/eval/run_eval_two_step_linear_helpers.R"))
source(file.path(here(), "src/models/analytic_model.R"))
source(file.path(here(), "src/plotting/plotting_analytic_model.R"))
args = commandArgs(trailingOnly=TRUE)
config_dir <- file.path(here(), "config/analytic_case_study")
if (length(args)==0){
  config_file_name <- "intercept_slope_config.yml"
} else{
  config_file_name <- args[1]
}
config <- yaml::read_yaml(file.path(config_dir, config_file_name), eval.expr=TRUE)

seed <- config$seed
set.seed(seed)
plot_dir <- get_experiment_dir_name("plots", "analytic_case", config$model)
if (!dir.exists(plot_dir)){
  dir.create(plot_dir, recursive=TRUE)
}
file.copy(file.path(config_dir, config_file_name), plot_dir)

dirs <- get_dirs()
stan_code_dir <- dirs[[1]]
model <- config$model
eval_analytic <- config$eval_analytic
eval_mcmc <- config$eval_mcmc
eval_e_post <- TRUE
eval_e_lik <- TRUE

stan_files <- get_stan_files(stan_code_dir, model)

omega_t <- config$omega_t
c_gt <- config$c_gt
y_t <- get_y_model(omega_t, c_gt, model)
mu_t0 <- config$mu_t0
sigma_t0 <- config$sigma_t0
mu_i0 <- config$mu_i0
sigma_i0 <- config$sigma_i0
sigma_i <- config$sigma_i

omega_i_limits <- c(config$omega_i_range[1], config$omega_i_range[2])
omega_i_range <- seq(from=config$omega_i_range[1], to=config$omega_i_range[2], config$omega_i_range[3])

plots_list <- list()
eval_i_step_df_list <- list()
for (sigma_a in config$sigma_a_values){
  for (omega_i_gt in config$omega_i_gt_values){
    
    #### T-Step ####
    if (eval_analytic == TRUE || eval_mcmc == TRUE){
      t_posterior_analytic_parameters <- get_t_posterior_analytic_parameters(omega_t, y_t, 
                                                                 mu_t0, sigma_t0,
                                                                 sigma_a, model)
      sigma_t1 <- t_posterior_analytic_parameters$sigma_t1
      mu_t1 <- t_posterior_analytic_parameters$mu_t1
      if (model == "slope_cubic" || model ==  "intercept_slope"){
        # Set integration bounds for c (check what happens if sigma_t1 has large
        # off-diagonal elements)
        config$numint_lower_c_e_post <- c(mu_t1[1] - 3*sqrt(sigma_t1[1,1]),
                                    mu_t1[2] - 3*sqrt(sigma_t1[2,2]))
        config$numint_upper_c_e_post <- c(mu_t1[1] + 3*sqrt(sigma_t1[1,1]),
                                    mu_t1[2] + 3*sqrt(sigma_t1[2,2]))
      }
      plot_analytic_t_posterior(model, config$numint_lower_c_e_post,
                                config$numint_upper_c_e_post, t_posterior_densities)
    }
    
    if (eval_mcmc == TRUE){
      t_model_fit <- get_t_posterior_mcmc(omega_t, y_t, mu_t0, sigma_t0,
                                         sigma_a, stan_files$t_model, seed,
                                         config$adapt_delta_t_step)
    }
    
    #### I-step (point) ####
    method = "point"
    y_i <- get_y_model(omega_i_gt, c_gt, model)
    
    i_prior_densities <- dnorm(omega_i_range, mean=mu_i0, sd=sigma_i0)
    if (eval_analytic == TRUE){
      print(plot_analytic_posterior_predictive(omega_i_range, mu_t1, sigma_a,
                                               sigma_t1, sigma_i, model))
      if (model == "slope" || model == "intercept" || model == "intercept_slope"){
        point_i_posterior_analytic_parameters <- get_point_i_posterior_analytic(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, model)
        mu_i1_point <- point_i_posterior_analytic_parameters$mu_i1_point
        sigma_i1_point <- point_i_posterior_analytic_parameters$sigma_i1_point
        posterior_densities_point <- dnorm(omega_i_range, mean=mu_i1_point, sd=sigma_i1_point)
      } else if (model == "cubic" || model == "slope_cubic"){
        posterior_densities_point <- get_point_i_posterior_numeric(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, model, omega_i_range)
      }
      (p_i_point_analytic <- plot_analytic_posterior(posterior_densities_point, omega_i_range, "point"))
    }
    if (eval_mcmc == TRUE){
      imodel_fit <- get_point_i_posterior_mcmc(y_i, mu_i0, sigma_i0, sigma_i, mu_t1,
                                         stan_files$i_model_point, seed,
                                         config$adapt_delta_point)
      
      p_point <- plot_mcmc_analytic_posterior(
        rename_variables(imodel_fit$draws("omega_i"), "point\n(mcmc)"="omega_i"),
        posterior_densities_point, omega_i_range, omega_i_limits, omega_i_gt, method)
      print(p_point)
      print(mu_i1_point- mean(imodel_fit$draws("omega_i")))
      print(mu_i1_point-sd(imodel_fit$draws("omega_i")))
    }
    
    #### I-step (E-Lik) ####
    if (eval_e_lik == TRUE){
      method <- "E-Lik\n(analytic)"
      
      if (eval_analytic == TRUE){
        source(file.path(here(), "src/models/analytic_model.R"))
        if (model == "slope" || model == "intercept_slope"){
          posterior_densities_e_lik <- get_e_lik_i_posterior_numeric_2(y_i, mu_i0, sigma_i0, sigma_i,
                                                     mu_t1, sigma_t1, model,
                                                     i_prior_densities, omega_i_range,
                                                     config$numint_lower_w_e_post, config$numint_upper_w_e_post,
                                                     config$r_tol)
        } else if (model=="quadratic" || model == "cubic" || model == "slope_cubic" || model=="intercept"){
          if (model=="quadratic" || model == "cubic" || model=="intercept"){
            config$numint_lower_c_e_post <- mu_t1-(3*sigma_t1)
            config$numint_upper_c_e_post <- mu_t1+(3*sigma_t1)
          }
          posterior_densities_e_lik <- get_e_lik_i_posterior_numeric(y_i, mu_i0, sigma_i0, sigma_i,
                                                     mu_t1, sigma_t1, model,
                                                     i_prior_densities, omega_i_range,
                                                     config$numint_lower_c_e_post, config$numint_upper_c_e_post,
                                                     config$numint_lower_w_e_post, config$numint_upper_w_e_post,
                                                     config$r_tol)
        }
        print(paste0("Integral over E_Lik: ", sum(posterior_densities_e_lik*(omega_i_range[2]-omega_i_range[1]))))
        plot_analytic_posterior(posterior_densities_e_lik, omega_i_range, method="E-Lik\n(analytic)")+
          geom_vline(xintercept = omega_i_gt)
      }
      
      if (eval_mcmc == TRUE){
        imodel_e_lik_fit <- get_e_lik_i_posterior_mcmc(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, sigma_t1,
                                            stan_files$i_model_e_lik, seed, model,
                                            config$adapt_delta_e_lik)
        
        p_e_lik <- plot_mcmc_analytic_posterior(
          rename_variables(imodel_e_lik_fit$draws("omega_i"), "E-Lik\n(mcmc)"="omega_i"),
          posterior_densities_e_lik, omega_i_range, omega_i_limits, omega_i_gt, method)
        print(p_e_lik)
        samples_analytic_posterior <- sample(omega_i_range, prob=posterior_densities_e_lik, size=1e4, replace=TRUE)
        print(mean(samples_analytic_posterior)- mean(imodel_e_lik_fit$draws("omega_i")))
        print(sd(samples_analytic_posterior)-sd(imodel_e_lik_fit$draws("omega_i")))
      }
    } else{
      posterior_densities_e_lik <- rep(0, length(omega_i_range))
    }
    
    #### I-step (E-Post) ####
    if (eval_e_post == TRUE){
      method <- "E-Post\n(analytic)"
      if (eval_analytic == TRUE){
        source(file.path(here(), "src/models/analytic_model.R"))
        posterior_densities_e_post <- get_e_post_i_posterior_numeric(y_i, mu_i0, sigma_i0, sigma_i,
                                                   mu_t1, sigma_t1, model,
                                                   i_prior_densities, omega_i_range,
                                                   config$numint_lower_c_e_post, config$numint_upper_c_e_post,
                                                   config$numint_lower_w_e_post, config$numint_upper_w_e_post,
                                                   config$r_tol)
        print(paste0("Integral over E_Post: ", sum(posterior_densities_e_post*(omega_i_range[2]-omega_i_range[1]))))
        plot_analytic_posterior(posterior_densities_e_post, omega_i_range, method="E-Post\n(analytic)")+
          geom_vline(xintercept = omega_i_gt)
      }
      
      if (eval_mcmc == TRUE){
        ### -------- MCMC ---------- ###
        imodel_e_post_fit <- get_e_post_i_posterior_mcmc(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, sigma_t1,
                                            stan_files$i_model_point, seed, model,
                                            n_c_skip=config$n_c_skip,
                                            config$adapt_delta_e_post)
        p_posterior_mix <- plot_mcmc_analytic_posterior(
          rename_variables(subset_draws(imodel_e_post_fit, variable="omega_i[1,1]"), "E-Post\n(mcmc)"="omega_i"),
          posterior_densities_e_post, omega_i_range, omega_i_limits, omega_i_gt, method)
        print(p_posterior_mix)
        ### ------------------------ ###
      }
    } else{
      posterior_densities_e_post <- rep(0, length(omega_i_range))
    }
    
    #### i-step (E-Log-Lik) ####
    method <- "E-Log-Lik\n(analytic)"
    if (eval_analytic == TRUE){
      source(file.path(here(), "src/models/analytic_model.R"))
      posterior_densities_e_log_lik <- get_e_log_lik_i_posterior_analytic(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, sigma_t1, model, i_prior_densities, omega_i_range)
      plot_analytic_posterior(posterior_densities_e_log_lik, omega_i_range, method="E-Log-Lik\n(analytic)")+
        geom_vline(xintercept = omega_i_gt)    }
    if (eval_mcmc == TRUE){
      ### -------- MCMC ---------- ###
      imodel_e_log_lik_fit <- get_e_log_lik_i_posterior_mcmc(y_i, mu_i0, sigma_i0, sigma_i, 
                                            mu_t1, sigma_t1, stan_files$i_model_e_log_lik, 
                                            seed, model, config$adapt_delta_e_log_lik)
      p_e_log_lik <- plot_mcmc_analytic_posterior(
        rename_variables(imodel_e_log_lik_fit$draws("omega_i"), "E-Log-Lik\n(mcmc)"="omega_i"),
        posterior_densities_e_log_lik, omega_i_range, omega_i_limits, omega_i_gt, method)
      print(p_e_log_lik)
      samples_analytic_posterior <- sample(omega_i_range, prob=posterior_densities_e_log_lik, size=1e4, replace=TRUE)
      print(mean(samples_analytic_posterior)- mean(imodel_e_log_lik_fit$draws("omega_i")))
      print(sd(samples_analytic_posterior)-sd(imodel_e_log_lik_fit$draws("omega_i")))
      ### ------------------------ ###
    }
    
    if (eval_analytic == TRUE){
      eval_i_step_df_non_tidy <- get_i_eval_df_non_tidy(omega_i_range, posterior_densities_point, posterior_densities_e_lik, posterior_densities_e_post, posterior_densities_e_log_lik, sigma_a, omega_i_gt)
      eval_i_step_df_list[[length(eval_i_step_df_list)+1]] <- eval_i_step_df_non_tidy
      print(paste0("Area under curve, point: ", sum(posterior_densities_point), " m-l: ", sum(posterior_densities_e_lik), " m-p: ", sum(posterior_densities_e_post)))
    }
    if (eval_mcmc == TRUE){
      p_all <- grid.arrange(p_point, p_e_lik, p_posterior_mix, p_e_log_lik, nrow=4)
      plots_list[[length(plots_list)+1]] <- p_all
    }
  }
  if (eval_mcmc == TRUE){
    p_all_posteriors_all_omega_is <- do.call(grid.arrange, c(plots_list, ncol=length(plots_list)))
    ggsave(file.path(plot_dir, paste0("mcmc_analytic_posteriors_sigma_a", sigma_a, ".pdf")), p_all_posteriors_all_omega_is, width=15, height=8)
    ggsave(file.path(plot_dir, paste0("mcmc_analytic_posteriors_sigma_a", sigma_a, ".png")), p_all_posteriors_all_omega_is, width=15, height=8)
    
  }
}

eval_i_step_df_non_tidy <- do.call(rbind, eval_i_step_df_list)
saveRDS(eval_i_step_df_non_tidy, file.path(plot_dir, "analytic_posteriors.rds"))
colors <- c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')

appender <- function(value){
  return(TeX(paste("$sigma_A = $", value)))
}
eval_i_step_df_non_tidy %>%
  select(omega_i, point, e_log_lik, e_lik, e_post, sigma_a, omega_i_gt)%>%
  rename(`Point`=point, `E-Log-Lik`=e_log_lik, `E-Lik`=e_lik,
         `E-Post`=e_post)%>%
  pivot_longer(c(`Point`, `E-Log-Lik`, `E-Lik`, `E-Post`), names_to = "method", values_to = "posterior")%>%
  ggplot()+
  geom_line(aes(x=omega_i, y=posterior, colour=method), size=1.5)+
  scale_color_manual(values = colors)+
  facet_wrap(vars(sigma_a), scales = "fixed", labeller = as_labeller(appender, 
                         default = label_parsed))+
  theme_bw()+
  xlim(-1, 1)+
  theme(legend.position = "bottom")+
  labs(x= TeX(r'($\omega_{I}$)'), y="I-posterior density")+
  theme(text = element_text(size = 41))
  
ggsave(file.path(plot_dir, paste0("analytic_posteriors.pdf")), width=30, height=10)
ggsave(file.path(plot_dir, paste0("analytic_posteriors.png")), width=30, height=10)
