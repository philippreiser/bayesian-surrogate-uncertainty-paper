plot_mcmc_analytic_posterior <- function(model_fit_draws, post_plot_analytic, 
                             x_plot, x_lim_plot, x_gt, method){
  mcmc_areas_ridges(model_fit_draws)+
    xlim(if (is.null(x_lim_plot)) c(NA, NA) else x_lim_plot)+
    geom_vline(xintercept = x_gt)+
    geom_ridgeline(data = data.frame(
      "omega_i" = x_plot,
      "posterior"=post_plot_analytic,
      "method"=method
    ),
    aes(x=omega_i, y=rep(1, length(posterior)),height=1.8*posterior),
    fill="transparent", linetype=2
    )
}

plot_analytic_posterior <- function(posterior_plot, omega_i_plot, method){
  ggplot(data = data.frame(
    "omega_i" = omega_i_plot,
    "posterior"=posterior_plot,
    "method"=method
  ))+
    geom_ridgeline(aes(x=omega_i,y=rep(0, length(omega_i)), height=posterior))+
    ylab(method)+
    theme(axis.title.y = element_text(angle = 0, hjust=1, vjust=0))
}

plot_analytic_posterior_predictive <- function(omega_i_plot, mu_s1, sigma_s,
                                               sigma_s1, sigma_i, model){
  get_pp_mean <- function(omega_i){
    if (model == "slope"){
      return(omega_i*mu_s1)
    } else if (model == "intercept_slope"){
      W_i <- get_design_matrix(omega_i, model)
      return(W_i%*%mu_s1)
    }
  }
  get_pp_025 <- function(omega_i){
    if (model == "slope"){
      pp_025 <- omega_i*mu_s1 - 2*sqrt(omega_i_plot^2*sigma_s1^2 + sigma_s^2)
    } else if (model == "intercept_slope"){
      W_i <- get_design_matrix(omega_i, model)
      pp_025 <- W_i%*%mu_s1 - 2*sqrt(rowSums(W_i%*%sigma_s1*W_i) + sigma_s^2)
    }
    return(pp_025)
  }
  get_pp_975 <- function(omega_i){
    if (model == "slope"){
      pp_975 <- omega_i*mu_s1 + 2*sqrt(omega_i_plot^2*sigma_s1^2 + sigma_s^2)
    } else if (model == "intercept_slope"){
      W_i <- get_design_matrix(omega_i, model)
      pp_975 <- W_i%*%mu_s1 + 2*sqrt(rowSums(W_i%*%sigma_s1*W_i) + sigma_s^2)
    }
    return(pp_975)
  }
  posterior_predictive = data.frame(
    "q500" = get_pp_mean(omega_i_plot),
    "q025" = get_pp_025(omega_i_plot),
    "q975" = get_pp_975(omega_i_plot),
    "w"    = omega_i_plot)
  posterior_predictive%>%
    ggplot(aes(x=w))+
    geom_line(aes(y=q500))+
    geom_line(aes(y = q025), lty = 2) +
    geom_line(aes(y = q975), lty = 2) +
    theme_classic()+
    ylab("y")
}

plot_analytic_t_posterior <- function(model, lower_bounds, upper_bounds, t_posterior_densities){
  if (model == "slope" || model == "cubic" || model == "intercept"){
    c_range <- seq(from=-10, to=6, 0.05)
  }
  if (model == "slope_cubic" || model ==  "intercept_slope"){
    c_range <- expand.grid(c_1 = seq(lower_bounds[1],
                                     upper_bounds[1],
                                     length.out=100),
                           c_2 = seq(lower_bounds[2],
                                     upper_bounds[2],
                                     length.out=100))
  } 
  if (model == "slope" || model == "cubic" || model=="intercept"){
    t_posterior_densities <- dnorm(c_range, mean=mu_t1, sd=sigma_t1)
  } else if (model == "slope_cubic" || model ==  "intercept_slope"){
    t_posterior_densities <- dmvnorm(c_range, mu_t1, sigma_t1)
    data.frame(c_1=c_range$c_1, c_2=c_range$c_2, t_posterior=t_posterior_densities)%>%
      ggplot()+
      geom_contour_filled(aes(c_1, c_2, z=t_posterior))
  }
}
