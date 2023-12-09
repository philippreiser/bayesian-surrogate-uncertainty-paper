library(ggplot2)
library(patchwork)
library(latex2exp)

colors <- c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')
plot_pp_1d <- function(smodel, fit_smodel, fct, l_poly_coeffs_mat, comb,
                       n_test=100, xmax=1, xmin=-1,
                       model_title="PCE", n_sim=0,
                       prior_only=0, hide_legend=FALSE,
                       title=NULL, sdata=NULL,
                       linewidth=1, propagate_sigma_a=FALSE
                       ){
  # Posterior predictions for [-1, 1]
  var_dim <- "w_sim"
  p <- sdata$d
  M <- 1
  x1 <- seq(from=-1, to=1,length.out=n_test)
  tmp_1 <- data.frame(x1 = x1) %>%
    mutate(
      N_sim = n_test,
      y_sim = fct(x1),
      w_sim = x1
    ) 
  sdata_1 <- list(
    N_sim = n_test,
    sigma_sim_lower = 0.0,
    d = p,
    M = M,
    w_sim = matrix(c(tmp_1$w_sim), ncol=M),
    y_sim = tmp_1$y_sim,
    l_poly_coeffs = t(l_poly_coeffs_mat),
    N_comb = nrow(comb),
    comb = comb,
    prior_only = prior_only
  )
  fct_str <- paste0(as.character(substitute(fct)))[1]
  fit_gq <- smodel$generate_quantities(fit_smodel, data = sdata_1, seed = 100, parallel_chains=2)
  mu <- as_draws_df(fit_gq$draws(c("mu_pred_sim"))) %>%
    as_draws_df() %>%
    as_tibble() %>%
    select(starts_with("mu")) %>%
    apply(2, quantile, c(0.05, 0.5, 0.95)) %>%
    t() %>% 
    data.frame(w_sim = tmp_1[[var_dim]], .)  %>% 
    gather(pct, y_sim, -w_sim)
  y_rep <- as_draws_df(fit_gq$draws(c("y_rep_sim"))) %>%
    as_draws_df() %>%
    as_tibble() %>%
    select(starts_with("y_rep")) %>%
    apply(2, quantile, c(0.05, 0.95)) %>%
    t() %>% 
    data.frame(w_sim = tmp_1[[var_dim]], .)  %>% 
    gather(pct_rep, y_sim, -w_sim)
  if (model_title == "PCE"){
    model_title = paste0("PCE d = ", p)
  } else{
    model_title = "T-Step"
  }
  pfit_1 <- ggplot() +
    geom_point(aes(x=w_sim, y=y_sim, color="Training Data"),
               data = data.frame(w_sim=sdata$w_sim, y_sim=sdata$y_sim), size=linewidth*2)+
    geom_line(aes_string(x=var_dim, y="y_sim"), color="gray", linetype = "dashed", data = tmp_1) +
    geom_line(aes(x=w_sim, y_sim, linetype = pct, color="Mean of Posterior Predictive"), data = mu,
              linewidth=linewidth)+
    #geom_ribbon(aes(x=w_sim, y_))
    scale_linetype_manual(values = c(2,1,2)) +
    labs(y = TeX(r'($y_{T}$)'), x= TeX(r'($\omega_{T}$)')) +
    ylim(-1.6, 1.6) +
    xlim(xmin, xmax)+
    guides(linetype = "none")
    #scale_color_manual(name=model_title,values=c(TeX(r'($D_T$)')='black', 'Mean of Posterior Predictive'='red',
    #                                             'Posterior Predictive'='#8b0000'))+
    
  if (propagate_sigma_a){
    pfit_1 <- pfit_1+
      geom_line(aes(x=w_sim, y_sim, linetype = pct_rep, 
                    color="Posterior Predictive"), 
                data = y_rep, linewidth=linewidth)+
      scale_color_manual(name=model_title,
                         values=c('Training Data'='black',
                                  'Mean of Posterior Predictive'='red',
                                  'Posterior Predictive'='#8b0000'))
  } else {
    pfit_1 <- pfit_1+
      scale_color_manual(name=model_title,
                         values = c('Mean of Posterior Predictive'='red',
                                    'Training Data'='black'))
  }
  if (hide_legend == TRUE){
    pfit_1 <- pfit_1 + theme(legend.position = "none")
  }
  if (!is.null(title)){
    pfit_1 <- pfit_1 + ggtitle(title)
  }
  pfit_1
}

post_process_t_step_pp_plot <- function(p_pp_surrs){
  for (idx in seq_along(p_pp_surrs)){
    p_pp_surrs[[idx]] <- p_pp_surrs[[idx]]+
      theme(text = element_text(size = 41),
            axis.text = element_text(size = 20),
            plot.title = element_blank())
  }
  return(p_pp_surrs)
}

plot_imodel_posterior <- function(imodel_fit, w_exp_gt, xmin=-1, xmax=1, title="I-Step", color="blue"){
  color_scheme_set(scheme = color)
  if ('post_warmup_draws' %in% names(imodel_fit)){
    imodel_draws <- imodel_fit$post_warmup_draws
  }
  else if ('draws' %in% names(imodel_fit)){
    imodel_draws <- imodel_fit$draws()
  }
  else {
    imodel_draws <- imodel_fit
  }
  plot <- mcmc_dens(imodel_draws, regex_pars="w_exp")+
    vline_at(w_exp_gt)+plot_bg(fill = "white")+xlim(xmin, xmax)+
    labs(title = title)
  plot
}

plot_imodel_metrics <- function(plot_dir, eval_df){
  colors <- c("yellow", "green", "blue", "cyan")
  ggplot(data=eval_df)+
    geom_point(mapping = aes(x = w_exp_gt, y = rmse, color = inf_method))+
    geom_smooth(mapping = aes(x = w_exp_gt, y = rmse, color = inf_method), method = "gam", se = TRUE)+
    scale_color_manual(values=colors)
  ggsave(file.path(plot_dir, paste0("rmse.png")))
  ggplot(data=eval_df)+
    geom_point(mapping = aes(x = w_exp_gt, y = mae, color = inf_method))+
    geom_smooth(mapping = aes(x = w_exp_gt, y = mae, color = inf_method), method = "gam", se = TRUE)+
    scale_color_manual(values=colors)
  ggsave(file.path(plot_dir, paste0("mae.png")))
  ggplot(data=eval_df)+
    geom_point(mapping = aes(x = w_exp_gt, y = bias, color = inf_method))+
    geom_smooth(mapping = aes(x = w_exp_gt, y = bias, color = inf_method), method = "gam", se = TRUE)+
    scale_color_manual(values=colors)
  ggsave(file.path(plot_dir, paste0("bias.png")))
  ggplot(data=eval_df)+
    geom_point(mapping = aes(x = w_exp_gt, y = s66, color = inf_method))+
    geom_smooth(mapping = aes(x = w_exp_gt, y = s66, color = inf_method), method = "gam", se = TRUE)+
    scale_color_manual(values=colors)
  ggsave(file.path(plot_dir, paste0("sharpness66.png")))
  ggplot(data=eval_df)+
    geom_point(mapping = aes(x = w_exp_gt, y = w_gt_in_ci66, color = inf_method))+
    geom_smooth(mapping = aes(x = w_exp_gt, y = w_gt_in_ci66, color = inf_method), method = "glm", method.args=list(family="binomial"), se = TRUE)+
    scale_color_manual(values=colors)
  ggsave(file.path(plot_dir, paste0("w_gt_in_ci66.png")))
}
plot_imodel_posterior_single <- function(posterior_df, w_exp_gt_query, n_sim_query){
  if (w_exp_gt_query == 0.05){
    x_lim <- c(0.04, 0.1)
  } else if (w_exp_gt_query == 0.1){
    x_lim <- c(0.07, 0.12)
  } else if (w_exp_gt_query == 0.2){
    x_lim <- c(0.15, 0.25)
  } else if (w_exp_gt_query == 0.3){
    x_lim <- c(0.2, 0.4)
  } else if (w_exp_gt_query == 0.4){
    x_lim <- c(0.3, 0.41)
  } else if (w_exp_gt_query == 0.5){
    x_lim <- c(0.3, 0.6)
  } else if (w_exp_gt_query == -0.05){
    x_lim <- c(-0.2, 0.05) # without sigma_sim: c(-0.15, 0)
  } else{
    x_lim <- NULL
  }
  plot_posterior_1 <- posterior_df %>%
    filter(w_exp_gt==w_exp_gt_query, n_sim==n_sim_query)%>%
    pivot_longer(c("mean", "median", "E-post", "E-log-lik", "E-lik"), names_to = "i-posterior")%>%
    
    ggplot()+
    scale_color_viridis_d()+
    geom_vline(xintercept = w_exp_gt_query)+
    geom_density(aes(value, color=`i-posterior`))+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    xlab(TeX(r'($\omega_{I}$)'))+
    ylab("posterior")+
    scale_color_manual(values = colors)
  if (!is.null(x_lim)){
    plot_posterior_1 <- plot_posterior_1+
      xlim(x_lim[1], x_lim[2])
  }
  plot_posterior_1
}
plot_sigma_exp_single <- function(posterior_df, w_exp_gt_query, n_sim_query, sigma_exp_gt){
  plot_posterior_1 <- posterior_df %>%
    filter(w_exp_gt==w_exp_gt_query, n_sim==n_sim_query)%>%
    rename("Point" = "median", "E-Post" = "E-post", 
           "E-Log-Lik" = "E-log-lik", "E-Lik" = "E-lik")%>%
    pivot_longer(c("Point", "E-Post", "E-Log-Lik", "E-Lik"), names_to = "I-Step")%>%
    ggplot()+
    scale_color_viridis_d()+
    geom_vline(xintercept = sigma_exp_gt)+
    geom_density(aes(value, color=`I-Step`))+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    xlab(TeX(r'($\sigma_{I}$)'))+
    ylab("posterior")
  plot_posterior_1
}
plot_imodel_posterios_ggplot_fig1 <- function(posterior_df, p_pp_surrs,
                                              w_exp_gts, n_sims){
  plot_posteriors <- list()
  for (n_sim in n_sims){
    for (w_exp_gt in w_exp_gts){
      plot_posterior <- plot_imodel_posterior_single(posterior_df, w_exp_gt,
                                                     n_sim)
      plot_posteriors[[length(plot_posteriors) + 1]] <- plot_posterior
    }
  }
  if (length(n_sims)==3){
    p_all_posteriors <- (p_pp_surrs[[1]] | plot_posteriors[[1]] | plot_posteriors[[2]] | plot_posteriors[[3]]) /
      (p_pp_surrs[[2]] | plot_posteriors[[4]] | plot_posteriors[[5]] | plot_posteriors[[6]]) /
      (p_pp_surrs[[3]] | plot_posteriors[[7]] | plot_posteriors[[8]] | plot_posteriors[[9]]) /
      plot_layout(guides = 'collect') &
      theme(legend.position = "bottom")
  } else if (length(n_sims)==2){
    p_all_posteriors <- (p_pp_surrs[[1]] | plot_posteriors[[1]] | plot_posteriors[[2]]) /
      (p_pp_surrs[[2]] | plot_posteriors[[3]] | plot_posteriors[[4]]) /
      plot_layout(guides = 'collect') &
      theme(legend.position = "bottom")
  }
  p_all_posteriors
}

plot_iposterior_densities <- function(posterior_df, p_pp_surrs, w_exp_gts,
                                      n_sims, linewidth=1.5){
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
    pivot_longer(c("Point", "E-Post", "E-Log-Lik", "E-Lik"), names_to = "I-Step")%>%
    ggplot()+
    scale_color_viridis_d()+
    geom_density(aes(value, color=`I-Step`), linewidth=linewidth)+
    #theme(axis.text.y=element_blank(),
    #      axis.ticks.y=element_blank())+
    facet_grid2(vars(n_sim), vars(w_exp_gt), 
                labeller = labeller(w_exp_gt = to_string_w_exp_gt,
                                    n_sim = to_string_n_sim),
                scales="free", independent="y")+
    geom_vline(aes(xintercept = w_exp_gt), linewidth=linewidth) +
    xlab(TeX(r'($\omega_{I}$)'))+
    ylab("I-posterior distribution")+
    scale_color_manual(values = colors)+
    theme(legend.position = "bottom")+
    theme(text = element_text(size = 41),
          axis.text = element_text(size = 20))
  if (surrogate_model == "true_model"){
    p_i <- p_i + 
      scale_x_facet(w_exp_gt == -0.05, limits = c(-0.1, 0.0))+
      scale_x_facet(w_exp_gt == 0.1, limits = c(0.04, 0.14))+
      scale_x_facet(w_exp_gt == 0.3, limits = c(0.25, 0.35))+
      scale_x_facet(w_exp_gt == 0.4, limits = c(0.1, 0.8))
  } else if (surrogate_model == "pce"){
    if (propagate_sigma_a)
    p_i <- p_i + 
      scale_x_facet(w_exp_gt == -0.05, limits = c(-0.2, 0.05))+ # 
      scale_x_facet(w_exp_gt == 0.1, limits = c(0.0, 0.3))+ # 
      scale_x_facet(w_exp_gt == 0.3, limits = c(0.2, 0.4))+
      scale_x_facet(w_exp_gt == 0.4, limits = c(0.1, 1)) #
  }
  p_pp_surrs_post <- post_process_t_step_pp_plot(p_pp_surrs)
  ((p_pp_surrs_post[[1]]/p_pp_surrs_post[[2]]/p_pp_surrs_post[[3]]) | p_i)+
    plot_layout(guides = 'collect', widths = c(1, 3)) &
    theme(legend.position = "bottom")
}

plot_sigma_exp_posteriors <- function(posterior_df, p_pp_surrs,
                                              w_exp_gts, n_sims, sigma_exp_gt){
  plot_posteriors <- list()
  for (n_sim in n_sims){
    for (w_exp_gt in w_exp_gts){
      plot_posterior <- plot_sigma_exp_single(posterior_df, w_exp_gt,
                                                     n_sim, sigma_exp_gt)
      plot_posteriors[[length(plot_posteriors) + 1]] <- plot_posterior
    }
  }
  if (length(n_sims)==3){
    p_all_posteriors <- (p_pp_surrs[[1]] | plot_posteriors[[1]] | plot_posteriors[[2]] | plot_posteriors[[3]]) /
      (p_pp_surrs[[2]] | plot_posteriors[[4]] | plot_posteriors[[5]] | plot_posteriors[[6]]) /
      (p_pp_surrs[[3]] | plot_posteriors[[7]] | plot_posteriors[[8]] | plot_posteriors[[9]]) /
      plot_layout(guides = 'collect') &
      theme(legend.position = "bottom")
  } else if (length(n_sims)==2){
    p_all_posteriors <- (p_pp_surrs[[1]] | plot_posteriors[[1]] | plot_posteriors[[2]]) /
      (p_pp_surrs[[2]] | plot_posteriors[[3]] | plot_posteriors[[4]]) /
      plot_layout(guides = 'collect') &
      theme(legend.position = "bottom")
  }
  p_all_posteriors
}

plot_samplingtime_over_nclusters_poly_degree_facett <- function(
    fit_diagnostics_df, chains, iter_warmup, iter_sampling, plot_dir=NULL,
    log_scale=FALSE, width=16, height=6, scale_sampling_time_e_post=FALSE
    ){
  if (scale_sampling_time_e_post){
    fit_diagnostics_df <- fit_diagnostics_df %>%
      mutate(
        sampling_and_warmup_time=if_else(method == "mix-post",
              warmup_time + sampling_time/n_clusters,sampling_and_warmup_time))
  }
  plot_time <- fit_diagnostics_df %>%
    mutate(method = replace(method, method=="mean", "Point"),
           method = replace(method, method=="mix-lik", "E-Lik"),
           method = replace(method, method=="mix-log-post", "E-Log-Lik"),
           method = replace(method, method=="mix-post", "E-Post"))%>%
    mutate(`I-posterior`=method, d=poly_degree)%>%
    group_by(`I-posterior`, n_clusters, d) %>%
    summarise(
      sampling_and_warmup_time_mean=mean(sampling_and_warmup_time)/60,
      sampling_and_warmup_time_std=sd(sampling_and_warmup_time)/60)%>%
    ggplot(aes(x=n_clusters, y=sampling_and_warmup_time_mean, color=`I-posterior`))+
    scale_color_manual(values=colors)+
    geom_line()+
    geom_point()+
    geom_linerange(aes(ymin=sampling_and_warmup_time_mean-sampling_and_warmup_time_std,
                       ymax=sampling_and_warmup_time_mean+sampling_and_warmup_time_std))+
    xlab(TeX(r'($K$)'))+
    ylab(paste0('Sampling time(', chains, ' chains; ', iter_warmup, '\nwarmup ', iter_sampling+iter_warmup,' iterations) [min]'))+
    facet_wrap(~d, labeller = labeller(.rows=label_both), nrow=1)+
    {if (log_scale) scale_y_continuous(trans="log10")}+
    theme(legend.position = "bottom")+
    theme(text = element_text(size = 20))
  if (!is.null(plot_dir)){
    if (log_scale){
      ggsave(file.path(plot_dir, "sampling_log_times.pdf"), plot_time, width=width, height=height)
    } else{
      ggsave(file.path(plot_dir, "sampling_times.pdf"), plot_time, width=width, height=height)
    }
  }
  return(plot_time)
}

plot_samplingtimeratio_over_nclusters_poly_degree_facett <- function(
    fit_diagnostics_df, chains, iter_warmup, iter_sampling, plot_dir=NULL,
    log_scale=FALSE, width=8, height=4){
  fit_diagnostics_ratio_df <- fit_diagnostics_df %>%
    mutate(`i-posterior`=method)%>%
    group_by(poly_degree)%>%
    mutate(base_time_mean = mean(sampling_and_warmup_time[`i-posterior`=="mean"]), base_time_std = sd(sampling_and_warmup_time[`i-posterior`=="mean"]))%>%
    mutate(time_ratio = sampling_and_warmup_time/base_time_mean)%>%
    ungroup()%>%
    group_by(`i-posterior`, n_clusters, poly_degree) %>%
    summarise(
      time_ratio_mean=mean(time_ratio),
      time_ratio_std=sd(time_ratio))
  diagnostics_df_simple <- fit_diagnostics_ratio_df %>%
    filter(`i-posterior`=="mix-post", poly_degree==2)
  diagnostics_df_simple.lm <- lm(formula = time_ratio_mean ~ n_clusters,
                                    data = diagnostics_df_simple)
  plot_time_ratio <- fit_diagnostics_ratio_df %>%
    ggplot(aes(x=n_clusters, y=time_ratio_mean, color=`i-posterior`))+
    scale_colour_viridis_d()+
    geom_point()+
    geom_linerange(aes(ymin=time_ratio_mean-time_ratio_std,
                       ymax=time_ratio_mean+time_ratio_std))+
    xlab(TeX(r'($K$)'))+
    ylab(paste0('Time(K)/Time(1)'))+
    stat_smooth(method=lm, formula = y ~ x)+
    facet_wrap(~poly_degree, labeller = labeller(.rows=label_both), nrow=1)+
    {if (log_scale) scale_y_continuous(trans="log10")}+
    theme(legend.position = "bottom")
  plot_time_slope <- fit_diagnostics_ratio_df %>%
    group_by(poly_degree, `i-posterior`) %>%
    do({
      mod = lm(time_ratio_mean ~ n_clusters, data = .)
      data.frame(Intercept = coef(mod)[1],
                 Slope = coef(mod)[2])
    })%>%
    mutate(Slope = ifelse(`i-posterior` == "mix-post", 1, Slope))%>%# normalized to slope of mix-post
    ggplot(aes(x=poly_degree, y=Slope, color=`i-posterior`))+
    {if (log_scale) scale_y_continuous(trans="log10")}+
    geom_point()+
    geom_line()+
    ylab("Relative slope of sampling\ntimes normalized to Mix-Post")+
    scale_colour_viridis_d()
    # theme(text = element_text(size = 20))

  if (!is.null(plot_dir)){
    if (log_scale){
      ggsave(file.path(plot_dir, "sampling_times_ratio_log_scale.pdf"), plot_time_ratio, width=width, height=height)
      ggsave(file.path(plot_dir, "sampling_times_slope_log_scale.pdf"), plot_time_slope, width=width, height=height)
    } else{
      ggsave(file.path(plot_dir, "sampling_times_ratio_scale.pdf"), plot_time_ratio, width=width, height=height)
      ggsave(file.path(plot_dir, "sampling_times_slope_scale.pdf"), plot_time_slope, width=width, height=height)
    }
  }
  return(plot_time_slope)
}
#' Creates an ecdf_diff-plot using facetting on method and nsim 
#' using as basis code from
#' https://github.com/hyunjimoon/SBC/blob/dbbed5dedd3ff0befe25db9bb27f2bf4d0c488eb/R/plot.R#L224
#'
#' @param x ecdf results tibble.
#' @param variables list of strings of the variables of interest contained in x.
#' @returns A ggplot object.
#' @examples
plot_ecdf_diff_facett_method_nsim <- function(x, methods, n_sims,
                                              surrogate_model,
                                              variables = NULL,
                                              K = NULL,
                                              gamma = NULL,
                                              prob = 0.95,
                                              size = 1,
                                              alpha = 0.33,
                                              ...,
                                              parameters = NULL) {
  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }
  ecdf_df_list <- list()
  limits_df_trans_list <- list()
  for (i in seq_along(n_sims)){
    for(j in seq_along(methods)){
      for (k in seq_along(propagate_sigma_t_vec)){
        x_idx <- x %>%
          filter(n_sim==n_sims[[i]], method==methods[[j]], 
                 propagate_sigma_t==propagate_sigma_t_vec[[k]])%>%
          select(-c(n_sim, method, propagate_sigma_t))
        ecdf_data <-
          data_for_ecdf_plots(x_idx, variables = variables,
                              prob = prob, K = K, gamma = gamma, ...)
        
        N <- ecdf_data$N
        K <- ecdf_data$K
        z <- ecdf_data$z
        
        ecdf_df <- dplyr::mutate(ecdf_data$ecdf_df, z_diff = ecdf - z, type = "sample ECDF")
        limits_df <- ecdf_data$limits_df
        limits_df_trans <- data.frame(
          x = c(0, rep(z[2:(K + 1)], each = 2)),
          ymax = limits_df$upper / N - c(rep(z[1:K], each = 2), 1),
          ymin = limits_df$lower / N - c(rep(z[1:K], each = 2), 1),
          type = "theoretical CDF"
        )
        ecdf_df$n_sim <- n_sims[[i]]
        ecdf_df$method <- methods[[j]]
        ecdf_df$propagate_sigma_t <- propagate_sigma_t_vec[[k]]
        limits_df_trans$n_sim <- n_sims[[i]]
        limits_df_trans$method <- methods[[j]]
        limits_df_trans$propagate_sigma_t <- propagate_sigma_t_vec[[k]]
        ecdf_df_list[[length(ecdf_df_list)+1]] <- ecdf_df
        limits_df_trans_list[[length(limits_df_trans_list)+1]] <- limits_df_trans
      }
    }
  }
  ecdf_dfs <- bind_rows(ecdf_df_list)
  limits_df_transs <- bind_rows(limits_df_trans_list)
  
  if (surrogate_model=="true_model"){
    row_variable <- "n_sim"
  } else{
    row_variable <- "propagate_sigma_t"
  }
  prop_sigma_a_labels <- list(
    'FALSE'=TeX('$\\theta=c$', output="character"),
    'TRUE'=TeX('$\\theta=(c, \\sigma_A)$', output="character")
  )
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
  ggplot(ecdf_dfs, aes(color = type, fill = type)) +
    geom_ribbon(
      data = limits_df_transs,
      aes(x = x, ymax = ymax, ymin = ymin),
      alpha = alpha,
      size = size) +
    geom_step(
      aes(x = z, y = z_diff)
    ) +
    scale_color_manual(
      name = "",
      values = rlang::set_names(
        c("skyblue1", "black"),
        c("theoretical CDF", "sample ECDF")),
      labels = c(
        "theoretical CDF" = expression(italic("theoretical CDF")),
        "sample ECDF" = expression(italic("sample ECDF"))
      )
    ) +
    scale_fill_manual(
      name = "",
      values = c("theoretical CDF" = "skyblue",
                 "sample ECDF" = "transparent"),
      labels = c(
        "theoretical CDF" = expression(italic("theoretical CDF")),
        "sample ECDF" = expression(italic("sample ECDF"))
      )
    ) +
    xlab(NULL) +
    ylab(NULL) +
    # TODO: case for true_model! (c("n_sim", "method"))
    facet_grid(c(row_variable, "method"), scales = "free_y",
               #labeller = labeller(.rows = label_both, .cols = label_value),
               labeller = as_labeller(appender, 
                                      default = label_parsed))+
    theme(legend.position = "bottom")
}

data_for_ecdf_plots <- function(x, ...,
                                prob = 0.95,
                                gamma = NULL,
                                K = NULL
) {
  UseMethod("data_for_ecdf_plots")
}


data_for_ecdf_plots.SBC_results <- function(x, variables = NULL,
                                            prob = 0.95,
                                            gamma = NULL,
                                            K = NULL,
                                            parameters = NULL) {
  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }
  
  data_for_ecdf_plots(x$stats, variables = variables, prob = prob, gamma = gamma, K = K)
}


data_for_ecdf_plots.data.frame <- function(x, variables = NULL,
                                           prob = 0.95,
                                           gamma = NULL,
                                           K = NULL,
                                           max_rank = x$max_rank,
                                           parameters = NULL) {
  
  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }
  
  if("parameter" %in% names(x)) {
    if(!("variable" %in% names(x))) {
      warning("The x parameter contains a `parameter` column, which is deprecated, use `variable` instead.")
      x$variable <- x$parameter
    }
  }
  
  if("dataset_id" %in% names(x)) {
    if(!("sim_id" %in% names(x))) {
      warning("The x parameter contains a `dataset_id` column, which is deprecated, use `sim_id` instead.")
      x$sim_id <- x$dataset_id
    }
  }
  
  
  if(!all(c("variable", "rank", "sim_id") %in% names(x))) {
    stop(SBC_error("SBC_invalid_argument_error",
                   "The stats data.frame needs a 'variable', 'rank' and 'sim_id' columns"))
  }
  
  stats <- x
  if(!is.null(variables)) {
    stats <- dplyr::filter(stats, variable %in% variables)
  }
  
  if(is.null(max_rank)) {
    stop("max_rank either has to be supplied explicitly or be a column in the data")
  }
  max_rank <- unique(max_rank)
  if(length(max_rank) > 1) {
    stop("Differing max_rank across variables is not supported yet.")
  }
  
  summ <- dplyr::summarise(dplyr::group_by(stats, variable), count = dplyr::n(), .groups = "drop")
  if(length(unique(summ$count)) > 1) {
    stop("Not all variables have the same number of simulations.")
  }
  
  rank <- dplyr::select(stats, sim_id, variable, rank)
  rank_matrix <- tidyr::pivot_wider(rank, names_from = "variable",
                                    values_from = "rank")
  rank_matrix <- as.matrix(dplyr::select(rank_matrix, -sim_id))
  
  
  data_for_ecdf_plots(rank_matrix, max_rank = max_rank, prob = prob,
                      gamma = gamma, K = K)
}

data_for_ecdf_plots.matrix <- function(x,
                                       max_rank,
                                       variables = NULL,
                                       prob = 0.95,
                                       gamma = NULL,
                                       K = NULL,
                                       size = 1,
                                       alpha = 0.33,
                                       parameters = NULL) {
  
  if(!is.null(parameters)) {
    warning("The `parameters` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- parameters
    }
  }
  
  ranks_matrix <- x
  if(any(!is.finite(ranks_matrix))) {
    stop("Ranks may only contain finite values")
  }
  
  if(!is.null(variables)) {
    ranks_matrix <- ranks_matrix[, variables]
  }
  
  pit <- ranks_to_empirical_pit(ranks_matrix, max_rank)
  N <- nrow(pit)
  if (is.null(K)) {
    K <- min(max_rank + 1, N)
  }
  if (is.null(gamma)) {
    gamma <- SBC:::adjust_gamma(
      N = N,
      L = 1,
      K = K,
      conf_level = prob
    )
  }
  limits_df <- as.data.frame(SBC:::ecdf_intervals(
    N = N,
    L = 1,
    K = K,
    gamma = gamma))
  z <- seq(0,1, length.out = K + 1)
  
  ecdf_vals <- apply(pit, 2, function(col) ecdf(col)(z))
  
  ecdf_df <- as.data.frame(ecdf_vals)
  ecdf_df$..z <- z
  ecdf_df <- tidyr::pivot_longer(ecdf_df, -..z, names_to = "variable", values_to = "ecdf")
  ecdf_df <- dplyr::rename(ecdf_df, z = ..z)
  
  structure(list(limits_df = limits_df, ecdf_df = ecdf_df, K = K, N = N, z = z),
            class = "SBC_ecdf_data")
}
ranks_to_empirical_pit <- function(ranks, n_posterior_samples){
  (1 + ranks) / (1 + n_posterior_samples)
}
