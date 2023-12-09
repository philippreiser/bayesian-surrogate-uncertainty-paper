library(cubature)
library(mvtnorm)

get_design_matrix <- function(w, model="slope_cubic"){
  if(model == "slope_cubic"){
  return(matrix(c(w,
           w^3),
         nrow=length(w), ncol=2))
  } else if(model == "intercept_slope"){
    return(matrix(c(rep(1, length(w)),
      w),
                nrow=length(w), ncol=2))
  } 
}

get_y_model <- function(omega_t, c_gt, model){
  if (model == "slope"){
    y_t <- omega_t * c_gt
  } else if(model == "quadratic"){
    y_t <- omega_t**2 * c_gt
  } else if(model == "cubic"){
    y_t <- omega_t**3 * c_gt
  } else if(model == "slope_cubic" || model == "intercept_slope"){
    W_sim <- get_design_matrix(omega_t, model)
    y_t <- W_sim %*% c_gt
  } else if(model == "intercept"){
    y_t <- omega_t + c_gt
  }
  y_t
}

## s-step ##
get_t_posterior_analytic_parameters <- function(omega_t, y_t, mu_t0, sigma_t0,
                                          sigma_a, model){
  if (model == "quadratic"){
    omega_t <- omega_t^2
  } else if (model == "cubic"){
    omega_t <- omega_t^3
  } else if (model == "slope_cubic" || model == "intercept_slope"){
    W_sim <- get_design_matrix(omega_t, model)
    Sigma_t1 <- solve(sigma_t0^(-2)*diag(2)+sigma_a^(-2) * t(W_sim)%*%W_sim)
    mu_t1 <- Sigma_t1%*%(sigma_t0^(-2)*diag(2)%*%mu_t0+sigma_a^(-2)*t(W_sim)%*%y_t)
    return(list(sigma_t1 = Sigma_t1, mu_t1 = mu_t1))
  } else if (model == "intercept"){
    sigma_t1 <- sqrt(1/(sigma_t0^(-2)+sigma_a^(-2)))
    mu_t1 <- sigma_t1^(2)*(sigma_t0^(-2)*mu_t0 + sigma_a^(-2)*(y_t-omega_t))
    return(list(sigma_t1 = sigma_t1, mu_t1 = mu_t1))
  }
  sigma_t1 <- sqrt(1/(sigma_t0^(-2)+sigma_a^(-2)*omega_t^2))
  mu_t1 <- sigma_t1^(2)*(sigma_t0^(-2)*mu_t0 + sigma_a^(-2)*omega_t * y_t)
  list(sigma_t1 = sigma_t1, mu_t1 = mu_t1)
}

get_t_posterior_mcmc <- function(omega_t, y_t, mu_t0, sigma_t0, sigma_a, 
                                      t_model_file, seed, adapt_delta=0.95){
  sdata <- list(
    N_sim = length(omega_t),
    omega_t = matrix(c(omega_t), ncol=1),
    y_t = c(y_t),
    prior_only = FALSE,
    sigma_a = sigma_a
  )
  smodel <- cmdstan_model(t_model_file)
  t_model_fit <- smodel$sample(
    data = sdata,
    seed = seed,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = 1000,
    iter_warmup = 1000,
    adapt_delta = adapt_delta
  )
  t_model_fit
}

## point i-step ##

get_idata <- function(mu_t1, sigma_t1, y_i, sigma_i, model){
  if (model == "slope" || model == "cubic" || model == "intercept"){
    c_samples_analytic <- matrix(rnorm(n=4000, mean=mu_t1, sd=sigma_t1))
  }
  else if (model == "slope_cubic" || model == "intercept_slope"){
    c_samples_analytic <- rmvnorm(n=4000, mean=mu_t1, sigma=sigma_t1)
  }
  idata_st = list(
    y_i = matrix(y_i),
    N_i = 1,
    N_measures = 1,
    c = c_samples_analytic,
    cluster_weights = rep(1, nrow(c_samples_analytic))/nrow(c_samples_analytic),
    N_clusters = nrow(c_samples_analytic),
    prior_only = FALSE,
    sigma_i = sigma_i
  )
}

get_point_i_posterior_analytic <- function(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, model){
  if (model == "quadratic"){
    mu_t1 <- mu_t1^2
  } else if (model == "cubic"){
    mu_t1 <- mu_t1^3
  } else if (model == "intercept"){
    sigma_i1_point <- sqrt(1/(sigma_i0^(-2)+sigma_i^(-2)))
    mu_i1_point <- sigma_i1_point^(2)*(sigma_i0^(-2)*mu_i0 + sigma_i^(-2)*(y_i-mu_t1))
    return(list(sigma_i1_point=sigma_i1_point, mu_i1_point=mu_i1_point))
  } else if (model == "intercept_slope"){
    sigma_i1_point <- sqrt(1/(sigma_i0^(-2)+sigma_i^(-2)*mu_t1[2]^2))
    mu_i1_point <- sigma_i1_point^(2)*(sigma_i0^(-2)*mu_i0 + sigma_i^(-2)*mu_t1[2]*(y_i-mu_t1[1]))
    return(list(sigma_i1_point=sigma_i1_point, mu_i1_point=mu_i1_point))
  }
  sigma_i1_point <- sqrt(1/(sigma_i0^(-2)+sigma_i^(-2)*mu_t1^2))
  mu_i1_point <- sigma_i1_point^(2)*(sigma_i0^(-2)*mu_i0 + sigma_i^(-2)*mu_t1 * y_i)
  list(sigma_i1_point=sigma_i1_point, mu_i1_point=mu_i1_point)
}

get_point_i_posterior_numeric <- function(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, model, omega_i_range){
  post_omega_i <- dnorm(omega_i_range, mu_i0, sigma_i0) * dnorm(y_i, get_y_model(omega_i_range, mu_t1, model), sigma_i)
  post_omega_i <- post_omega_i / sum(post_omega_i)/(omega_i_range[2]-omega_i_range[1])
  post_omega_i
}

get_i_step_point_flat_prior_numeric <- function(y_i, sigma_i, mu_t1, model, omega_i_range){
  post_omega_i <- dnorm(y_i, get_y_model(omega_i_range, mu_t1, model), sigma_i)
  post_omega_i <- post_omega_i / sum(post_omega_i)
  post_omega_i
}

get_point_i_posterior_mcmc <- function(y_i, mu_i0, sigma_i0, sigma_i, mu_t1,
                                 imodel_file, seed, adapt_delta){
  idata = list(
    y_i = matrix(y_i),
    N_i = length(y_i),
    N_measures = 1,
    c = c(mu_t1),
    prior_only = FALSE,
    sigma_i = sigma_i
  )
  imodel <- cmdstan_model(imodel_file)
  imodel_fit <- imodel$sample(
    data = idata,
    seed = seed,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = 1000,
    iter_warmup = 1000,
    adapt_delta = adapt_delta
  )
  imodel_fit
}

## i-step ml ##

get_e_lik_i_posterior_mcmc <- function(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, sigma_t1,
                               i_model_file_e_lik, seed, model,
                               adapt_delta=0.95){
  idata <- get_idata(mu_t1, sigma_t1, y_i, sigma_i, model)
  imodel_st <- cmdstan_model(i_model_file_e_lik)
  imodel_st_fit <- imodel_st$sample(
    data = idata,
    seed = seed,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = 2000,
    iter_warmup = 1000,
    adapt_delta = adapt_delta,
    init = function() list(
      omega_i=-0.5)
  )
  imodel_st_fit
}

get_e_lik_i_posterior_numeric <- function(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, sigma_t1, model, i_prior_densities, omega_i_range,
                                  numint_lower_c, numint_upper_c,
                                  numint_lower_w, numint_upper_w,
                                  r_tol=1e-5){
  if (model == "slope" || model == "cubic" || model == "intercept"){
    integrand <- function(c, w_i) {dnorm(y_i, get_y_model(w_i, c, model), sigma_i)*dnorm(c, mean=mu_t1, sd=sigma_t1)}
    norm_integrand <- function(x) {
      dnorm(y_i, get_y_model(x[1], c(x[2]), model), sigma_i)*dnorm(c(x[2]), mu_t1, sigma_t1)*dnorm(x[1], mu_i0, sigma_i0)
    }
  }
  else if (model == "slope_cubic"|| model == "intercept_slope"){
    integrand <- function(c, w_i) {dnorm(y_i, get_y_model(w_i, c, model), sigma_i)*dmvnorm(c, mu_t1, sigma_t1)}
    # x = [w_i, c_1, c_2]
    norm_integrand <- function(x) {
      dnorm(y_i, get_y_model(x[1], c(x[2], x[3]), model), sigma_i)*dmvnorm(c(x[2], x[3]), mu_t1, sigma_t1)*dnorm(x[1], mu_i0, sigma_i0)
    }
  }
  integration_e_lik <- function(w_i) {cubintegrate(integrand, lower = numint_lower_c, upper = numint_upper_c, method = "cuhre", w_i=w_i, relTol=r_tol)$integral}
  i_likelihood_densities_list <- sapply(omega_i_range, integration_e_lik)
  # Uncomment to get proper normalization by numberical integration
  # normalization <- cubintegrate(norm_integrand, c(numint_lower_w, numint_lower_c), c(numint_upper_w, numint_upper_c), method = "cuhre")$integral
  # i_posterior_densities <-  i_prior_densities * i_likelihood_densities_list / normalization
  i_posterior_densities <-  i_prior_densities * i_likelihood_densities_list
  
  i_posterior_densities/sum(i_posterior_densities)/(omega_i_range[2]-omega_i_range[1])
}



## i-step mp ##

get_e_post_i_posterior_mcmc <- function(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, sigma_t1,
                               i_model_file_point, seed, model, n_c_skip=5,
                               adapt_delta=0.95){
  imodel <- cmdstan_model(i_model_file_point)
  idata <- get_idata(mu_t1, sigma_t1, y_i, sigma_i, model)
  fits_csv_files <- c()
  csv_dir <- file.path(here(), "fitted_models/_imodel_fits_cmdstan")
  dir.create(csv_dir, showWarnings = FALSE)
  # for (i in (1:nrow(idata_st$c))){
  for (i in (seq(from = 1, to = 1000, by = n_c_skip))){
    idata_tmp <- idata
    idata_tmp$c <- c(idata$c[i, ])
    idata_tmp$c_0 <- idata$c_0[[i]]
    imodel_fit_tmp <- imodel$sample(
      data = idata_tmp,
      seed = seed,
      chains = 4,
      parallel_chains = 4,
      iter_sampling = 1000,
      iter_warmup = 1000,
      output_dir=csv_dir,
      adapt_delta = adapt_delta
    )
    fits_csv_files <- append(fits_csv_files, imodel_fit_tmp$output_files())
  }
  imodel_mt_fit_data <- read_cmdstan_csv(fits_csv_files)
  imodel_mt_fit <- imodel_mt_fit_data$post_warmup_draws
  unlink(csv_dir, recursive = TRUE)
  imodel_mt_fit
}

get_e_post_i_posterior_numeric <- function(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, sigma_t1, model, i_prior_densities, omega_i_range,
                                  numint_lower_c, numint_upper_c,
                                  numint_lower_w, numint_upper_w,
                                  r_tol=1e-5){
  i_prior_densities <- dnorm(omega_i_range, mean=mu_i0, sd=sigma_i0)
  norm_integrand_e_post <- function(w, c) {dnorm(y_i, get_y_model(w, c, model), sigma_i)*dnorm(w, mu_i0, sigma_i0)}
  if (model == "slope"){
    normalization_e_post <- function(c) dnorm(y_i, c*mu_i0, sqrt(c^2*sigma_i0^2+sigma_i^2))
  } else if (model == "intercept_slope"){
    normalization_e_post <- function(c) dnorm(y_i, c%*%t(get_design_matrix(mu_i0, model)), sqrt(t(c)%*%c*sigma_i0^2+sigma_i^2))
  } else{
    normalization_e_post <- function(c) {cubintegrate(norm_integrand_e_post,
                                                  lower=numint_lower_w,
                                                  upper=numint_upper_w,
                                                  method = "cuhre",
                                                  c=c, relTol=r_tol)$integral}
  }
  if (model == "slope" || model == "cubic" || model=="intercept"){
    integrand_e_post <- function(c, w_i) {
      dnorm(y_i, get_y_model(w_i, c, model), sigma_i)*dnorm(c, mu_t1, sigma_t1)/normalization_e_post(c)
    }
  }
  else if (model == "slope_cubic" || model == "intercept_slope"){
    integrand_e_post <- function(c, w_i) {dnorm(y_i, get_y_model(w_i, c, model), sigma_i)*dmvnorm(c, mu_t1, sigma_t1)/normalization_e_post(c)}
  }
  integration_e_post <- function(w_i) {cubintegrate(integrand_e_post,
                                                lower = numint_lower_c,
                                                upper = numint_upper_c,
                                                method = "cuhre", w_i=w_i,
                                                relTol=r_tol)$integral}
  # integration_e_post <- function(w_i) {integrate(integrand_e_post,
  #                                            lower = numint_lower_c,
  #                                            upper = numint_upper_c,
  #                                            w_i=w_i)$value}
  likelihood_e_post <- sapply(omega_i_range, integration_e_post)
  posterior_e_post <- likelihood_e_post * i_prior_densities
  posterior_e_post
}

## i-step mll ##
get_e_log_lik_i_posterior_analytic <- function(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, sigma_t1, model, i_prior_densities, omega_i_range){
  ## Prior ##
  i_prior_densities <- dnorm(omega_i_range, mean=mu_i0, sd=sigma_i0)
  if (model == "quadratic"){
    mu_t1 <- mu_t1^2
  } else if (model == "cubic"){
    mu_t1 <- mu_t1^3
  }
  # outcome of analytic derivation of E-Log-Lik with Normal Conjugate Model
  if (model == "slope"){
    sigma_i1 <-sqrt(1/(sigma_i0^(-2) + sigma_i^(-2)*(mu_t1^2 + sigma_t1^2)))
    mu_i1 <- sigma_i1^2*(sigma_i0^(-2)*mu_i0 + sigma_i^(-2)*mu_t1*y_i)
  } else if (model == "intercept"){
    sigma_i1 <- sqrt(1/(sigma_i0^(-2)+sigma_i^(-2)))
    mu_i1 <- sigma_i1^(2)*(sigma_i0^(-2)*mu_i0 + sigma_i^(-2)*(y_i-mu_t1))
  } else if (model == "intercept_slope"){
    sigma_i1 <- sqrt(1/(sigma_i0^(-2) + sigma_i^(-2)*(mu_t1[2]^2 + sigma_t1[2,2])))
    mu_i1 <- sigma_i1^2*(sigma_i0^(-2)*mu_i0 + sigma_i^(-2)*(mu_t1[2]*y_i-sigma_t1[1,2]-mu_t1[1]*mu_t1[2]))
  }
  
  i_posterior_densities <- dnorm(omega_i_range, mu_i1, sigma_i1)
}

get_e_log_lik_i_posterior_mcmc <- function(y_i, mu_i0, sigma_i0, sigma_i, mu_t1, sigma_t1,
                                i_model_file_e_log_lik, seed, model, adapt_delta=0.95){
  idata <- get_idata(mu_t1, sigma_t1, y_i, sigma_i, model)
  imodel_e_log_lik <- cmdstan_model(i_model_file_e_log_lik)
  imodel_e_log_lik_fit <- imodel_e_log_lik$sample(
    data = idata,
    seed = seed,
    chains = 4,
    parallel_chains = 4,
    iter_sampling = 1000,
    iter_warmup = 1000,
    adapt_delta = adapt_delta
  )
}

