library(truncnorm)
library(SobolSequence)
library(dplyr)
library(here)
source(file.path(here(),"src", "utils","helpers_pce.R"))
source(file.path(here(), "src/utils/utils.R"))

get_1d_sdata <- function(n_sim, poly_degree, M, fct, sigma_sim=0.001, prior_only=0,
                         tmp=NULL, fixed_w_sim=FALSE, pce_idx=NULL){
  pce_vars <- get_pce_vars(poly_degree, M, pce_idx)
  l_poly_coeffs_mat <- pce_vars[[1]]
  comb <- pce_vars[[2]]
  if (is.null(tmp) & !fixed_w_sim){
    # tmp <- data.frame(x1 = runif(n_sim, -1, 1)) %>%
    tmp <- data.frame(x1 = seq(-1, 1, length.out=n_sim)) %>%
      mutate(
        N_sim = n_sim,
        y_sim = fct(x1)+ rnorm(n_sim, sd=sigma_sim),
        w_sim_1 = x1,
      )
  } else if (fixed_w_sim){
    tmp <- data.frame(x1 = get_w_sims_logistic()[1:n_sim]) %>%
      mutate(
        N_sim = n_sim,
        y_sim = fct(x1)+ rnorm(n_sim, sd=sigma_sim),
        w_sim_1 = x1,
      )
  }
  
  sdata <- list(
    N_sim = n_sim,
    sigma_sim_lower = 0.0,
    d = poly_degree,
    M = M,
    w_sim = matrix(c(tmp$w_sim_1), ncol=M),
    y_sim = tmp$y_sim,
    l_poly_coeffs = t(l_poly_coeffs_mat),
    comb = comb,
    N_comb = nrow(comb),
    prior_only = prior_only
  )
  sdata
}

get_training_data_indices <- function(n_all, n_query){
  # adapted python code from 
  # https://stackoverflow.com/questions/69790561/what-is-the-name-of-this-iteration-criterion-and-is-it-implemented-for-python
  if (n_all < n_query) {
    stop("n_all must be greater or equal n_query")
  }
  lst <- c(1:n_all)
  result <- c(lst[1], lst[length(lst)])
  ranges <- list(c(1, length(lst)-1))
  
  idx <- 1
  while (length(result)<n_query){
    start <- ranges[[idx]][[1]]
    stop <- ranges[[idx]][2]
    if (start < stop){
      middle <- floor((start + stop)/2)
      result[length(result)+1] <- lst[middle+1]
      ranges[[length(ranges)+1]] <- c(start, middle)
      ranges[[length(ranges)+1]] <- c(middle+1, stop)
    }
    idx <- idx + 1
  }
  return(result)
}

get_1d_sdatasets <- function(n_sims, poly_degree, M, fct, sigma_sim){
  sdatasets <- list()
  sdata <- get_1d_sdata(129, poly_degree, M, fct, sigma_sim, prior_only = 0,
                        fixed_w_sim = FALSE)
  for (n_sim in n_sims){
    indices <- get_training_data_indices(129, n_sim)
    sdata_tmp <- sdata
    sdata_tmp$N_sim <- n_sim
    sdata_tmp$y_sim <- sdata$y_sim[indices]
    sdata_tmp$w_sim <- matrix(sdata$w_sim[indices], ncol=M)
    sdatasets[[length(sdatasets)+1]] <- sdata_tmp
  }
  sdatasets
}

get_1d_idata <- function(n_exp, w_1_value, fct, sigma_exp=0.01){
  tmp_e <- data.frame(V1=rep(w_1_value, n_exp)) %>%
    rename(x1 = V1) %>%
    mutate(
      y_exp = fct(x1)+rnorm(n_exp, sd=sigma_exp),
      w_exp_1 = x1
      )
  idata <- list(
    w_exp_1 = matrix(tmp_e$w_exp_1, ncol=n_exp),
    y_exp = matrix(tmp_e$y_exp, ncol=n_exp)
  )
  idata
}

censor_w_exp <- function(w_exp){
  mask_1 <- (w_exp < 0.001)
  w_exp[mask_1] <- 0.001
  mask_2 <- (w_exp > 0.999)
  w_exp[mask_2] <- 0.999
  w_exp
}

get_w_exp <- function(number, mean=0.5, sd=0.18, lower=NULL, upper=NULL){
  if ( (!is.null(lower)) & (!is.null(upper)) ) {
    w_exp_1 <- rtruncnorm(number, lower, upper, mean, sd)
    w_exp_2 <- rtruncnorm(number, lower, upper, mean, sd)
  }
  else {
    w_exp_1 <- rnorm(number, mean, sd)
    w_exp_2 <- rnorm(number, mean, sd)
  }
  data.frame(w_exp_1=w_exp_1, w_exp_2=w_exp_2)
}

get_w_exp_unif <- function(number, min=0.0, max=1.0){
  w_exp_1 <- runif(number, min, max)
  w_exp_2 <- runif(number, min, max)
  data.frame(w_exp_1=w_exp_1, w_exp_2=w_exp_2)
}

i_step_1d_generator_single <- function(N, poly_degree, M,
                                       smodel_fit, func,
                                       pce_idx=NULL,
                                       method="mean",
                                       number_draws=25,
                                       upper=1, lower=-1,
                                       fixed_w_exp=NULL,
                                       fixed_y_exp=NULL,
                                       sigma_exp=NULL,
                                       rate_sigma_exp=30,
                                       prior_only=0,
                                       kmeans.nstart=25,
                                       kmeans.iter.max=100){
  # method = ["mean", "median_",
  #           "weighted_cluster_draws",
  #           "all_draws"]
  # number_draws: number of random draws from mcmc chain or number of clusters
  pce_vars <- get_pce_vars(poly_degree, M, pce_idx)
  l_poly_coeffs_mat <- pce_vars[[1]]
  comb <- pce_vars[[2]]
  cluster_weights <- c(1)
  # If w_exp are provided, use them
  if (is.null(fixed_w_exp)){
    w_exp <- get_w_exp(1, mean=0, sd=0.5, lower=lower, upper=upper)$w_exp_1
  }
  else{
    w_exp <- fixed_w_exp
  }
  if (is.null(sigma_exp)){
    sigma_exp <- runif(1, min=0., max=0.05)
  }
  idata <- get_1d_idata(N, w_exp, func, sigma_exp)
  if (is.null(fixed_y_exp)){
    y_exp <- idata$y_exp
  }
  else{
    y_exp <- fixed_y_exp
  }
  # process joint posterior distribution (c, sigma_s) from S-Step
  c_median <- c(apply(smodel_fit$draws("c", format="matrix"), 2, median))
  c_0_median <- c(median(smodel_fit$draws("c_0", format="matrix")))
  sigma_sim_median <- c(apply(smodel_fit$draws("sigma_sim", format="matrix"), 2, median))
  c_mean <- c(apply(smodel_fit$draws("c", format="matrix"), 2, mean))
  c_0_mean <- c(mean(smodel_fit$draws("c_0", format="matrix")))
  sigma_sim_mean <- c(apply(smodel_fit$draws("sigma_sim", format="matrix"), 2, mean))
  if (method == "mean"){
    c <- c(apply(smodel_fit$draws("c", format="matrix"), 2, mean))
    c_0 <- c(mean(smodel_fit$draws("c_0", format="matrix")))
    sigma_sim <- c(apply(smodel_fit$draws("sigma_sim", format="matrix"), 2, mean))
  } else if (method == "median"){
    c <- c(apply(smodel_fit$draws("c", format="matrix"), 2, median))
    c_0 <- c(median(smodel_fit$draws("c_0", format="matrix")))
    sigma_sim <- c(apply(smodel_fit$draws("sigma_sim", format="matrix"), 2, median))
  } else if (method == "weighted_cluster_draws" || method == "single_trial"){
    cluster_draws <- get_cluster_draws(smodel_fit$draws("c_0", format="matrix"),
                                       smodel_fit$draws("c", format="matrix"),
                                       smodel_fit$draws("sigma_sim", format="matrix"),
                                       number_centers=number_draws,
                                       nstart=kmeans.nstart, iter.max=kmeans.iter.max)
    sigma_sim <- cluster_draws$sigma_sim
    c_0 <- cluster_draws$c_0
    c <- cluster_draws$c
    cluster_weights <- cluster_draws$cluster_weights
  } else if (method == "all_draws" || method == "single_trial_all"){
    c <- smodel_fit$draws("c", format="matrix")
    c_0 <- as.vector(smodel_fit$draws("c_0", format="matrix"))
    sigma_sim <- as.vector(smodel_fit$draws("sigma_sim", format="matrix"))
    cluster_weights <- rep(1, length(c_0))/length(c_0)
  } else if (method == "thin"){
    thin <- 40
    c <- thin_draws(smodel_fit$draws("c", format="matrix"), thin=thin)
    c_0 <- as.vector(thin_draws(smodel_fit$draws("c_0", format="matrix"), thin=thin))
    sigma_sim <- as.vector(thin_draws(smodel_fit$draws("sigma_sim", format="matrix"), thin=thin))
    cluster_weights <- rep(1, length(c_0))/length(c_0)
  }
  
  list(
    variables = list(
      sigma_exp = sigma_exp,
      w_exp = matrix(c(w_exp), ncol=M)
    ),
    generated = list(
      d = poly_degree,
      M = M,
      y_exp = y_exp,
      l_poly_coeffs = t(l_poly_coeffs_mat),
      comb = comb,
      N_comb = nrow(comb),
      N_exp = nrow(idata$w_exp_1),
      N_measures = ncol(idata$w_exp_1),
      c = c,
      c_0 = c_0,
      cluster_weights = cluster_weights,
      N_clusters = length(cluster_weights),
      prior_only = prior_only,
      w_exp_init = w_exp,
      sigma_sim = sigma_sim,
      c_median = c_median,
      c_0_median = c_0_median,
      sigma_sim_median = sigma_sim_median,
      c_mean = c_mean,
      c_0_mean = c_0_mean,
      sigma_sim_mean = sigma_sim_mean
    )
  )
}
