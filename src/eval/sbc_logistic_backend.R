# library(SBC)
library(cmdstanr)
library(here)
source(file.path(here(), "src/utils/utils.R"))

combine_args <- function(args1, args2) {
  if(is.null(names(args1)) || is.null(names(args2))) {
    c(args1, args2)
  } else {
    shared <- intersect(names(args1), names(args2))
    shared <- setdiff(shared, "")
    for(s in shared) {
      args1[[s]] <- args2[[s]]
    }
    c(args1, args2[!(names(args2) %in% shared)])
  }
}

SBC_backend_multitrial <- function(model, ...) {
  stopifnot(inherits(model, "CmdStanModel"))
  if(length(model$exe_file()) == 0) {
    stop("The model has to be already compiled, call $compile() first.")
  }
  args <- list(...)
  unacceptable_params <- c("data", "parallel_chains ", "cores", "num_cores")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }
  structure(list(model = model, args = args), class = "SBC_backend_multitrial")
}

SBC_fit.SBC_backend_multitrial <- function(backend, generated, cores) {
  csv_dir <- file.path(here(), "fitted_models", strsplit(tempdir(), "/")[[1]][3])
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  fits_csv_files <- c()
  for (i in (1:nrow(generated$c))){
    generated_tmp <- generated
    generated_tmp$c <- c(generated_tmp$c[i, ])
    generated_tmp$c_0 <- generated_tmp$c_0[[i]]
    generated_tmp$sigma_sim <- generated_tmp$sigma_sim[[i]]
    fit <- do.call(backend$model$sample,
                   combine_args(backend$args,
                                list(
                                  data = generated_tmp,
                                  parallel_chains = cores,
                                  output_dir=csv_dir,
                                  init = function() list(
                                    w_exp=generated$w_exp_init)
                                )))
    fits_csv_files <- append(fits_csv_files, fit$output_files())
  }
  fit_all <- read_cmdstan_csv(fits_csv_files)
  unlink(csv_dir, recursive = TRUE)
  if (length(unique(generated$cluster_weights)) == 1){
    return(fit_all$post_warmup_draws)
  }
  else {
    fit_all$cluster_weights <- generated$cluster_weights
    return(fit_all)
  }
}

SBC_fit_to_draws_matrix.draws_array<- function(fit) {
  as_draws_matrix(fit)
}

SBC_fit_to_draws_matrix.list<- function(fit) {
  draws_w_exp <- subset_draws(fit$post_warmup_draws, variable=c("w_exp[1,1]", "sigma_exp"))
  chains <- max(fit$metadata$id)
  iter_sampling <- fit$metadata$iter_sampling
  cluster_weights <- matrix(rep(fit$cluster_weights, each=chains*iter_sampling), ncol=ncol(draws_w_exp))
  cluster_weights_prob <- cluster_weights/sum(cluster_weights)
  draws_w_exp_weighted <- resample_draws(draws_w_exp, weights = cluster_weights_prob, ndraws=chains*iter_sampling)
  as_draws_matrix(draws_w_exp_weighted)
}

## Single-trial
SBC_backend_singletrial <- function(model, ...) {
  stopifnot(inherits(model, "CmdStanModel"))
  if(length(model$exe_file()) == 0) {
    stop("The model has to be already compiled, call $compile() first.")
  }
  args <- list(...)
  unacceptable_params <- c("data", "parallel_chains ", "cores", "num_cores")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }
  structure(list(model = model, args = args), class = "SBC_backend_singletrial")
}

SBC_fit.SBC_backend_singletrial <- function(backend, generated, cores) {
  fit <- do.call(backend$model$sample,
                 combine_args(backend$args,
                              list(
                                data = generated,
                                parallel_chains = cores,
                                init = function() list(
                                  w_exp=generated$w_exp_init)
                              )))
  fit
}

## Mean
SBC_backend_mean <- function(model, ...) {
  stopifnot(inherits(model, "CmdStanModel"))
  if(length(model$exe_file()) == 0) {
    stop("The model has to be already compiled, call $compile() first.")
  }
  args <- list(...)
  unacceptable_params <- c("data", "parallel_chains ", "cores", "num_cores")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }
  structure(list(model = model, args = args), class = "SBC_backend_mean")
}

SBC_fit.SBC_backend_mean <- function(backend, generated, cores) {
  generated$c <- generated$c_mean
  generated$c_0 <- generated$c_0_mean
  generated$sigma_sim <- generated$sigma_sim_mean
  fit <- do.call(backend$model$sample,
                 combine_args(backend$args,
                              list(
                                data = generated,
                                parallel_chains = cores,
                                init = function() list(
                                  w_exp=generated$w_exp_init)
                              )))
  fit
}

## Median
SBC_backend_median <- function(model, ...) {
  stopifnot(inherits(model, "CmdStanModel"))
  if(length(model$exe_file()) == 0) {
    stop("The model has to be already compiled, call $compile() first.")
  }
  args <- list(...)
  unacceptable_params <- c("data", "parallel_chains ", "cores", "num_cores")
  if(any(names(args) %in% unacceptable_params)) {
    stop(paste0("Parameters ", paste0("'", unacceptable_params, "'", collapse = ", "),
                " cannot be provided when defining a backend as they need to be set ",
                "by the SBC package"))
  }
  structure(list(model = model, args = args), class = "SBC_backend_median")
}

SBC_fit.SBC_backend_median <- function(backend, generated, cores) {
  generated$c <- generated$c_median
  generated$c_0 <- generated$c_0_median
  generated$sigma_sim <- generated$sigma_sim_median
  fit <- do.call(backend$model$sample,
                 combine_args(backend$args,
                              list(
                                data = generated,
                                parallel_chains = cores,
                                init = function() list(
                                  w_exp=generated$w_exp_init)
                              )))
  fit
}
