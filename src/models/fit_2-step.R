library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(purrr)
library(orthopolynom)
library(SobolSequence)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(here)
library(cluster)
source(file.path(here(),"src/utils/helpers_pce.R"))
source(file.path(here(),"src/utils/true_models_helpers.R"))
source(file.path(here(),"src/data_generation/get_data.R"))
source(file.path(here(), "src/plotting/plotting.R"))

get_smodel_fit <- function(sdata, file_smodel, fct, plot_pp=TRUE,
                           model_title="PCE", hide_legend=FALSE,
                           title=NULL, iter_sampling=1000, chains=4,
                           adapt_delta=0.8, linewidth=1,
                           propagate_sigma_a=FALSE
                           ){
  smodel <- cmdstan_model(file_smodel)
  smodel_fit <- smodel$sample(
    data = sdata,
    seed = 100,
    chains = chains,
    parallel_chains = chains,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    refresh=0
  )
  smodel_draws <- smodel_fit$draws("c", format="matrix")
  if (plot_pp == TRUE){
    pp <- plot_pp_1d(smodel, smodel_fit, fct, t(sdata$l_poly_coeffs), 
                     sdata$comb, n_test=100, model_title=model_title, 
                     n_sim=sdata$N_sim, prior_only=sdata$prior_only,
                     hide_legend=hide_legend, title=title, sdata=sdata,
                     linewidth=linewidth,
                     propagate_sigma_a=propagate_sigma_a)
  }
  list(smodel_fit, pp)
}

sample_draws_weighted_mp <- function(imodel_fit, idata, chains, iter_sampling,
                                     variables=c("w_exp[1,1]", "sigma_exp")){
  draws_w_exp <- subset_draws(imodel_fit$post_warmup_draws, variable=variables)
  cluster_weights <- matrix(rep(idata$cluster_weights, each=chains*iter_sampling), ncol=ncol(draws_w_exp))
  cluster_weights_prob <- cluster_weights/sum(cluster_weights)
  draws_w_exp_weighted <- resample_draws(draws_w_exp, weights = cluster_weights_prob, ndraws=chains*iter_sampling)
  return(draws_w_exp_weighted)
}

get_imodel_fit <- function(idata, file_imodel, smodel_fit, sdata, method="mean",
                           number_draws=25, refresh=NULL, chains=4,
                           iter_sampling=1000, iter_warmup=1000, seed=101,
                           nweigted_samples=10000, adapt_delta=0.8, init=NULL){
  # method = ["mean", "single_draws", "multi_draws" or "cluster_draws"]
  imodel <- cmdstan_model(file_imodel)
  if (method == "mean" | method == "single_trial" | method == "single_trial_all"){
    imodel_fit <- imodel$sample(
      data = idata,
      seed = seed,
      chains = chains,
      parallel_chains = chains,
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup,
      refresh = refresh,
      adapt_delta = adapt_delta,
      init = init
    )
  }
  else if (method == "multi_draws" | method == "cluster_draws" |
           method == "weighted_cluster_draws" | method == "kmeanspp_draws"){
    fits_csv_files <- c()
    csv_dir <- file.path(here(), "fitted_models/_imodel_fits_cmdstan", strsplit(tempdir(), "/")[[1]][3])
    dir.create(csv_dir, showWarnings = FALSE, recursive=TRUE)
    library(future)
    plan(multisession, gc=TRUE, workers=10)
    for (i in (1:nrow(idata$c))){
      idata_tmp <- idata
      idata_tmp$c <- c(idata$c[i, ])
      idata_tmp$c_0 <- idata$c_0[[i]]
      idata_tmp$sigma_sim <- idata$sigma_sim[[i]]
      imodel_fit_tmp <- imodel$sample(
        data = idata_tmp,
        seed = seed,
        chains = chains,
        parallel_chains = chains,
        iter_sampling = iter_sampling,
        iter_warmup = iter_warmup,
        output_dir=csv_dir,
        refresh = refresh,
        adapt_delta = adapt_delta,
        init = init
      )
      fits_csv_files <- append(fits_csv_files, imodel_fit_tmp$output_files())
    }
    imodel_fit <- read_cmdstan_csv(fits_csv_files)
    unlink(csv_dir, recursive = TRUE)
  }
  if (method == "weighted_cluster_draws"){
    draws_weighted <- sample_draws_weighted_mp(imodel_fit, idata, chains, iter_sampling,
                                         variables=c("w_exp[1,1]", "sigma_exp"))
    imodel_fit <- draws_weighted # TODO: convert to CmdStanMCMC object
  }
  imodel_fit
}


