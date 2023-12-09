library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(purrr)
library(orthopolynom)
library(SobolSequence)
library(here)
Rcpp::sourceCpp(file.path(here(),"src", "utils","PCE_helpers.cpp"))

# code adapted from Paul-Christian Bürkner
# https://github.com/paul-buerkner/Bayesian-sparse-PCE
get_pce_vars <- function(p, M, idx=NULL){
  #idx <- c(38, 66, 51, 166, 68, 168, 21, 41, 61, 70, 60, 31, 232, 170, 230, 75,
  #         234, 65, 33, 39, 6, 266, 35, 100, 102) # Indices from bayes sparse PCE
  l_polys <- legendre.polynomials(p, normalized=TRUE)
  l_poly_coeffs <- polynomial.coefficients(l_polys)
  l_poly_coeffs_mat <- matrix(0, p+1, p+1)
  for (i in 1:(p+1)){
    for (j in 1:length(l_poly_coeffs[[i]])){
      l_poly_coeffs_mat[i, j] = l_poly_coeffs[[i]][j]
    }
    if (i+1<=p+1){
      #l_poly_coeffs[[i]][seq(i+1,p+1)] <- 0
      l_poly_coeffs_mat[i, seq(i+1, p+1)] <- 0
    }
  }
  # Comb & Poly selection
  comb <- poly_idx_cpp(p, M)
  # first column is the constant polynomial f0
  comb <- comb[-1, , drop = FALSE]
  rownames(comb) <- seq_len(nrow(comb))
  if (!is.null(idx)) {
    # select only desired polynomials
    comb <- comb[idx, , drop = FALSE]
  }
  list(l_poly_coeffs_mat, comb)
}
