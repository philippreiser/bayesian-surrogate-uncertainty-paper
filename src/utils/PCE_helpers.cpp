#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>
using namespace Rcpp;
using namespace R;

// Code from Paul-Christian Bürkner:
// https://github.com/paul-buerkner/Bayesian-sparse-PCE/blob/main/PCE_helpers.cpp
// generates multi-indices of multivariate polynomial base
// uses graded lexicographic ordering (P. 156, Sullivan)
// adapted from Matlab code of Ilja Kröker
// @param p maximal order of polynomials
// @param M number of input variables
// [[Rcpp::export]]
NumericMatrix poly_idx_cpp(int p, int M) {
  double P = R::choose(double ((p + M)*1.0), double (M*1.0));
  // NumericMatrix out = na_matrix(P, M);
  NumericMatrix out(P, M);
  NumericVector tA(M);
  int l = 1;
  int pmax = pow(p + 1, M);
  for (int i = 1; i < pmax + 1; i++) {
    int ri = i;
    for (int d = 0; d < M; d++) {
      int md = pow(p + 1, M - d - 1);
      int val = floor(ri / md);
      tA[d] = val;
      ri = ri - val * md;
    }
    if (sum(tA) <= p) {
      out(l, _) = tA;
      l = l + 1;
    }
  }
  return out;
}



