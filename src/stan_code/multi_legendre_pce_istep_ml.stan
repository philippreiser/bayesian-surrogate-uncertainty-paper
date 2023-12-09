//#include pce_helpers.stan

functions {
  int pow_(int a, int b){
    int out = 1;
    for (i in 1:b){
      out = a * out;
    }
    return out;
  }
  // generates multi-indices of multivariate polynomial base
  // uses graded lexicographic ordering (P. 156, Sullivan)
  // adapted from Matlab code of Ilja Kröker
  // and C++ code of Paul Bürkner
  // @param p maximal order of polynomials
  // @param M number of input variables
  array[,] int get_poly_idx(int p, int M){
    int P = choose(p+M, M); 
    array[P, M] int out = rep_array(0, P, M);
    array[M] int tA = rep_array(0, M);
    int l = 2;
    int pmax = pow_(p+1, M);
    for (i in 1:pmax){
      int ri = i;
      for (d in 1:M){
        int md = pow_(p + 1, M - d);
        int val = ri %/% md;
        tA[d] = val;
        ri = ri - val * md;
      }
      if (sum(tA) <= p){
        out[l, :] = tA;
        l = l+1;
      }
    }
    return out;
  }

  matrix get_poly_order_d(vector w, int d){
    matrix[rows(w), d+1] w_p = rep_matrix(1, rows(w), d+1);
    for (i in 1:d){
      w_p[:, i+1] = w_p[:, i] .* w;
    }
    return w_p;
  }

  matrix get_PCE(matrix w_sim, int d, matrix l_poly_coeffs, array[,] int comb, int N_comb){
    // TODO: poly (for legendre, hermite, ...), scale
    int N = rows(w_sim);
    int M = cols(w_sim);
    array[M] matrix[N, d+1] poly;
    for (m in 1:M){
      matrix[N, d+1] w_sim_i_p = get_poly_order_d(w_sim[:, m], d);
      poly[m] = w_sim_i_p * l_poly_coeffs;
    }
    matrix[N, N_comb] out = rep_matrix(1., N, N_comb);
    for (i in 1:N_comb){
      for (j in 1:M){
          out[:, i] = out[:, i] .* to_vector(poly[j,:, comb[i, j]+1]); 
      }
    }
    return out;
  }
}


data {
  int<lower=1> N_exp;          // Number of exp (measured) values
  int<lower=1> N_measures;      // Number of measurements per w_exp
  int<lower=1> N_clusters;
  int<lower=1> M;               // dimension of input variable
  array[N_exp] vector[N_measures] y_exp; // Output exp (measured) variables
  int<lower=1> d;               // Degree of polynomials

  matrix[d+1, d+1] l_poly_coeffs; // coefficients of legendrepolynomials
  int<lower=1> N_comb; // Number of selected polynomials
  array[N_comb, M] int comb; // polynomial selection
  vector[N_clusters] c_0;                                      // coefficient for constant "polynomial"
  matrix[N_clusters, N_comb] c;                                   // coefficients of non-constant polynomials
  vector[N_clusters] cluster_weights;
  int prior_only;  // should the likelihood be ignored?
}

parameters {
  real<lower=0, upper=0.05> sigma_exp;
  matrix<lower=-1, upper=1>[N_exp, M] w_exp;
}

transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += uniform_lpdf(sigma_exp | 0, 0.05);
  for (n in 1:N_exp){
    for (m in 1:M){
      lprior += normal_lpdf(w_exp[n, m] | 0, 0.5);
    }
  }
}

model {
  // Observational model  
  if (!prior_only) {
    real llik = 0;
    matrix[N_exp, N_comb] w_exp_pce = get_PCE(w_exp, d, l_poly_coeffs, comb, N_comb);
    for (n in 1:N_exp) {
      vector[N_clusters] cluster_lliks;
      for (k in 1:N_clusters) {
        vector[N_exp] mu_y_exp =  c_0[k] + w_exp_pce*to_vector(c[k, :]);
        cluster_lliks[k] = log(cluster_weights[k]) + normal_lpdf(to_vector(y_exp[n, :]) | mu_y_exp[n], sigma_exp);
      }
      llik += log_sum_exp(cluster_lliks);
    }
    target += llik;
  }
  // Prior model
  target += lprior;
}
