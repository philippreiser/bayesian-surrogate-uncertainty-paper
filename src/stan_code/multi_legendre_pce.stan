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
  int<lower=1> N_sim;           // Number of input/output simulation pairs
  int<lower=1> M;               // dimension of input variable
  vector[N_sim] y_sim;          // Output simulation variables
  int<lower=1> d;               // Degree of polynomials
  matrix[N_sim, M] w_sim;       // Input simulation variables

  real sigma_sim_lower;         // Lower bound for sigma_sim
  matrix[d+1, d+1] l_poly_coeffs; // coefficients of legendrepolynomials
  int<lower=1> N_comb; // Number of selected polynomials
  array[N_comb, M] int comb; // polynomial selection
}

transformed data {
  matrix[N_sim, N_comb] W_sim = get_PCE(w_sim, d, l_poly_coeffs, comb, N_comb);
}

parameters {
  real c_0;                                      // coefficient for constant "polynomial"
  vector[N_comb] c;                                   // coefficients of non-constant polynomials
  real<lower=sigma_sim_lower> sigma_sim;         // sigma for simulation values
}


model {
  // Prior model
  c_0 ~ normal(0, 5);
  c ~ normal(0, 5);
  sigma_sim ~ exponential(1);
  
  // Observational model
  y_sim ~ normal_id_glm(W_sim, c_0, c, sigma_sim);
}

generated quantities {
  vector[N_sim] mu_pred_sim;
  mu_pred_sim = W_sim*c + c_0;

  array[N_sim] real y_rep_sim;
  y_rep_sim = normal_rng(mu_pred_sim, sigma_sim);

  vector[N_sim] log_lik_sim;
  for (n in 1:N_sim) log_lik_sim[n] = normal_lpdf(y_sim[n] | mu_pred_sim[n], sigma_sim);
}

