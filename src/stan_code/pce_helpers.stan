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
