library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)
library(RcppDist)

# -------------------------------------------------------------------------

sourceCpp(code = '
  #include <RcppArmadillo.h>
  #include <RcppArmadilloExtensions/sample.h>
  #include <RcppEigen.h>
  #include <RcppDist.h>
  
  // [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppDist)]]
  
  using namespace Rcpp;
  using namespace arma;
  
  static bool exist(IntegerVector z, int j);
  static int K_z_j(IntegerVector z, int j);
  static NumericVector VN(int Kmax, int N, double gamma, double lambda = 1);
  static IntegerVector unique_z(IntegerVector x);
  static int num_k(IntegerVector x, int x_0);
  static double w_ldirich(NumericVector x, NumericVector eta, double w_x);
  static double lprior_s(int s, int n, double alpha_s, double beta_s);
  static double w_ldirich_full(NumericVector x, NumericVector eta, double w_x);
  static double logamma(double x, double shape, double rate);
  static double lbern(double x, double p);
  static List z_int(int K_max, int m);
  
  // [[Rcpp::export]]
  bool exist(IntegerVector z, int j) {
    return sum(z == z[j]) > 1;
  }
  
  // [[Rcpp::export]]
  int K_z_j(IntegerVector z, int j) {
    std::unordered_set<int> z_j;
    
    for (int jj = 0; jj < z.length(); jj++) {
      if (jj != j) {
        z_j.insert(z[jj]); 
      }
    }
    
    return z_j.size();
  }
  
  // https://arxiv.org/pdf/2312.08324
  // [[Rcpp::export]]
  NumericVector VN(int Kmax, int N, double gamma, double lambda) {
    int iters = std::max(N + 100, 1000);
    
    NumericVector vn(Kmax);
    for(int k = 1; k < Kmax + 1; k++) {
      double r = -datum::inf;
      
      for(int t = k; t < iters + 1; t++) {
        double b = 0;
        
        vec s1 = arma::linspace(t - k + 1, t, k);
        b += sum(log(s1));
        
        vec s2 = arma::linspace(t * gamma, t * gamma + N - 1, N);
        b += -sum(log(s2));
        
        double s3 = R::dpois(t - 1, lambda, false);
        b += s3;
        
        double m = std::max(b, r);
        r = log(exp(r - m) + exp(b - m)) + m;
      }
      vn(k - 1) = r;
    }
    return vn;
  }
  
  // [[Rcpp::export]]
  IntegerVector unique_z(IntegerVector x) {
    std::set<int> unique_vals(x.begin(), x.end());
    return IntegerVector(unique_vals.begin(), unique_vals.end());
  }

  // [[Rcpp::export]]
  int num_k(IntegerVector x, int x_0) {
    return sum(x == x_0);
  }
  
  // [[Rcpp::export]]
  double w_ldirich(NumericVector x, NumericVector eta, double w_x) {
    return w_x * sum((eta - 1) * log(x));
  }
  
  // [[Rcpp::export]]
  double lprior_s(int s, int n, double alpha_s, double beta_s) {
    double pt1 = -lgamma(s + 1) - lgamma(n - s + 1);
    double pt2 = lgamma(s + alpha_s) + lgamma(n - s + beta_s);
  
    return pt1 + pt2;
  }

  // [[Rcpp::export]]
  double w_ldirich_full(NumericVector x, NumericVector eta, double w_x) {
    double pt1 = sum(eta);
    double pt2 = lgamma(pt1);          
    
    double pt3 = sum(lgamma(eta));         

    double pt4 = sum((eta - 1) * log(x));  

    double pt5 = w_x * (pt2 - pt3 + pt4);

    return pt5;
  }

  // [[Rcpp::export]]
  double logamma(double x, double shape, double rate) {
    return (shape - 1) * log(x) - rate * x;
  }

  // [[Rcpp::export]]
  double lbern(double x, double p) {
    return x * log(p) + (1 - x) * log(1 - p);
  }

  // [[Rcpp::export]]
  List z_int(int K_max, int m) {
    IntegerVector z = sample(K_max, m, true) - 1;

    IntegerVector z_tmp = sort_unique(z);
    int K = z_tmp.size();

    z = match(z, z_tmp) - 1;

    return List::create(Named("z") = z, Named("K") = K);
  }  





  // [[Rcpp::export]]
  List bacon_mfm(NumericMatrix A, NumericMatrix L, double w_A, double w_L, bool doub_dirich, bool ddirch_A, bool ddirch_L, int Kmax_0, bool est_sr, bool est_s, bool est_r, double alpha_s, double beta_s, int iter) {
    NumericMatrix& A_L = ddirch_A ? A : L;
    int m = A_L.nrow();
    int n = A_L.ncol();
    
    int K_max = std::min(Kmax_0, m);
    
    double a_alpha = 0.001;
    double b_alpha = 0.001;
    
    double omega = 0.5;
    
    double tau_theta = 1;
    double tau_lambda = 1;
    double a_theta_tau = 0.1;
    double b_theta_tau = 0.1;
    double a_lambda_tau = 100;
    double b_lambda_tau = 100;
    double a_theta = 0.001;  
    double b_theta = 0.001;
    double a_lambda = 0.001;  
    double b_lambda = 0.001; 
    
    
    
    
    
    int i, ii, j, k, t, prog = 10;
    double accept_z = 0, accept_alpha = 0, accept_theta = 0, accept_lambda = 0;
    double burn_theta = 0, burn_lambda = 0;
    int ind_K = -1, s_0, r_0, sr_0, ind_sr, ind_prob_sr, sr_new, s_new, r_new, ind_s, ind_r, ind_new;
    double theta_ki_tmp, theta_ki_new, lambda_ki_tmp, lambda_ki_new, r_theta = 0, r_lambda = 0;
    double post_prob_s_tmp = 0, post_prob_r_tmp = 0, post_prob_theta_tmp = 0, post_prob_lambda_tmp = 0;
    
    
    
    List z_0 = z_int(K_max, m);
    IntegerVector z = z_0["z"];    
    int K = as<int>(z_0["K"]);
    
    
    double alpha = 1;
    

    // Create the spaces to store the results 
    IntegerMatrix z_store(iter, m);
    NumericVector K_store(iter * 2 * m);
    NumericVector alpha_store(iter);
    IntegerMatrix s_store(iter, m);
    IntegerMatrix r_store(iter, m);
    List Theta_store(iter);
    List Lambda_store(iter);
    NumericVector map_z_store(iter);
    NumericVector map_s_store(iter);
    NumericVector map_r_store(iter);
   
    
    
    

  
    
    
    
    
    
    
    
    NumericMatrix A_1(m, n);
    NumericMatrix L_1(m, n);
    A_1 = A;
    L_1 = L;
    
    
    
    
    
    
    
    
    
    
    
    
    
   
    


    
    IntegerVector int_sr(n * 2);
    for (sr_0 = 0; sr_0 < n * 2; sr_0++) {
      int_sr(sr_0) = sr_0;
    }
    NumericVector prob_sr(n * 2);
    
    
    IntegerVector sr_int(n * 2);
    for (sr_0 = 0; sr_0 < n * 2; sr_0++) {
      sr_int(sr_0) = sr_0;
    }
    
    

    IntegerVector s(m);
    IntegerVector s_int = seq(0, n - 1);
    NumericVector prob_s(n);
    
    
    IntegerVector r(m);
    NumericVector r_int = {0, 1};
    NumericVector prob_r(2);
    
    
    
    NumericMatrix Theta(K, n);
    NumericMatrix Lambda(K, n);
    std::fill(Theta.begin(), Theta.end(), 1.0);
    std::fill(Lambda.begin(), Lambda.end(), 1.0);
    
    
    
   
    NumericVector A1_new(n);
    NumericVector L1_new(n);
    
    
    

    
    NumericVector Theta_new(n);
    NumericVector Lambda_new(n);
    
    
    

    NumericVector post_prob_s(iter);
    NumericVector post_prob_r(iter);
    NumericVector post_prob_theta(iter);
    NumericVector post_prob_lambda(iter);
    
    
    NumericVector v_m = VN(K_max + 1, m, 1);
    
    
    
    
    
    
    
    int burn = static_cast<int>(std::floor(iter / 2));
    
    
    for (t = 0; t < iter; t++) {
      Rcpp::Rcout << "t = " << t << endl;
      
      double temp = pow(0.999, t + 1);
      
      
      for (j = 0; j < m; j++) {
        if (exist(z, j)) {
          int K_new = K + 1;
          if (K_new > K_max) {
            continue;
          }
          int zj_new = K_new;
        
          int K_zj = K_z_j(z, j);
          int zj_tmp = z(j);
        
          NumericVector theta_tmp = Rcpp::rgamma(n, a_theta, 1/b_theta);
          NumericVector lambda_tmp = Rcpp::rgamma(n, a_lambda, 1/b_lambda);
        
          double r_z0 = log(alpha) - log(m - 1 + K_zj * alpha);
        
          r_z0 += w_ldirich_full(A_1(j, _), theta_tmp, w_A) - w_ldirich_full(A_1(j, _), Theta(zj_tmp, _), w_A);
          if (doub_dirich) {
            r_z0 += w_ldirich_full(L_1(j, _), lambda_tmp, w_L) - w_ldirich_full(L_1(j, _), Lambda(zj_tmp, _), w_L);
          }
        
          double r_z = std::max(v_m(K_zj) - v_m(K_zj - 1) + r_z0, -999999.0);
        
          if (r_z >= log(R::runif(0, 1))) {
            K = K_new;
            z(j) = zj_new;
            
            NumericMatrix Theta_nummat(K, n);
            std::copy(Theta.begin(), Theta.end(), Theta_nummat.begin());
            Theta_nummat(K - 1, _) = theta_tmp;
            Theta = Theta_nummat;
            
            NumericMatrix Lambda_nummat(K, n);
            std::copy(Lambda.begin(), Lambda.end(), Lambda_nummat.begin());
            Lambda_nummat(K - 1, _) = lambda_tmp;
            Lambda = Lambda_nummat;
            
            if (t > burn) {
              accept_z = accept_z + 1;
            }
          }
        }
        else {
          IntegerVector z_tmp = z;
          z_tmp.erase(z_tmp.begin() + j);
          IntegerVector znew_int = unique_z(z_tmp);
          int Knew_int = znew_int.size();
          
          NumericVector prob_z(Knew_int);
          for (k = 0; k < Knew_int; k++) {
            prob_z(k) = std::max(log(num_k(z_tmp, znew_int(k)) + alpha) - log(m - 1 + Knew_int * alpha), -999999.0);
          }
          
          prob_z = (prob_z - max(prob_z)) / temp;
          prob_z = exp(prob_z) / sum(exp(prob_z));
          
          int zj_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(znew_int, 1, TRUE, prob_z));
          int zj_tmp = z(j);
          
          double r_z0 = log(m - 1 + Knew_int * alpha) - log(alpha);
          
          r_z0 += w_ldirich_full(A_1(j, _), Theta(zj_new, _), w_A) - w_ldirich_full(A_1(j, _), Theta(zj_tmp, _), w_A);
          if (doub_dirich) {
            r_z0 += w_ldirich_full(L_1(j, _), Lambda(zj_new, _), w_L) - w_ldirich_full(L_1(j, _), Lambda(zj_tmp, _), w_L);
          }
          
          double r_z = std::max(v_m(Knew_int - 1) - v_m(Knew_int) + r_z0, -999999.0);
          
          if (r_z >= log(R::runif(0, 1))) {
            z(j) = zj_new;
            z = ifelse(z > zj_tmp, z - 1, z);
            
            arma::mat theta_mat(Theta.begin(), K, n, false);
            theta_mat.shed_row(zj_tmp);
            
            arma::mat lambda_mat(Lambda.begin(), K, n, false);
            lambda_mat.shed_row(zj_tmp);
            
            K -= 1;
            
            Theta = NumericMatrix(K, n);
            std::copy(theta_mat.begin(), theta_mat.end(), theta_mat.begin());
            
            Lambda = NumericMatrix(K, n);
            std::copy(lambda_mat.begin(), lambda_mat.end(), lambda_mat.begin());
            
            if (t > burn) {
              accept_z += 1;
            }
          }
        }
      
        ind_K += 1;
        K_store(ind_K) = K;
      }
      
      
      
      for (j = 0; j < m; j++) {
        if (exist(z, j)) {
          IntegerVector z_int = seq(0, K - 1);
          
          IntegerVector z_tmp = z;
          z_tmp.erase(z_tmp.begin() + j);
          
          NumericVector prob_k(K);
          for (k = 0; k < K; k++) {
            double prob_k0 = log(num_k(z_tmp, z_int(k)) + alpha) - log(m - 1 + K * alpha);
            
            prob_k0 += w_ldirich_full(A_1(j, _), Theta(k, _), w_A);
            if (doub_dirich) {
              prob_k0 += w_ldirich_full(L_1(j, _), Lambda(k, _), w_L);
            }
            
            prob_k(k) = std::max(prob_k0, -999999.0);
          }
          
          prob_k = (prob_k - max(prob_k)) / temp;
          prob_k = exp(prob_k) / sum(exp(prob_k));
          
          z(j) = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(z_int, 1, TRUE, prob_k));
        }
      
        ind_K += 1;
        K_store(ind_K) = K;
      } 
      
      
      
      double eta = R::rbeta(alpha + 1, m);

      double pi_eta = exp(log(a_alpha + K - 1) - log(a_alpha + m * (b_alpha - log(eta)) + K - 1));
  
      double u = R::runif(0, 1);
  
      if (u <= pi_eta) {
        alpha = R::rgamma(a_alpha + K, 1 / (b_alpha - log(eta)));
      }
      else {
        alpha = R::rgamma(a_alpha + K - 1, 1 / (b_alpha - log(eta)));
      }
      
      
      
      
      
      
      if (est_sr) {
        for (j = 1; j < m; j++) {
          for(r_0 = 0; r_0 < 2; r_0++) {
            for (s_0 = 0; s_0 < n; s_0++) {
              for (i = 0; i < n; i++) {
                for (ii = 0; ii < n; ii++) {
                  ind_sr = ((s_0 + i * (1 - 2 * r_0)) + n) % n;
                  if (ind_sr == ii) {
                    A1_new(i) = A(j, ii);
                    L1_new(i) = L(j, ii);
                  }
                }
              }
              
              ind_prob_sr = s_0 + n * r_0;
              prob_sr(ind_prob_sr) = w_ldirich(A1_new, Theta(z(j), _), w_A);
              if (doub_dirich) {
                prob_sr(ind_prob_sr) += w_ldirich(L1_new, Lambda(z(j), _), w_L);
              }
              
              prob_sr(ind_prob_sr) += lprior_s(s_0, n - 1, alpha_s, beta_s);
              prob_sr(ind_prob_sr) += std::max(lbern(r_0, omega), -999999.0);
            }
          }
          
          prob_sr = (prob_sr - max(prob_sr)) / temp;
          prob_sr = exp(prob_sr) / sum(exp(prob_sr));
          
          sr_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(sr_int, 1, TRUE, prob_sr));
          
          if (sr_new < n) {
            r(j) = 0;
          }
          else {
            r(j) = 1;
          }
          
          s(j) = sr_new - n * r(j);
        }
      }
      
      
      
      
      
      
      
      
      if (est_s) {
        for (j = 1; j < m; j++) {
          for (s_0 = 0; s_0 < n; s_0++) {
            for (i = 0; i < n; i++) {
              for (ii = 0; ii < n; ii++) {
                ind_s = ((s_0 + i * (1 - 2 * r(j))) + n) % n;
                if (ind_s == ii) {
                  A1_new(i) = A(j, ii);
                  L1_new(i) = L(j, ii);
                }
              }
            }
            
            prob_s(s_0) = w_ldirich(A1_new, Theta(z(j), _), w_A);
            if (doub_dirich) {
              prob_s(s_0) += w_ldirich(L1_new, Lambda(z(j), _), w_L);
            }
            
            prob_s(s_0) += std::max(lprior_s(s_0, n - 1, alpha_s, beta_s), -999999.0);
          }
          
          prob_s = (prob_s - max(prob_s)) / temp;
          prob_s = exp(prob_s) / sum(exp(prob_s));
          
          s_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(s_int, 1, TRUE, prob_s));
          s(j) = s_new;
        }
      }
      
      
      
      
      
      
      
      
      if (est_r) {
        for (j = 1; j < m; j++) {
          for (r_0 = 0; r_0 < 2; r_0++) {
            for (i = 0; i < n; i++) {
              for (ii = 0; ii < n; ii++) {
                ind_r = ((s(j) + i * (1 - 2 * r_0)) + n) % n;
                if (ind_r == ii) {
                  A1_new(i) = A(j, ii);
                  L1_new(i) = L(j, ii);
                }
              }
            }
            
            prob_r(r_0) = w_ldirich(A1_new, Theta(z(j), _), w_A);
            if (doub_dirich) {
              prob_r(r_0) += w_ldirich(L1_new, Lambda(z(j), _), w_L);
            }
            
            prob_r(r_0) += std::max(lbern(r_0, omega), -999999.0);
          }
          
          prob_r = (prob_r - max(prob_r)) / temp;
          prob_r = exp(prob_r) / sum(exp(prob_r));
          
          r_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(r_int, 1, TRUE, prob_r));
          r(j) = r_new;
        }
      }
      
      
      
      
      
      
      
      
      for (j = 1; j < m; j++) {
        for (i = 0; i < n; i++) {
          for (ii = 0; ii < n; ii++) {
            int ind_new = ((s(j) + i * (1 - 2 * r(j))) + n) % n;  
            if (ind_new == ii) {
              A_1(j, i) = A(j, ii);
              L_1(j, i) = L(j, ii);
            }
          }
        }
      }
      
      
      
      
      
      
      
      
      
      for (k = 0; k < K; k++) {
        Theta_new = Theta(k, _);
        
        for (i = 0; i < n; i++) {
          r_theta = 0;
        
          theta_ki_tmp = Theta(k, i);
          
          theta_ki_new = exp(r_truncnorm(log(theta_ki_tmp), tau_theta, log(a_theta_tau), log(b_theta_tau)));
          Theta_new(i) = theta_ki_new;
          
          for (j = 0; j < m; j++) {
            if (z(j) == k) {
              r_theta += w_ldirich_full(A_1(j, _), Theta_new, w_A);
              r_theta -= w_ldirich_full(A_1(j, _), Theta(k, _), w_A);
            }
          }
          
          r_theta += logamma(theta_ki_new, a_theta, b_theta);
          r_theta -= logamma(theta_ki_tmp, a_theta, b_theta);
          
          if (t > burn) {
              burn_theta += 1;
          }
          
          if (r_theta >= log(R::runif(0, 1))) {
            Theta(k, i) = theta_ki_new;
            
            if (t > burn) {
              accept_theta += 1;
            }
          }
        }
      }
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      if (doub_dirich) {
        for (k = 0; k < K; k++) {
          Lambda_new = Lambda(k, _);
          
          for (i = 0; i < n; i++) {
            r_lambda = 0;
          
            lambda_ki_tmp = Lambda(k, i);
            
            lambda_ki_new = exp(r_truncnorm(log(lambda_ki_tmp), tau_lambda, log(a_lambda_tau), log(b_lambda_tau)));
            Lambda_new(i) = lambda_ki_new;
            
            for (j = 0; j < m; j++) {
              if (z(j) == k) {
                r_lambda += w_ldirich_full(L_1(j, _), Lambda_new, w_L);
                r_lambda -= w_ldirich_full(L_1(j, _), Lambda(k, _), w_L);
              }
            }
            
            r_lambda += logamma(lambda_ki_new, a_lambda, b_lambda);
            r_lambda -= logamma(lambda_ki_tmp, a_lambda, b_lambda);
            
            if (t > burn) {
              burn_lambda += 1;
            }
            
            if (r_lambda >= log(R::runif(0, 1))) {
              Lambda(k, i) = lambda_ki_new;
              
              if (t > burn) {
                accept_lambda += 1;
              }
            }
          }
        } 
      }
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      post_prob_s_tmp = 0;
      post_prob_r_tmp = 0;
      for (j = 1; j < m; j++) {
        post_prob_s_tmp += w_ldirich(A_1(j, _), Theta(z(j), _), w_A);
        post_prob_r_tmp += w_ldirich(A_1(j, _), Theta(z(j), _), w_A);
        
        if (doub_dirich) {
          post_prob_s_tmp += w_ldirich(L_1(j, _), Lambda(z(j), _), w_L);
          post_prob_r_tmp += w_ldirich(L_1(j, _), Lambda(z(j), _), w_L);
        }
        
        post_prob_s_tmp += lprior_s(s(j), n - 1, alpha_s, beta_s);
        post_prob_r_tmp += lbern(r(j), omega);
      }
      
      
      
      
      
      post_prob_s(t) = post_prob_s_tmp;
      post_prob_r(t) = post_prob_r_tmp;
    
      
      
      
      post_prob_theta_tmp = 0;
      post_prob_lambda_tmp = 0;
      for (j = 0; j < m; j++) {
        post_prob_theta_tmp += w_ldirich_full(A_1(j, _), Theta(z(j), _), w_A);
        
        for (i = 0; i < n; i++) {
          post_prob_theta_tmp += logamma(Theta(z(j), i), a_theta, b_theta);
        }
        
        if (doub_dirich) {
          post_prob_lambda_tmp += w_ldirich_full(L_1(j, _), Lambda(z(j), _), w_L);
          
          for (i = 0; i < n; i++) {
            post_prob_lambda_tmp += logamma(Lambda(z(j), i), a_lambda, b_lambda);
          }
        }
      }
      
      
      
      post_prob_theta(t) = post_prob_theta_tmp;
      if (doub_dirich) {
        post_prob_lambda(t) = post_prob_lambda_tmp;
      }
      
      
      
      
      
      
      
      z_store(t, _) = z;
      alpha_store(t) = alpha;
      s_store(t, _) = s;
      r_store(t, _) = r;
      Theta_store[t] = Theta;
      Lambda_store[t] = Lambda;
      
      
      
      
      if(((t * 100)/(iter - 1)) == prog) {
        Rcout << prog << "%" << std::endl;
        prog += 10;
      }
      
      
      
    }
    
    
    accept_z = accept_z/(burn * 2 * m);
    accept_alpha = accept_alpha/burn;
    accept_theta = accept_theta/burn_theta;
    if (doub_dirich) {
        accept_lambda = accept_lambda/burn_lambda;
    }
    
    

    return List::create(Named("z_store") = z_store, Named("accept_z") = accept_z,
                        Named("K_store") = K_store,
                        Named("alpha_store") = alpha_store,
                        Named("s_store") = s_store, Named("post_prob_s") = post_prob_s,
                        Named("r_store") = r_store, Named("post_prob_r") = post_prob_r,
                        Named("Lambda_store") = Lambda_store, Named("post_prob_theta") = post_prob_theta, Named("accept_lambda") = accept_lambda,
                        Named("Theta_store") = Theta_store, Named("post_prob_lambda") = post_prob_lambda, Named("accept_theta") = accept_theta);
  
    
    
    
  }



')
