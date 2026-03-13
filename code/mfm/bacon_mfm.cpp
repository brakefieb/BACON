#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
#include <RcppDist.h>

// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppDist)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double w_ldirich(NumericVector x, NumericVector eta, double w_x) {
    return w_x * sum((eta - 1.0) * log(x));
}

// [[Rcpp::export]]
double w_ldirich_full(NumericVector x, NumericVector eta, double w_x) {
    int n = x.size();
    double sum_eta = 0.0;
    double sum_lgamma_eta = 0.0;
    double kernel = 0.0;

    for (int i = 0; i < n; i++) {
        sum_eta += eta[i];
        sum_lgamma_eta += std::lgamma(eta[i]);
        kernel += (eta[i] - 1.0) * log(x[i]);
    }

    double lgamma_sum_eta = std::lgamma(sum_eta);
    
    return w_x * (lgamma_sum_eta - sum_lgamma_eta + kernel);
}

// [[Rcpp::export]]
bool exist(IntegerVector z, int j) {
    int m = z.size();
  
    for (int j_idx = 0; j_idx < m; j_idx++) {
        if (j_idx != j && z[j_idx] == z[j]) {
            return true;
        }
    }
    return false;
}

// [[Rcpp::export]]
IntegerVector unique_z(IntegerVector z) {
    IntegerVector res = clone(z);
    
    std::sort(res.begin(), res.end());
    auto it = std::unique(res.begin(), res.end());
    res.erase(it, res.end());
    
    return res;
}

// [[Rcpp::export]]
int num_k(IntegerVector z, int k) {
    int m = z.size();
    int num = 0;
  
    for (int j = 0; j < m; j++) {
        if (z[j] == k) {
            num++;
        }
    }
  
    return num;
}

// [[Rcpp::export]]
NumericVector V_n(int K_max, int n, double gamma, double lambda) {
    int iter = std::max(n + 100, 1000);
    NumericVector v_n(K_max);

    for(int k = 1; k < K_max + 1; k++) {
        double r = -datum::inf;
      
        for(int t = k; t < iter + 1; t++) {
            double b = 0.0;
        
            arma::vec s_1 = arma::linspace(t - k + 1, t, k);
            b += sum(log(s_1));
        
            arma::vec s_2 = arma::linspace(t * gamma, t * gamma + n - 1, n);
            b -= sum(log(s_2));
        
            double s_3 = R::dpois(t - 1, lambda, true);
            b += s_3;
        
            double m = std::max(b, r);
            r = log(exp(r - m) + exp(b - m)) + m;
        }
        v_n(k - 1) = r;
    }

    return v_n;
}

// [[Rcpp::export]]
double lprior_s(int s, int n, double alpha_s, double beta_s) {
    double res = -std::lgamma(s + 1.0) - std::lgamma(n - s + 1.0);
    res += std::lgamma(s + alpha_s) + std::lgamma(n - s + beta_s);
  
    return res;
}

// [[Rcpp::export]]
double lbern(int x, double p) {
    return x * log(p) + (1.0 - x) * log(1.0 - p);
}

// [[Rcpp::export]]
double logamma(double x, double shape, double rate) {
    return (shape - 1.0) * log(x) - rate * x;
}

// [[Rcpp::export]]
Rcpp::List bacon_mfm(NumericMatrix A, NumericMatrix L, double w_A, double w_L, 
                     IntegerVector z, int K_max, double alpha, 
                     bool est_sr, bool est_s, bool est_r, double alpha_s, double beta_s, double omega, 
                     double a_theta, double b_theta, double a_lambda, double b_lambda, double tau_theta, double tau_lambda, 
                     int iter) {
    int m = (w_A > 0.0) ? A.nrow() : L.nrow();
    int n = (w_A > 0.0) ? A.ncol() : L.ncol();
    int K = max(z) + 1;

    int t, j, i, k;
    int j_idx, k_idx, s_0, r_0;
    int ind_z = -1, ind_K = -1;

    NumericVector v_m = V_n(K_max + 1, m, 1.0, 1.0);

    IntegerVector sr_int = seq(0, n * 2 - 1);
    IntegerVector s_int = seq(0, n - 1);
    IntegerVector r_int = {0, 1};
    
    IntegerVector s(m);
    IntegerVector r(m);
    for (j = 1; j < m; j++) {
        s(j) = floor(R::runif(0, n));
        r(j) = floor(R::runif(0, 2));
    }

    NumericMatrix A_1 = (w_A > 0.0) ? NumericMatrix(m, n) : NumericMatrix();
    NumericMatrix L_1 = (w_L > 0.0) ? NumericMatrix(m, n) : NumericMatrix();
    if (w_A > 0.0) { A_1(0, _) = A(0, _); }
    if (w_L > 0.0) { L_1(0, _) = L(0, _); }
    for (j = 1; j < m; j++) {
        for (i = 0; i < n; i++) {
            int ind_0 = ((s(j) + i * (1 - 2 * r(j))) + n) % n;
                
            if (w_A > 0.0) {
                A_1(j, i) = A(j, ind_0);
            }
            if (w_L > 0.0) {
                L_1(j, i) = L(j, ind_0);
            }
        }
    }

    NumericMatrix Theta = (w_A > 0.0) ? NumericMatrix(K_max, n) : NumericMatrix();
    NumericMatrix Lambda = (w_L > 0.0) ? NumericMatrix(K_max, n) : NumericMatrix();
    if (w_A > 0.0) {
        for (k = 0; k < K; k++) {
            Theta(k, _) = rep(1.0, n);
        }
    }
    if (w_L > 0.0) {
        for (k = 0; k < K; k++) {
            Lambda(k, _) = rep(1.0, n);
        }
    }

    IntegerMatrix z_store(iter, m);
    IntegerVector K_store(iter * m);
    IntegerMatrix s_sr_store(iter, m);
    IntegerMatrix r_sr_store(iter, m);
    IntegerMatrix s_store(iter, m);
    IntegerMatrix r_store(iter, m);
    arma::cube Theta_store(K_max, n, iter);
    arma::cube Lambda_store(K_max, n, iter);

    NumericVector post_prob_sr(iter);
    NumericVector post_prob_s(iter);
    NumericVector post_prob_r(iter);
    NumericVector r_ppi(m);
    NumericVector post_prob_theta(iter);
    NumericVector post_prob_lambda(iter);

    arma::cube prob_sr_store(m, n * 2, iter);
    arma::cube prob_s_store(m, n, iter);
    NumericMatrix prob_r_store(iter, m);

    double accept_z = 0.0, accept_theta = 0.0, accept_lambda = 0.0;
    int burn_theta = 0, burn_lambda = 0;

    int burn = iter / 2, prog = 10;
                    
    for (t = 0; t < iter; t++) {
        double temp = pow(0.999, t);

        for (j = 0; j < m; j++) {
            if (exist(z, j)) {
                int K_curr = K; 
                int K_new = K_curr + 1;

                if (K_new > K_max) { 
                    goto store_K; 
                }
                
                int zj_new = K_curr;
                int zj_tmp = z(j);    

                NumericVector theta_tmp;
                if (w_A > 0.0) {
                    theta_tmp = Rcpp::rgamma(n, a_theta, 1.0/b_theta); 
                }
                NumericVector lambda_tmp;
                if (w_L > 0.0) {
                    lambda_tmp = Rcpp::rgamma(n, a_lambda, 1.0/b_lambda); 
                }
                
                double r_z0 = 0.0;
                if (w_A > 0.0) {
                    r_z0 += w_ldirich_full(A_1(j, _), theta_tmp, w_A);
                    r_z0 -= w_ldirich_full(A_1(j, _), Theta(zj_tmp, _), w_A);
                }
                if (w_L > 0.0) {
                    r_z0 += w_ldirich_full(L_1(j, _), lambda_tmp, w_L);
                    r_z0 -= w_ldirich_full(L_1(j, _), Lambda(zj_tmp, _), w_L);
                }
                
                r_z0 += v_m(K_curr) - v_m(K_curr - 1); 
                r_z0 += log(alpha) - log(m - 1.0 + K_curr * alpha);

                double r_z = std::max(r_z0, -999999.0);
                
                if (r_z >= log(R::runif(0, 1))) {
                    z(j) = zj_new;
                    K = K_new;
            
                    if (w_A > 0.0) {
                        Theta(K_curr, _) = theta_tmp;
                    }
                    if (w_L > 0.0) {
                        Lambda(K_curr, _) = lambda_tmp;
                    }
                    
                    if (t > burn) {
                        accept_z += 1;
                    }
                }
            }
            else {
                if (K > K_max) { 
                    goto store_K; 
                }

                IntegerVector z_tmp = z;
                z_tmp.erase(z_tmp.begin() + j); 

                IntegerVector znew_int = unique_z(z_tmp);
                int Knew_int = K - 1; 

                NumericVector prob_z(Knew_int);
                for (k = 0; k < Knew_int; k++) {
                    int k_idx = znew_int(k);
                    double prob_z0 = log(num_k(z_tmp, k_idx) + alpha);
                    prob_z(k) = std::max(prob_z0, -999999.0);
                }
                
                prob_z = (prob_z - max(prob_z)) / temp; 
                prob_z = exp(prob_z) / sum(exp(prob_z));
                
                int zj_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(znew_int, 1, TRUE, prob_z));
                int zj_tmp = z(j); 

                double r_z0 = 0.0;
                if (w_A > 0.0) {
                    r_z0 += w_ldirich_full(A_1(j, _), Theta(zj_new, _), w_A);
                    r_z0 -= w_ldirich_full(A_1(j, _), Theta(zj_tmp, _), w_A);
                }
                if (w_L > 0.0) {
                    r_z0 += w_ldirich_full(L_1(j, _), Lambda(zj_new, _), w_L);
                    r_z0 -= w_ldirich_full(L_1(j, _), Lambda(zj_tmp, _), w_L);
                }

                r_z0 += v_m(Knew_int - 1) - v_m(Knew_int); 
                r_z0 += log(m - 1.0 + Knew_int * alpha) - log(alpha);

                double r_z = std::max(r_z0, -999999.0);

                if (r_z >= log(R::runif(0, 1))) {
                    z(j) = zj_new; 
                    for (j_idx = 0; j_idx < m; j_idx++) {
                        if (z(j_idx) > zj_tmp) {
                            z(j_idx) -= 1;
                        }
                    }

                    K = Knew_int;
                
                    for (k_idx = zj_tmp; k_idx < Knew_int; k_idx++) {
                        if (w_A > 0.0) {
                            Theta(k_idx, _) = Theta(k_idx + 1, _);
                        }
                        if (w_L > 0.0) {
                            Lambda(k_idx, _) = Lambda(k_idx + 1, _);
                        }
                    }
                    if (w_A > 0.0) { Theta(Knew_int, _) = rep(0.0, n); }
                    if (w_L > 0.0) { Lambda(Knew_int, _) = rep(0.0, n); }
                    
                    if (t > burn) {
                        accept_z += 1;
                    }
                }
            }

            store_K:
                ind_K += 1;
                K_store(ind_K) = K;
        }

        IntegerVector m_j_k(K);
        for (j = 0; j < m; j++) {
            m_j_k[z[j]]++;
        }
        IntegerVector z_int = seq(0, K - 1);
        for (j = 0; j < m; j++) {
            if (exist(z, j)) {
                m_j_k[z[j]]--; 
                
                NumericVector prob_k(K);
                for (k = 0; k < K; k++) {
                    double prob_k0 = log(m_j_k[k] + alpha);

                    if (w_A > 0.0) {
                        prob_k0 += w_ldirich_full(A_1(j, _), Theta(k, _), w_A);
                    }
                    if (w_L > 0.0) {
                        prob_k0 += w_ldirich_full(L_1(j, _), Lambda(k, _), w_L);
                    }
                    
                    prob_k(k) = std::max(prob_k0, -999999.0);
                }
            
                prob_k = (prob_k - max(prob_k)) / temp; 
                prob_k = exp(prob_k) / sum(exp(prob_k));
                
                int z_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(z_int, 1, TRUE, prob_k));
                
                z(j) = z_new;
                m_j_k[z_new]++; 
            }
        }
        ind_z += 1;
        z_store(ind_z, _) = z;

        if (est_sr) {
            NumericVector A1_new(n);
            NumericVector L1_new(n);

            for (j = 1; j < m; j++) {
                NumericVector prob_sr(n * 2);

                for(r_0 = 0; r_0 < 2; r_0++) {
                    for (s_0 = 0; s_0 < n; s_0++) {
                        for (i = 0; i < n; i++) {
                            int ind_sr = ((s_0 + i * (1 - 2 * r_0)) + n) % n;

                            if (w_A > 0.0) {
                                A1_new(i) = A(j, ind_sr);
                            }
                            if (w_L > 0.0) {
                                L1_new(i) = L(j, ind_sr);
                            }
                        }
                        
                        int ind_prob_sr = s_0 + n * r_0;

                        if (w_A > 0.0) {
                            prob_sr(ind_prob_sr) += w_ldirich(A1_new, Theta(z(j), _), w_A);
                        }
                        if (w_L > 0.0) {
                            prob_sr(ind_prob_sr) += w_ldirich(L1_new, Lambda(z(j), _), w_L);
                        }
                        
                        prob_sr(ind_prob_sr) += lprior_s(s_0, n - 1, alpha_s, beta_s);
                        prob_sr(ind_prob_sr) += lbern(r_0, omega);
                        prob_sr(ind_prob_sr) = std::max(prob_sr(ind_prob_sr), -999999.0);
                    }
                }
                
                prob_sr = (prob_sr - max(prob_sr)) / temp;
                prob_sr = exp(prob_sr) / sum(exp(prob_sr));

                prob_sr_store.slice(t).row(j) = as<arma::rowvec>(prob_sr);
                
                int sr_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(sr_int, 1, TRUE, prob_sr));
                
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
            NumericVector A1_new(n);
            NumericVector L1_new(n);

            for (j = 1; j < m; j++) {
                NumericVector prob_s(n);

                for (s_0 = 0; s_0 < n; s_0++) {
                    for (i = 0; i < n; i++) {
                        int ind_s = ((s_0 + i * (1 - 2 * r(j))) + n) % n;

                        if (w_A > 0.0) {
                            A1_new(i) = A(j, ind_s);
                        }
                        if (w_L > 0.0) {
                            L1_new(i) = L(j, ind_s);
                        }
                    }

                    if (w_A > 0.0) {
                        prob_s(s_0) += w_ldirich(A1_new, Theta(z(j), _), w_A);
                    }
                    if (w_L > 0.0) {
                        prob_s(s_0) += w_ldirich(L1_new, Lambda(z(j), _), w_L);
                    }

                    prob_s(s_0) += lprior_s(s_0, n - 1, alpha_s, beta_s);
                    prob_s(s_0) = std::max(prob_s(s_0), -999999.0);
                }
                
                prob_s = (prob_s - max(prob_s)) / temp;
                prob_s = exp(prob_s) / sum(exp(prob_s));

                prob_s_store.slice(t).row(j) = as<arma::rowvec>(prob_s);
                
                int s_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(s_int, 1, TRUE, prob_s));
                
                s(j) = s_new;
            }
        }

        if (est_r) {
            NumericVector A1_new(n);
            NumericVector L1_new(n);

            for (j = 1; j < m; j++) {
                NumericVector prob_r(2);

                for (r_0 = 0; r_0 < 2; r_0++) {
                    for (i = 0; i < n; i++) {
                        int ind_r = ((s(j) + i * (1 - 2 * r_0)) + n) % n;

                        if (w_A > 0.0) {
                            A1_new(i) = A(j, ind_r);
                        }
                        if (w_L > 0.0) {
                            L1_new(i) = L(j, ind_r);
                        }
                    }

                    if (w_A > 0.0) {
                        prob_r(r_0) += w_ldirich(A1_new, Theta(z(j), _), w_A);
                    }
                    if (w_L > 0.0) {
                        prob_r(r_0) += w_ldirich(L1_new, Lambda(z(j), _), w_L);
                    }

                    prob_r(r_0) += lbern(r_0, omega);
                    prob_r(r_0) = std::max(prob_r(r_0), -999999.0);
                }
                
                prob_r = (prob_r - max(prob_r)) / temp;
                prob_r = exp(prob_r) / sum(exp(prob_r));

                prob_r_store(t, j) = prob_r(1);
                
                int r_new = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(r_int, 1, TRUE, prob_r));
                
                r(j) = r_new;
            }
        }

        for (j = 1; j < m; j++) {
            for (i = 0; i < n; i++) {
                int ind_new = ((s(j) + i * (1 - 2 * r(j))) + n) % n;
                
                if (w_A > 0.0) {
                    A_1(j, i) = A(j, ind_new);
                }
                if (w_L > 0.0) {
                    L_1(j, i) = L(j, ind_new);
                }
            }
        }

        if (w_A > 0.0) {
            for (k = 0; k < K; k++) {
                NumericVector Theta_new = Theta(k, _);
                
                for (i = 0; i < n; i++) {
                    double theta_ki_tmp = Theta(k, i);
                
                    double theta_ki_new = exp(R::rnorm(log(theta_ki_tmp), tau_theta));
                    Theta_new(i) = theta_ki_new;
                
                    double r_theta = 0.0;
                    for (j = 0; j < m; j++) {
                        if (z(j) == k) {
                            r_theta += w_ldirich_full(A_1(j, _), Theta_new, w_A);
                            r_theta -= w_ldirich_full(A_1(j, _), Theta(k, _), w_A);
                        }
                    }
                
                    r_theta += logamma(theta_ki_new, a_theta, b_theta);
                    r_theta -= logamma(theta_ki_tmp, a_theta, b_theta);

                    r_theta += log(theta_ki_new) - log(theta_ki_tmp);
                
                    if (r_theta >= log(R::runif(0, 1))) {
                        Theta(k, i) = theta_ki_new;
                    
                        if (t > burn) {
                            accept_theta += 1;
                        }
                    }
                    else {
                        Theta_new(i) = theta_ki_tmp;
                    }
                }
            }

            if (t > burn) {
                burn_theta += K * n;
            }
        }

        if (w_L > 0.0) {
            for (k = 0; k < K; k++) {
                NumericVector Lambda_new = Lambda(k, _);
                
                for (i = 0; i < n; i++) {
                    double lambda_ki_tmp = Lambda(k, i);
                
                    double lambda_ki_new = exp(R::rnorm(log(lambda_ki_tmp), tau_lambda));
                    Lambda_new(i) = lambda_ki_new;
                
                    double r_lambda = 0.0;
                    for (j = 0; j < m; j++) {
                        if (z(j) == k) {
                            r_lambda += w_ldirich_full(L_1(j, _), Lambda_new, w_L);
                            r_lambda -= w_ldirich_full(L_1(j, _), Lambda(k, _), w_L);
                        }
                    }
                
                    r_lambda += logamma(lambda_ki_new, a_lambda, b_lambda);
                    r_lambda -= logamma(lambda_ki_tmp, a_lambda, b_lambda);

                    r_lambda += log(lambda_ki_new) - log(lambda_ki_tmp);
                
                    if (r_lambda >= log(R::runif(0, 1))) {
                        Lambda(k, i) = lambda_ki_new;
                    
                        if (t > burn) {
                            accept_lambda += 1;
                        }
                    }
                    else {
                        Lambda_new(i) = lambda_ki_tmp;
                    }
                }
            }

            if (t > burn) {
                burn_lambda += K * n;
            }
        }

        double post_prob_sr_0 = 0.0;
        if (est_sr || est_s || est_r) {
            for (j = 1; j < m; j++) {
                if (w_A > 0.0) {
                    post_prob_sr_0 += w_ldirich(A_1(j, _), Theta(z(j), _), w_A);
                }
                if (w_L > 0.0) {
                    post_prob_sr_0 += w_ldirich(L_1(j, _), Lambda(z(j), _), w_L);  
                }
            }    
        }

        double post_prob_s_0 = 0.0;
        if (est_sr || est_s) {
            for (j = 1; j < m; j++) {     
                post_prob_s_0 += lprior_s(s(j), n - 1, alpha_s, beta_s);
            }
        }

        double post_prob_r_0 = 0.0;
        if (est_sr || est_r) {
            for (j = 1; j < m; j++) {      
                post_prob_r_0 += lbern(r(j), omega);
            }
        }

        if (est_r) {
            if (t > burn) {
                for (j = 0; j < m; j++) {
                    r_ppi(j) += r(j); 
                }
            } 
        }

        double post_prob_theta_0 = 0.0;
        if (w_A > 0.0) {
            for (j = 0; j < m; j++) {
                post_prob_theta_0 += w_ldirich_full(A_1(j, _), Theta(z(j), _), w_A);
                
                for (i = 0; i < n; i++) {
                    post_prob_theta_0 += logamma(Theta(z(j), i), a_theta, b_theta);
                }
            }
        }

        double post_prob_lambda_0 = 0.0;
        if (w_L > 0.0) {
            for (j = 0; j < m; j++) {
                post_prob_lambda_0 += w_ldirich_full(L_1(j, _), Lambda(z(j), _), w_L);
                
                for (i = 0; i < n; i++) {
                    post_prob_lambda_0 += logamma(Lambda(z(j), i), a_lambda, b_lambda);
                }
            }
        }

        if (est_sr) {
            s_sr_store(t, _) = s;
            r_sr_store(t, _) = r;
        }
        if (est_s) { s_store(t, _) = s; }
        if (est_r) { r_store(t, _) = r; }
        if (w_A > 0.0) { Theta_store.slice(t) = as<arma::mat>(Theta); }
        if (w_L > 0.0) { Lambda_store.slice(t) = as<arma::mat>(Lambda); }

        if (est_sr || est_s || est_r) { post_prob_sr(t) = post_prob_sr_0; }
        if (est_s) { post_prob_s(t) = post_prob_s_0; }
        if (est_r) { post_prob_r(t) = post_prob_r_0; }
        if (w_A > 0.0) { post_prob_theta(t) = post_prob_theta_0; }
        if (w_L > 0.0) { post_prob_lambda(t) = post_prob_lambda_0; }

        if(((t * 100) / (iter - 1)) == prog) {
            Rcout << prog << "%" << std::endl;
            prog += 10;
        }
    }

    if (est_r) { r_ppi = r_ppi / burn; }

    accept_z = accept_z / (burn * m);
    if (w_A > 0.0) { accept_theta = accept_theta / burn_theta; }
    if (w_L > 0.0) { accept_lambda = accept_lambda / burn_lambda; }

    return Rcpp::List::create(Rcpp::Named("z_store") = z_store, Rcpp::Named("K_store") = K_store,
                              Rcpp::Named("s_sr_store") = s_sr_store, Rcpp::Named("r_sr_store") = r_sr_store, Rcpp::Named("s_store") = s_store, Rcpp::Named("r_store") = r_store,
                              Rcpp::Named("Theta_store") = Theta_store, Rcpp::Named("Lambda_store") = Lambda_store,
                              Rcpp::Named("post_prob_sr") = post_prob_sr, Rcpp::Named("post_prob_s") = post_prob_s, Rcpp::Named("post_prob_r") = post_prob_r, Rcpp::Named("r_ppi") = r_ppi, Rcpp::Named("post_prob_theta") = post_prob_theta, Rcpp::Named("post_prob_lambda") = post_prob_lambda,
                              Rcpp::Named("prob_sr_store") = prob_sr_store, Rcpp::Named("prob_s_store") = prob_s_store, Rcpp::Named("prob_r_store") = prob_r_store,
                              Rcpp::Named("accept_z") = accept_z, Rcpp::Named("accept_theta") = accept_theta, Rcpp::Named("accept_lambda") = accept_lambda);
}