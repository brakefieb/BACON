library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(RcppEigen)
library(compositions)
library(R.matlab)

sourceCpp("./code/mfm/bacon_mfm.cpp")

load("./data/tri_mesh/tri_mesh.RData")

K_0 <- 10
A_clr <- clr(A)
L_clr <- clr(L)
X_combined <- cbind(A_clr, L_clr)
init_kmeans <- kmeans(X_combined, centers = K_0, nstart = 25)
init_labels <- init_kmeans$cluster
z_0 <- init_labels - 1

K_max <- 30
iter <- 1000
burn <- iter/2

alpha <- c(10, 30, 50)
for (a in alpha) {
  res <- bacon_mfm(A, L, w_A = 1, w_L = 0, 
                   z_0, K_max = K_max, alpha = a, 
                   est_sr = TRUE, est_s = TRUE, est_r = TRUE, alpha_s = 1, beta_s = 1, omega = 0.5, 
                   a_theta = 0.001, b_theta = 0.001, a_lambda = 0.001, b_lambda = 0.001, tau_theta = 0.1, tau_lambda = 0.1, 
                   iter = iter)
  
  z_store <- res$z_store
  max_K <- max(res$K_store)
  prob_z <- matrix(0, m, max_K)
  for(j in 1:m) {
    for(k in 1:max_K) {
      prob_z[j, k] <- mean((z_store[(burn + 1):iter, j] + 1) == k)
    }
  }
  z_map <- as.vector(apply(prob_z, 1, function(x) which(x == max(x))[1]))
  
  writeMat(paste0("./res/tri_mesh/z_est_", a, ".mat"), groups = z_map)
  
  rm(res)
  
  print(a)
}
