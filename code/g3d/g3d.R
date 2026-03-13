library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(RcppEigen)
library(compositions)
library(mclust)
library(R.matlab)

sourceCpp("./code/mfm/bacon_mfm.cpp")

load("./data/g3d/g3d.RData")

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

res <- bacon_mfm(A, L, w_A = 1, w_L = 1, 
                 z_0, K_max = K_max, alpha = 10, 
                 est_sr = TRUE, est_s = TRUE, est_r = TRUE, alpha_s = 1, beta_s = 1, omega = 0.5, 
                 a_theta = 0.001, b_theta = 0.001, a_lambda = 0.001, b_lambda = 0.001, tau_theta = 0.1, tau_lambda = 0.1, 
                 iter = iter)

#save(res, file = "./res/g3d/res.RData")

load("./res/g3d/res.RData")

z_store <- res$z_store
max_K <- max(res$K_store)
prob_z <- matrix(0, m, max_K)
for(j in 1:m) {
  for(k in 1:max_K) {
    prob_z[j, k] <- mean((z_store[(burn + 1):iter, j] + 1) == k)
  }
}
z_hat <- as.vector(apply(prob_z, 1, function(x) which(x == max(x))[1]))
K_hat <- length(unique(z_hat))

z_2 <- ifelse(z %in% c(0, 1), 0,
              ifelse(z %in% c(2, 3), 1, NA))

adjustedRandIndex(z_hat, z_2)

Theta_store <- res$Theta_store
theta_burn <- Theta_store[, , (burn + 1):iter]
valid_iters <- apply(theta_burn, 3, function(iteration_matrix) {
  num_active <- sum(rowSums(iteration_matrix) > 0)
  return(num_active == K_hat)
})
theta_filter <- theta_burn[1:K_hat, , valid_iters, drop = FALSE]
theta_hat <- apply(theta_filter, c(1, 2), mean)
#Remove low range of motion cluster
theta_hat <- theta_hat[-1, ]

output_path <- "./res/g3d/actionValuesnew.mat"

kickLeftVals  <- c(theta_hat[1, 1], theta_hat[1, 5], theta_hat[1, 4], theta_hat[1, 3], theta_hat[1, 2])
kickRightVals <- theta_hat[1, ]
punchLeftVals <- c(theta_hat[2, 1], theta_hat[2, 5], theta_hat[2, 4], theta_hat[2, 3], theta_hat[2, 2])
punchRightVals <- theta_hat[2, ]

writeMat(output_path, 
         kickLeftVals = kickLeftVals, 
         kickRightVals = kickRightVals, 
         punchLeftVals = punchLeftVals, 
         punchRightVals = punchRightVals)
