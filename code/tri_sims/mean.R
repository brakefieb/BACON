library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(RcppEigen)
library(compositions)
library(mclust)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# -------------------------------------------------------------------------

sourceCpp("./code/mfm/bacon_mfm.cpp")

# -------------------------------------------------------------------------

load("./data/tri_sims/tri_sims_sr.RData")

ind <- 5
A <- A2_tensor[, , ind]
L <- L2_tensor[, , ind]

# -------------------------------------------------------------------------

K_0 <- 10
# 1. Transform data
# A and L are your m x n matrices
A_clr <- clr(A)
L_clr <- clr(L)
# 2. Combine (scaling if necessary)
X_combined <- cbind(A_clr, L_clr)
# 3. Cluster for Initial Values
init_kmeans <- kmeans(X_combined, centers = K_0, nstart = 25)
# 4. Extract Initial Parameters for Bayesian Algo
init_labels <- init_kmeans$cluster
z_0 <- init_labels - 1

K_max <- 30
iter <- 5000
burn <- iter/2

res <- bacon_mfm(A, L, w_A = 1, w_L = 0, 
                 z_0, K_max = K_max, alpha = 30, 
                 est_sr = TRUE, est_s = TRUE, est_r = TRUE, alpha_s = 1, beta_s = 1, omega = 0.5, 
                 a_theta = 0.001, b_theta = 0.001, a_lambda = 0.001, b_lambda = 0.001, tau_theta = 0.1, tau_lambda = 0.1, 
                 iter = iter)

# -------------------------------------------------------------------------

load("./res/tri_sims/res.RData")

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

Theta_store <- res$Theta_store
theta_burn <- Theta_store[, , (burn + 1):iter]
valid_iters <- apply(theta_burn, 3, function(iteration_matrix) {
  num_active <- sum(rowSums(iteration_matrix) > 0)
  return(num_active == K_hat)
})
theta_filter <- theta_burn[1:K_hat, , valid_iters, drop = FALSE]
theta_hat <- apply(theta_filter, c(1, 2), mean)

# -------------------------------------------------------------------------

source("./code/tri_sims/funcs.R")

# -------------------------------------------------------------------------

load("./res/tri_sims/theta_hat.RData")

a_lower <- rep(0, n)
b_upper <- rep(1, n)

theta_mean <- matrix(0, K_hat, n)
theta_var <- matrix(0, K_hat, n)
lambda <- matrix(0, K_hat, n)
pc_list <- list()

T <- 5000
a_1 <- rtdirichlet(n = T, eta = theta_hat[1, ], a = a_lower, b = b_upper)
a_2 <- rtdirichlet(n = T, eta = theta_hat[2, ], a = a_lower, b = b_upper)
a_3 <- rtdirichlet(n = T, eta = theta_hat[3, ], a = a_lower, b = b_upper)

theta_mean[1, ] <- colMeans(a_1)
theta_mean[2, ] <- colMeans(a_2)
theta_mean[3, ] <- colMeans(a_3)

theta_rad <- theta_mean*(n - 2)*pi

theta_var[1, ] <- apply(a_1, 2, var)
theta_var[2, ] <- apply(a_2, 2, var)
theta_var[3, ] <- apply(a_3, 2, var)
  
lambda[1, ] <- sin(theta_rad[1, ]) / sum(sin(theta_rad[1, ]))
lambda[1, ] <- c(lambda[1, 3], lambda[1, -3])
lambda[1, ] <- lambda[1, ] / sum(lambda[1, ]) 

lambda[2, ] <- sin(theta_rad[2, ]) / sum(sin(theta_rad[2, ]))
lambda[2, ] <- c(lambda[2, 3], lambda[2, -3])
lambda[2, ] <- lambda[2, ] / sum(lambda[2, ])

lambda[3, ] <- sin(theta_rad[3, ]) / sum(sin(theta_rad[3, ]))
lambda[3, ] <- c(lambda[3, 3], lambda[3, -3])
lambda[3, ] <- lambda[3, ] / sum(lambda[3, ])

pc_list[[1]] <- t(v(theta_rad[1, ], lambda[1, ], 0))
pc_list[[2]] <- t(v(theta_rad[2, ], lambda[2, ], 0))
pc_list[[3]] <- t(v(theta_rad[3, ], lambda[3, ], 0))

# -------------------------------------------------------------------------

K <- length(pc_list)

theta_var_list <- list(
  theta_var[1, ],      
  theta_var[2, ],    
  theta_var[3, ]     
)

# Combined data containers
combined_polygon_data <- data.frame()
combined_vertex_data <- data.frame() # Specifically for the points
# Loop over groups and build data
for (k in c(2, 1, 3)) {
  pc_matrix <- pc_list[[k]]
  n <- nrow(pc_matrix) - 1  # number of unique vertices
  
  # Polygon outline (The skeleton)
  poly_data <- data.frame(
    x = pc_matrix[, 1],
    y = pc_matrix[, 2],
    Group = factor(k)
  )
  
  # Vertex data (The angles)
  vertex_data <- data.frame(
    x = pc_matrix[1:n, 1],
    y = pc_matrix[1:n, 2],
    Variance = theta_var_list[[k]],
    Group = factor(k)
  )
  
  combined_polygon_data <- rbind(combined_polygon_data, poly_data)
  combined_vertex_data <- rbind(combined_vertex_data, vertex_data)
}

# Set shared color scale limits
common_limits <- range(unlist(theta_var_list))
my_custom_labels <- c(
  "2" = "1",
  "1" = "2",
  "3" = "3"
)

p <- ggplot() +
  # Triangle sides in black
  geom_path(data = combined_polygon_data, aes(x = x, y = y, group = Group), 
            color = "black", linewidth = 1) + # (Kept this as linewidth to prevent ggplot warnings)
  
  # Vertices (angles) colored by Variance
  geom_point(data = combined_vertex_data, aes(x = x, y = y, color = Variance), 
             size = 4) +
  
  # Your original settings
  scale_color_gradient(low = "blue", high = "red", limits = common_limits, name = "Variance") +
  coord_fixed(ratio = 1) +
  
  # ONLY the text mapping is changed here
  facet_grid(. ~ Group, labeller = labeller(Group = my_custom_labels)) +
  
  scale_x_continuous(name = "x") +
  scale_y_continuous(name = "y") +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line()
  )

p
