library(mclust)  

# -------------------------------------------------------------------------

#z

# Compute MAP and PPM estimates after burn-in
map_estimate <- calc_MAP(iter_mat)  # MAP estimate
PPM <- calc_PPM(iter_mat)  # Pairwise Probability Matrix
ppm_estimate <- calc_PPM_estimate(iter_mat, PPM)  # Best clustering based on PPM

# Compute Adjusted Rand Index (ARI) for both estimates
ari_map <- adjustedRandIndex(map_estimate, est_0)
ari_ppm <- adjustedRandIndex(ppm_estimate, est_0)

# -------------------------------------------------------------------------

#s or r

# Compute MAP and PPM estimates after burn-in
map_estimate <- calc_MAP(iter_mat)  # MAP estimate

# Compute Adjusted Rand Index (ARI) for both estimates
ari_map <- adjustedRandIndex(map_estimate, est_0)

# -------------------------------------------------------------------------

#Theta and Lambda

theta <- res$Theta_store
lambda <- res$Lambda_store
z_store <- res$z_store

# Compute burn-in threshold
burnin_index <- floor(length(theta) / 2) + 1

# Apply burn-in
post_burnin_theta <- theta[burnin_index:length(theta)]
post_burnin_lambda <- lambda[burnin_index:length(lambda)]
post_burnin_z_store <- z_store[burnin_index:nrow(z_store), , drop = FALSE]  
# Identify valid post-burn-in iterations where K matches K_final
K_final <- length(unique(map_estimate)) #or PPM
valid_indices <- which(sapply(1:length(post_burnin_theta), function(u) nrow(post_burnin_theta[[u]]) == K_final))

# Filter theta, lambda, and z_store using valid post-burn-in indices
filtered_theta <- post_burnin_theta[valid_indices]
filtered_lambda <- post_burnin_lambda[valid_indices]
filtered_z_store <- post_burnin_z_store[valid_indices, , drop = FALSE]

# Apply label-switching correction
corrected <- switch_label(filtered_z_store, K_final, filtered_theta, filtered_lambda)
filtered_z_store <- corrected$z_store
filtered_theta <- corrected$Theta_store
filtered_lambda <- corrected$Lambda_store

# Compute MAP Estimate After Label Switching 
map_estimate <- apply(filtered_z_store, 2, function(x) { 
  as.numeric(names(which.max(table(x))))  
})

# Cpmpute Posterior Mean for Theta and Lambda
theta_mean <- apply(simplify2array(filtered_theta), c(1, 2), mean)
lambda_mean <- apply(simplify2array(filtered_lambda), c(1, 2), mean)

theta_mean_shape <- mean_k(theta_mean)
theta_var_cov <- var_cov_k(theta_mean)
theta_var <- t(apply(theta_var_cov, 3, diag))

lambda_mean_shape <- mean_k(lambda_mean)
lambda_var_cov <- var_cov_k(lambda_mean)
lambda_var <- t(apply(lambda_var_cov, 3, diag))
