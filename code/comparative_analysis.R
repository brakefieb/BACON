load("transformed_shapes.RData")

m_z <- rep(50, 3)
z <- c(rep(0, m_z[1]), rep(1, m_z[2]), rep(2, m_z[3]))

# -------------------------------------------------------------------------

#https://github.com/jd-strait/LCESA

library(R.matlab)

mat_array <- simplify2array(transformed_shapes)

writeMat("E:/LCESA-master/transformed_shapes.mat", transformed_shapes = mat_array)

# -------------------------------------------------------------------------

#https://github.com/estfernan/SAFARI

library(SAFARI)
library(sp)

resolution <- 1000
feats_list <- list()
chains_list <- list()
for (i in 1:length(transformed_shapes)) {
  vertices <- transformed_shapes[[i]]
  
  # Determine the bounding box of the polygon
  x_range <- range(vertices[,1])
  y_range <- range(vertices[,2])
  
  # Define resolution: adjust resolution to control the grid density
  x_seq <- seq(from = x_range[1], to = x_range[2], length.out = resolution)
  y_seq <- seq(from = y_range[1], to = y_range[2], length.out = resolution)
  
  # Create a grid of (x, y) coordinates
  grid <- expand.grid(x = x_seq, y = y_seq)
  
  # Check each grid point: returns 0 (outside), 1 (boundary), or 2 (inside)
  inside <- point.in.polygon(grid$x, grid$y, vertices[,1], vertices[,2])
  
  # Convert to binary: set any point with a value > 0 as 1 (inside) and 0 (outside)
  binary_matrix <- matrix(as.integer(inside > 0), nrow = resolution, ncol = resolution)
  
  img_seg <- binary.segmentation(binary_matrix,
                                 id = 10,
                                 filter = 10,
                                 k = 3,
                                 categories = c("geometric", "boundary", "topological"))
  
  feats_list[[i]] <- img_seg$desc
  
  chains_list[[i]] <- img_seg$plg.chains
  
  rm(img_seg)
  
  print(i)
}

feats_trimmed <- lapply(feats_list, function(x) x[, -c(1,2)])
feats_matrix <- do.call(rbind, feats_trimmed)
feats <- feats_matrix[, -c(28)]
features_scaled <- feats %>% mutate_all(scale)
save(features_scaled, file = "features_scaled.RData")

chains_extracted <- lapply(chains_list, `[[`, "1")
save(chains_extracted, file = "chains_extracted.RData")

# -------------------------------------------------------------------------

#https://github.com/kevinwjin/SAFARI-cluster-analysis

library(dplyr)
library(compositions)
library(ggfortify)
library(ggplot2)
library(mclust)
library(factoextra)

#PCA

pca_scaled <- prcomp(features_scaled)

# Create a data frame with PCA results and cluster information
pca_df <- data.frame(features_scaled, Cluster = factor(z + 1))

p_0 <- autoplot(pca_scaled, data = pca_df, colour = "Cluster", 
                frame = TRUE, frame.type = "norm", frame.colour = "Cluster") +
  scale_color_manual(values = c("#2AA02B", "#FF7F0F", "#2077B4")) +
  scale_fill_manual(values = c("#2AA02B", "#FF7F0F", "#2077B4")) + theme_minimal()

ggsave(p_0, filename = "tri_sims_pca.png", bg = "white")

#k-Means

max_k <- 10 
z_0 <- 3 #true # of clusters
n <- nrow(features_scaled)

kmeans_mat <- matrix(nrow = max_k, ncol = n, byrow = TRUE)
# Generate k-means clustering from scaled features
for (k in 1:max_k) {
  kmeans_mat[k, ] <- kmeans(
    x = features_scaled,
    centers = k
    # nstart = 20,
    # iter.max = 30
  )[["cluster"]]
}
kmeans_scaled <- as.data.frame(kmeans_mat)

# Calculate ARI values for k-means over k = 1:k_max
kmeans_ari <- matrix(nrow = max_k, ncol = 1, byrow = TRUE)
for (i in 1:max_k) {
  # ARI for k-means clustering of scaled features
  kmeans_ari[i, 1] <- adjustedRandIndex(
    unlist(kmeans_scaled[i, ]),
    z
  )
}
kmeans_ari <- as.data.frame(kmeans_ari)
names(kmeans_ari) <- c("ARI_scaled")
k_values <- 1:max_k
kmeans_ari <- mutate(kmeans_ari, k_values = as.numeric(row.names(kmeans_ari)))

#Hierarchical

# Generate hierarchical clusters
hier_scaled_tree <- hclust(dist(features_scaled, method = "euclidean"),
                           method = "complete" # Complete linkage
)

hier_mat <- matrix(nrow = max_k, ncol = n, byrow = TRUE)
# Generate hierarchical clustering from scaled features
for (k in 1:max_k) {
  hier_mat[k, ] <- cutree(hier_scaled_tree, k = k)
}
hier_scaled <- as.data.frame(hier_mat)

# Calculate ARI values for hierarchical over k = 1:max_k 
hier_ari <- matrix(nrow = max_k, ncol = 1, byrow = TRUE)
for (i in 1:max_k) {
  # ARI for hierarchical clustering of scaled features
  hier_ari[i, 1] <- adjustedRandIndex(
    unlist(hier_scaled[i, ]),
    z
  )
}
hier_ari <- as.data.frame(hier_ari)
names(hier_ari) <- c("ARI_scaled")
k_values <- 1:max_k
hier_ari <- mutate(hier_ari, k_values = as.numeric(row.names(hier_ari)))

#GMM

gmm_mat <- matrix(nrow = max_k, ncol = n, byrow = TRUE)
# Generate GMM clustering from scaled features
for (k in 1:max_k) {
  gmm_mat[k, ] <- Mclust(features_scaled, G = k)$classification
}
gmm_scaled <- as.data.frame(gmm_mat)

# Calculate ARI values for GMM 
gmm_ari <- matrix(nrow = max_k, ncol = 1, byrow = TRUE)
for (i in 1:max_k) {
  # ARI for GMM clustering of scaled features
  gmm_ari[i, 1] <- adjustedRandIndex(
    unlist(gmm_scaled[i, ]),
    z
  )
}
gmm_ari <- as.data.frame(gmm_ari)
names(gmm_ari) <- c("ARI_scaled")
k_values <- 1:max_k
gmm_ari <- mutate(gmm_ari, k_values = as.numeric(row.names(gmm_ari)))

# Visualize accuracy of all clustering methods
accuracy <- right_join(kmeans_ari, gmm_ari, by = "k_values")
accuracy <- right_join(accuracy, hier_ari, by = "k_values")
names(accuracy) <- c(
  "kmeans_ari_scaled",
  "k_values",
  "gmm_ari_scaled",
  "hier_ari_scaled"
)

plot <- ggplot(accuracy, aes(x = k_values)) +
  geom_vline(xintercept = z_0, color = "black", linetype = "solid") +
  geom_smooth(aes(y = kmeans_ari_scaled, color = "k-Means"), se = FALSE) +
  geom_point(aes(y = kmeans_ari_scaled, color = "k-Means")) +
  geom_smooth(aes(y = hier_ari_scaled, color = "Hierarchical"), se = FALSE) +
  geom_point(aes(y = hier_ari_scaled, color = "Hierarchical")) +
  geom_smooth(aes(y = gmm_ari_scaled, color = "GMM"), se = FALSE) +
  geom_point(aes(y = gmm_ari_scaled, color = "GMM")) +
  labs(x = "Cluster", y = "ARI", color = "Method") +
  scale_x_continuous(breaks = 1:max_k) +
  scale_y_continuous(limits = c(min(accuracy[, -c(2)]), 1)) +  # Sets y-axis from 0 to 1
  scale_color_manual(
    values = c("k-Means" = "#2AA02B", "GMM" = "#FF7F0F", "Hierarchical" = "#2077B4"),
    breaks = c("k-Means", "GMM", "Hierarchical")
  ) +
  theme_minimal()

ggsave(plot, filename = "tri_sims_clus_meth.png", bg = "white")
