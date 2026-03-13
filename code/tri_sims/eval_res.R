library(compositions)
library(dplyr)
library(factoextra)
library(ggfortify)
library(ggplot2)
library(mclust)
library(R.matlab)
library(SAFARI)
library(sp)

# -------------------------------------------------------------------------

source("./code/tri_sims/funcs.R")

# -------------------------------------------------------------------------

load("./data/tri_sims/tri_sims_sr.RData")

ind <- 5
A <- A2_tensor[, , ind]
L <- L2_tensor[, , ind]

s <- s_matrix[, ind]
r <- r_matrix[, ind]

# -------------------------------------------------------------------------

# Process and store all m observations
transformed_shapes <- list()
for (i in 1:m) {
  pc_closed <- v(A[i, ] * (n - 2) * pi, L[i, ], 0)
  pc_closed[, n + 1] <- pc_closed[, 1]
  
  rownames(pc_closed) <- c("x", "y")
  pc_closed <- t(pc_closed)
  
  # Remove the duplicate closing vertex to work with unique vertices
  unique_points <- pc_closed[-nrow(pc_closed), ]
  
  # --- Parameters ---
  s_1 <- s[i]  # New starting vertex index (0, 1, or 2); here s != 0 so vertex 2 becomes our new start
  r_1 <- r[i]  # Orientation flag: 0 keeps the order, 1 reverses the order of the subsequent vertices
  
  # --- Reorder Vertices ---
  # Get the new order (0-indexed) and convert to R's 1-indexing
  vertex_order <- get_vertex_labels(s_1, r_1)
  order_indices <- vertex_order + 1
  new_order <- unique_points[order_indices, ]
  
  # --- Translate ---
  # Shift the polygon so that the new starting vertex moves to (0,0)
  translation_vector <- new_order[1, ]
  translated_points <- sweep(new_order, 2, translation_vector, FUN = "-")
  
  # --- Rotate ---
  # Compute the vector for the first edge (from (0,0) to the second vertex)
  first_edge <- translated_points[2, ]
  # Calculate its angle (in radians)
  angle <- atan2(first_edge[2], first_edge[1])
  # Create a rotation matrix to rotate by -angle so the edge aligns with the positive x-axis
  rotation_matrix <- matrix(c(cos(-angle), -sin(-angle),
                              sin(-angle),  cos(-angle)),
                            ncol = 2, byrow = TRUE)
  # Apply the rotation to all vertices
  rotated_points <- t(rotation_matrix %*% t(translated_points))
  
  # --- Re-close the Polygon ---
  # Append the first vertex to the end to maintain the closed polygon structure
  final_polygon <- rbind(rotated_points, rotated_points[1, ])
  
  colnames(final_polygon) <- c("x", "y")
  
  # Store the first observation as reference
  if (i == 1) {
    transformed_shapes[[i]] <- final_polygon  # Keep reference shape
  } else {
    transformed_shapes[[i]] <- random_transform(final_polygon)  # Apply transformations
  }
  
  print(i)
}

# -------------------------------------------------------------------------

load("./data/tri_sims/transformed_tri.RData")

resolution <- 100
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

# Remove the first two elements from each 1x31 matrix in the list
feats_trimmed <- lapply(feats_list, function(x) x[, -c(1,2)])
# Combine all the trimmed matrices into one matrix (sample size by 29)
feats_matrix <- do.call(rbind, feats_trimmed)

# -------------------------------------------------------------------------

load("./data/tri_sims/feats_mat.RData")

feats <- feats_matrix[, -c(28)]
features_scaled <- feats %>% mutate_all(scale)

# -------------------------------------------------------------------------

# PCA

# Perform PCA
pca_scaled <- prcomp(features_scaled)
# Create a data frame with PCA results and cluster information
pca_df <- data.frame(features_scaled, Cluster = factor(z + 1))

p_0 <- autoplot(pca_scaled, data = pca_df, colour = "Cluster", 
                frame = TRUE, frame.type = "norm", frame.colour = "Cluster") +
  scale_color_manual(values = c("#2AA02B", "#FF7F0F", "#2077B4")) +
  scale_fill_manual(values = c("#2AA02B", "#FF7F0F", "#2077B4")) + theme_minimal()

p_0

# -------------------------------------------------------------------------

# k-means 

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

# Hierarchical

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

# Calculate ARI values for hierarchical over k = 1:max_k (Ground truth: k = 6)
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

# GMM

gmm_mat <- matrix(nrow = max_k, ncol = n, byrow = TRUE)

# Generate GMM clustering from scaled features
for (k in 1:max_k) {
  gmm_mat[k, ] <- Mclust(features_scaled, G = k)$classification
}
gmm_scaled <- as.data.frame(gmm_mat)

# Calculate ARI values for GMM over k = 1:100 (Ground truth: k = 70)
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

# -------------------------------------------------------------------------

source("./code/fmm/functions.R")

source("./code/fmm/bacon.R")

# -------------------------------------------------------------------------

K_vals <- 2:10

ari_sims_ppm <- numeric(length(K_vals) + 1)
ari_sims_map <- numeric(length(K_vals) + 1)
ind <- 1
for (i in K_vals) {
  ind <- ind + 1
  
  res <- bacon(L, A, i, weight_L = 0.75, weight_A = 1, estimate.s = TRUE, estimate.r = TRUE, iter = 2000, burn = 1000)
  z_ppm <- res$cluster
  z_map <- res$z_map
  
  ari_sims_ppm[ind] <- adjustedRandIndex(z_ppm, z)
  ari_sims_map[ind] <- adjustedRandIndex(z_map, z)
  
  rm(res)
  
  print(i)
}

res_sims <- cbind(1:10, ari_sims_ppm, ari_sims_map)

# -------------------------------------------------------------------------

load("./res/tri_sims/acc.RData")

load("./res/tri_sims/fmm.RData")

accuracy_new <- accuracy
accuracy_new$fmm <- res_sims[, 3]

z_0 <- 3
max_k <- 10

# -------------------------------------------------------------------------

mat_data <- readMat("./res/tri_sims/lcesa.mat")

z_lcesa <- mat_data$true.labels
pred <- mat_data$predicted.labels

# -------------------------------------------------------------------------

load("./res/tri_sims/wa.RData")

res_vec <- res_sims[res_sims$w == 500 & res_sims$w_l == 0, 5]

# -------------------------------------------------------------------------

# 1. Define values
lcesa_val <- adjustedRandIndex(pred, z_lcesa)
mfm_val <- mean(res_vec)

# 2. Plot
plot <- ggplot(accuracy_new, aes(x = k_values)) +
  geom_vline(xintercept = z_0, color = "black", linetype = "solid") +
  
  # --- CURVES ---
  geom_smooth(aes(y = fmm, color = "BACON FMM"), se = FALSE) +
  geom_point(aes(y = fmm, color = "BACON FMM")) +
  
  geom_smooth(aes(y = kmeans_ari_scaled, color = "k-Means"), se = FALSE) +
  geom_point(aes(y = kmeans_ari_scaled, color = "k-Means")) +
  
  geom_smooth(aes(y = gmm_ari_scaled, color = "GMM"), se = FALSE) +
  geom_point(aes(y = gmm_ari_scaled, color = "GMM")) +
  
  geom_smooth(aes(y = hier_ari_scaled, color = "Hierarchical"), se = FALSE) +
  geom_point(aes(y = hier_ari_scaled, color = "Hierarchical")) +
  
  # --- SINGLE POINTS ---
  geom_point(aes(x = z_0, y = mfm_val, color = "BACON MFM"), size = 3) +
  geom_point(aes(x = z_0, y = lcesa_val, color = "LCESA"), size = 3) +
  
  # --- LABELS ---
  labs(x = "Cluster", y = "ARI", color = "Method") +
  scale_x_continuous(breaks = 1:max_k) +
  scale_y_continuous(limits = c(0, 1)) + 
  
  # --- MANUAL COLOR SELECTION ---
  scale_color_manual(
    # Change the hex codes below to pick your specific colors
    values = c(
      "BACON MFM"    = "#F8766D",  # Default Red/Coral
      "LCESA"        = "#B79F00",  # Default Olive
      "BACON FMM"    = "#00BA38",  # Default Green
      "k-Means"      = "#00BFC4",  # Default Teal
      "GMM"          = "#619CFF",  # Default Blue
      "Hierarchical" = "#F564E3"   # Default Pink
    ),
    # This controls the order in the legend
    breaks = c(
      "BACON MFM", 
      "LCESA", 
      "BACON FMM", 
      "k-Means", 
      "GMM", 
      "Hierarchical"
    )
  ) +
  theme_minimal()

plot
