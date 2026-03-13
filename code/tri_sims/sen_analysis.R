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

# -------------------------------------------------------------------------

K_0 <- 10

K_max <- 30
iter <- 1000
burn <- iter/2

w <- c(1, 10, 50, 100, 500)

w_a <- c(1)
w_l <- c(0, 0.25, 0.5, 0.75, 1)
reps <- 30
combin <- expand.grid(rep = 1:reps, w_l = w_l, w_a = w_a, w = w)

total_rows <- nrow(combin)
ari_sims_map <- numeric(total_rows)
ind <- 0
for (i in 1:length(w)) {
  A <- A2_tensor[, , i]
  L <- L2_tensor[, , i]
  
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
  
  for (j in 1:length(w_l)) {
    for (k in 1:reps) {
      ind <- ind + 1
      
      res <- bacon_mfm(A, L, w_A = combin$w_a[ind], w_L = combin$w_l[ind], 
                       z_0, K_max = K_max, alpha = 30, 
                       est_sr = TRUE, est_s = TRUE, est_r = TRUE, alpha_s = 1, beta_s = 1, omega = 0.5, 
                       a_theta = 0.001, b_theta = 0.001, a_lambda = 0.001, b_lambda = 0.001, tau_theta = 0.1, tau_lambda = 0.1, 
                       iter = iter)
      
      z_store <- res$z_store
      max_K <- max(res$K_store)
      prob_z <- matrix(0, m, max_K)
      for(x in 1:m) {
        for(y in 1:max_K) {
          prob_z[x, y] <- mean((z_store[(burn + 1):iter, x] + 1) == y)
        }
      }
      z_map <- as.vector(apply(prob_z, 1, function(x) which(x == max(x))[1]))
      
      ari_sims_map[ind] <- adjustedRandIndex(z_map, z)
      
      rm(res)
    }
  }
  
  print(i)
}

res_sims <- cbind(combin, ari_sims_map)

# -------------------------------------------------------------------------

w <- c(1, 10, 50, 100, 500)

w_l <- c(1)
w_a <- c(0, 0.25, 0.5, 0.75, 1)
reps <- 30
combin <- expand.grid(rep = 1:reps, w_a = w_a, w_l = w_l, w = w)

total_rows <- nrow(combin)
ari_sims_map <- numeric(total_rows)
ind <- 0
for (i in 1:length(w)) {
  A <- A2_tensor[, , i]
  L <- L2_tensor[, , i]
  
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
  #adjustedRandIndex(z, init_labels)
  z_0 <- init_labels - 1
  
  for (j in 1:length(w_a)) {
    for (k in 1:reps) {
      ind <- ind + 1
      
      res <- bacon_mfm(A, L, w_A = combin$w_a[ind], w_L = combin$w_l[ind], 
                       z_0, K_max = K_max, alpha = 30, 
                       est_sr = TRUE, est_s = TRUE, est_r = TRUE, alpha_s = 1, beta_s = 1, omega = 0.5, 
                       a_theta = 0.001, b_theta = 0.001, a_lambda = 0.001, b_lambda = 0.001, tau_theta = 0.1, tau_lambda = 0.1, 
                       iter = iter)
      
      z_store <- res$z_store
      max_K <- max(res$K_store)
      prob_z <- matrix(0, m, max_K)
      for(x in 1:m) {
        for(y in 1:max_K) {
          prob_z[x, y] <- mean((z_store[(burn + 1):iter, x] + 1) == y)
        }
      }
      z_map <- as.vector(apply(prob_z, 1, function(x) which(x == max(x))[1]))
      
      ari_sims_map[ind] <- adjustedRandIndex(z_map, z)
      
      rm(res)
    }
  }
  
  print(i)
}

res_sims <- cbind(combin, ari_sims_map)

# -------------------------------------------------------------------------

load("./res/tri_sims/wl.RData")

orange_shades <- rev(brewer.pal(n = 5, name = "Oranges"))
names(orange_shades) <- c("0", "0.25", "0.5", "0.75", "1")

w_a_labels <- c("0", "0.25", "0.5", "0.75", "1")

# Make w and w_a factors
res_sims$w <- factor(res_sims$w) # x-axis categories

res_sims$w_a <- factor(
  res_sims$w_a,
  levels = c(0, 0.25, 0.5, 0.75, 1),
  labels = c("0", "0.25", "0.5", "0.75", "1")
)

p <- ggplot(res_sims, aes(x = w, y = ari_sims_map, fill = w_a)) +
  geom_boxplot(
    width = 0.6,
    position = position_dodge(width = 0.9),
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = orange_shades,
    labels = w_a_labels,
    name = expression(w[a])
  ) + coord_cartesian(ylim = c(0, 1)) + # Added this line to zoom the y-axis
  labs(
    x = "c",
    y = "ARI"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_blank(),
    legend.position = "right"
  )

p

# -------------------------------------------------------------------------

# 1. LOAD AND PATCH DATA --------------------------------------------------

load("./res/tri_sims/wl.RData")

unstable_values <- res_sims %>%
  mutate(w = as.numeric(as.character(w))) %>% # Force 'w' to be numeric for joining
  filter(w_a == 1, w_l == 1) %>%
  select(rep, w, ari_sims_map) 

load("./res/tri_sims/wa.RData")

# Patch res_sims
res_sims <- res_sims %>%
  mutate(w = as.numeric(as.character(w))) %>% # Force 'w' to be numeric here too
  left_join(unstable_values, by = c("rep", "w"), suffix = c("", "_new")) %>%
  mutate(ari_sims_map = ifelse(w_a == 1 & w_l == 1, ari_sims_map_new, ari_sims_map)) %>%
  select(-ari_sims_map_new)

# 2. DATA PREPARATION FOR PLOTTING -----------------------------------------

# Now that the join is done, convert 'w' to a factor for the x-axis labels
res_sims$w <- factor(res_sims$w, levels = c(1, 10, 50, 100, 500)) 

# Convert w_l to factor for the fill colors
res_sims$w_l <- factor(res_sims$w_l, levels = c(0, 0.25, 0.5, 0.75, 1))

# Define Green shades
green_shades <- rev(brewer.pal(n = 5, name = "Greens"))
names(green_shades) <- c("0", "0.25", "0.5", "0.75", "1")

# 3. PLOT GENERATION -------------------------------------------------------

p <- ggplot(res_sims, aes(x = w, y = ari_sims_map, fill = w_l)) +
  geom_boxplot(
    width = 0.6,
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  scale_fill_manual(
    values = green_shades,
    name = expression(w[l])
  ) + 
  scale_y_continuous(
    # Horizontal numbered lines (0, 0.25, 0.5, 0.75, 1)
    breaks = seq(0, 1, by = 0.25),      
    # One horizontal unnumbered line exactly in the middle
    minor_breaks = seq(0.125, 0.875, by = 0.25) 
  ) + 
  coord_cartesian(ylim = c(0, 1)) + 
  labs(
    x = "c",
    y = "ARI"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    # HORIZONTAL GRID: Major is numbered, minor is the "in-between" line
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.minor.y = element_line(color = "gray95"), 
    
    # VERTICAL GRID: Show only the major line for each 'c' value
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor.x = element_blank() # Removes any extra vertical lines
  )

p
