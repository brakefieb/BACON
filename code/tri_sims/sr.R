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

s <- s_matrix[, ind]
r <- r_matrix[, ind]

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
                 est_sr = FALSE, est_s = TRUE, est_r = TRUE, alpha_s = 1, beta_s = 1, omega = 0.5, 
                 a_theta = 0.001, b_theta = 0.001, a_lambda = 0.001, b_lambda = 0.001, tau_theta = 0.1, tau_lambda = 0.1, 
                 iter = iter)

# -------------------------------------------------------------------------

load("./res/tri_sims/res.RData")

# -------------------------------------------------------------------------

prob_s <- matrix(0, m, n)
for(j in 2:m) {
  for(s_0 in 1:n) {
    prob_s[j, s_0] <- mean((res$s_store[(burn + 1):iter, j] + 1) == s_0)
  }
}
s_map <- as.vector(apply(prob_s[-1, ], 1, function(x) which(x == max(x))[1])) - 1

adjustedRandIndex(s_map[(sum(z == 0) + 1):(m - 1)], s[52:150])

prob_r <- matrix(0, m, 2)
for(j in 2:m) {
  for(r_0 in 1:2) {
    prob_r[j, r_0] <- mean((res$r_store[(burn + 1):iter, j] + 1) == r_0)
  }
}
r_map <- as.vector(apply(prob_r[-1, ], 1, function(x) which(x == max(x))[1])) - 1

adjustedRandIndex(r_map[(sum(z == 0) + 1):(m - 1)], r[52:150])

# -------------------------------------------------------------------------

s_burn <- res$prob_s_store[(sum(z == 0) + 2):m, , (burn + 1):iter]

s_new <- s_map[(sum(z == 0) + 1):(m - 1)]

s_probs_all_iters <- matrix(0, 99, burn)
for(j in 1:99) {
  state_index <- s_new[j] + 1
  s_probs_all_iters[j, ] <- s_burn[j, state_index, ]
}
mean_probs <- rowMeans(s_probs_all_iters)

plot_data <- data.frame(
  obs_index = 1:99,
  state = s[52:150],
  prob = mean_probs
)
plot_data <- plot_data[order(plot_data$state), ]
plot_data$plot_x <- 1:99 
plot_data$state_factor <- factor(plot_data$state, 
                                 levels = sort(unique(plot_data$state)))
plot_data <- plot_data %>%
  group_by(state_factor) %>%
  mutate(local_x = row_number()) %>%
  ungroup()

# 1. Create the palette (5 shades)
green_palette <- brewer.pal(n = 5, name = "Greens")
# 2. Pick the 3 darkest shades (indices 3, 4, 5) so they show up well
selected_greens <- green_palette[3:5]

# 3. Plot
p <- ggplot(plot_data, aes(x = local_x, y = prob)) +
  geom_point(aes(color = state_factor), alpha = 0.7, size = 2) +
  geom_hline(yintercept = 1/3, linetype = "dashed") +
  facet_wrap(~state_factor, scales = "free_x") + 
  
  theme_bw() + 
  
  # Apply the 3 specific green colors here
  scale_color_manual(values = selected_greens) +
  
  scale_x_continuous(breaks = seq(0, max(plot_data$local_x), by = 10)) +
  
  labs(
    # The fixed expression with big parentheses
    y = expression(bar(pi) * group("(", hat(s)[j]^MAP == s[j] * " | " * "·", ")")), 
    color = expression(s[j]),
    x = ""
  ) +
  ylim(0, 1) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p

# -------------------------------------------------------------------------

r_prob_burn <- res$prob_r_store[(burn + 1):iter, -1] 
r_prob_burn <- t(r_prob_burn)  

r_subset <- r[-1]

m_new <- m - 1
n_post_burn <- burn
r_probs_all_iters <- matrix(0, m_new, n_post_burn)
for(j in 1:m_new) {
  if(r_subset[j] == 1) {
    # If assigned state is 1, use the stored probability
    r_probs_all_iters[j, ] <- r_prob_burn[j, ]
  } else {
    # If assigned state is 0, use 1 - P(r=1)
    r_probs_all_iters[j, ] <- 1 - r_prob_burn[j, ]
  }
}
mean_probs <- rowMeans(r_probs_all_iters)
mean_probs <- mean_probs[(sum(z == 0) + 1):(m - 1)]

plot_data <- data.frame(
  obs_index = 1:99,
  state = r[52:150],
  prob = mean_probs
)
plot_data <- plot_data %>%
  arrange(state) %>%
  group_by(state) %>%
  mutate(local_x = row_number()) %>%
  ungroup() %>%
  mutate(state_factor = factor(state, levels = c(0, 1), labels = c("0", "1")))

# 1. Load the full palette (9 shades)
orange_palette <- brewer.pal(n = 9, name = "Oranges")
# 2. Pick two distinct colors (e.g., Index 5 is medium, Index 9 is very dark)
# This ensures both are visible and distinguishable.
two_oranges <- orange_palette[c(5, 9)]

# 3. Plot
p2 <- ggplot(plot_data, aes(x = local_x, y = prob)) +
  geom_point(aes(color = state_factor), alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  facet_wrap(~state_factor, scales = "free_x") + 
  
  theme_bw() + 
  
  # Apply the 2 specific green colors here
  scale_color_manual(values = two_oranges) +
  
  scale_x_continuous(breaks = seq(0, max(plot_data$local_x), by = 10)) +
  
  labs(
    # Updated expression for r
    y = expression(bar(pi) * group("(", hat(r)[j]^MAP == r[j] * " | " * "·", ")")), 
    color = expression(r[j]),
    x = ""
  ) +
  ylim(0, 1) +
  
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p2
