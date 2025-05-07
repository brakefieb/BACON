library(ggplot2)
library(dplyr)
library(RColorBrewer)

# -------------------------------------------------------------------------

w_l <- c(1)
w_a <- c(0, 0.5, 1)
reps <- 10
combin <- expand.grid(w_l = w_l, rep = 1:reps, w_a = w_a, w = w)

ari_map_sims <- c()
ind <- 0
for (i in 1:length(w)) {
  A <- A2_tensor[, , i]
  L <- L2_tensor[, , i]
  
  for (j in 1:length(w_a)) {
    for (k in 1:reps) {
      ind <- ind + 1
      
      res <- bacon_mfm(A = A, L = L, w_a = combin[ind, 3], w_l = combin[ind, 1])
      iter_mat <- res$z_store
      
      # Compute MAP and PPM estimates after burn-in
      map_estimate <- calc_MAP(iter_mat)  # MAP estimate
      #PPM <- calc_PPM(iter_mat)  # Pairwise Probability Matrix
      #ppm_estimate <- calc_PPM_estimate(iter_mat, PPM)  # Best clustering based on PPM
      
      # Compute Adjusted Rand Index (ARI) for both estimates
      ari_map_sims[ind] <- adjustedRandIndex(map_estimate, z)
      #ari_ppm <- adjustedRandIndex(ppm_estimate, z)
    }
  }
  
  print(i)
}

res_sims <- cbind(combin, ari_map_sims)

# Green shades using RColorBrewer
green_shades <- brewer.pal(n = 3, name = "Greens")
names(green_shades) <- c("0", "0.5", "1")
w_a_labels <- c("0", "0.5", "1")

# Plot
p <- ggplot(plot_data, aes(x = w, y = ARI, fill = w_a)) +
  geom_boxplot(position = position_dodge(width = 0.9), width = 0.6) +
  scale_fill_manual(values = green_shades,
                    labels = w_a_labels,
                    name = expression(w[a])) +
  labs(x = "c", y = "ARI") +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    legend.position = "right"
  )

# Save
ggsave(p, file = "mfm_sen.png", width = 10, height = 5, units = "in", dpi = 300)

# -------------------------------------------------------------------------

# Specify
n_0 <- 71

# Combine
r_j <- c(rep(0, n_0), rep(1, n_0))
pi_rj <- c(pi_rj_0, pi_rj_1)
index <- 1:(n_0 + n_0)
df <- data.frame(index = index, r_j = r_j, pi_rj = pi_rj)

# Plot
p1 <- ggplot(df, aes(x = index, y = pi_rj)) +
  geom_segment(aes(xend = index, yend = 0), color = "black") +
  geom_point(color = "black") +
  scale_x_continuous(breaks = df$index, labels = df$r_j) +
  labs(
    x = expression(r[j]),
    y = bquote(pi(r[j] ~ "|" ~ "\u00B7"))  
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(p1, file = "post_r.png", width = 15, height = 5, units = "in", dpi = 300)

# -------------------------------------------------------------------------

#Specify
n_s <- 23

# Combine labels and posterior values
s_j <- c(rep(0, n_s), rep(1, n_s), rep(2, n_s))
pi_sj <- c(pi_sj_0, pi_sj_1, pi_sj_2)
index <- 1:(3 * n_s)
df <- data.frame(index = index, s_j = s_j, pi_sj = pi_sj)

# Plot stem-style posterior plot
p2 <- ggplot(df, aes(x = index, y = pi_sj)) +
  geom_segment(aes(xend = index, yend = 0), color = "black") +
  geom_point(color = "black") +
  scale_x_continuous(breaks = df$index, labels = df$s_j) +
  labs(
    x = expression(s[j]),
    y = bquote(pi(s[j] ~ "|" ~ "\u00B7"))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save the plot
ggsave(p2, file = "post_s.png", width = 10, height = 5, units = "in", dpi = 300)

# -------------------------------------------------------------------------

# Set shared color scale limits
common_limits <- range(unlist(lambda_var_list))

# Final plot
p <- ggplot() +
  geom_polygon(data = combined_polygon_data, aes(x = x, y = y, group = Group), fill = NA) +
  geom_segment(data = combined_side_data,
               aes(x = x, y = y, xend = xend, yend = yend, color = Variance),
               size = 1.5) +
  scale_color_gradient(low = "blue", high = "red", limits = common_limits, name = "Variance") +
  coord_fixed(ratio = 1) +
  facet_grid(. ~ Group) +
  scale_x_continuous(name = "x") +
  scale_y_continuous(name = "y") +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.title = element_text(),
    axis.text = element_text(),
    axis.ticks = element_line()
  )

# Save the plot
ggsave("mean_shape.png", p, width = 10, height = 3, units = "in", dpi = 300)
