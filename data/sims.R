library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(RcppEigen)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# -------------------------------------------------------------------------

sourceCpp(code = '
  #include <RcppArmadillo.h>
  
  // [[Rcpp::depends(RcppArmadillo)]]
  using namespace arma;
  
  // [[Rcpp::export]]
  double qtbeta1(double p, double alpha, double beta, double a, double b) {
    double Fa = R::pbeta(a, alpha, beta, true, false);
    double Fb = R::pbeta(b, alpha, beta, true, false);
    double u = Fa + p*(Fb - Fa);
    double quan = R::qbeta(u, alpha, beta, true, false);
    return quan;
  }
  
  // [[Rcpp::export]]
  mat rtdirichlet1(vec eta, vec a, vec b) {
    int k = eta.size();
    vec u = randu(k);
    vec x = zeros(k);
    x[k - 2] = qtbeta1(u[k - 2], eta[k - 2], sum(eta) - eta[k - 2], std::max(a[k - 2], 1 - sum(b) + b[k - 2]), std::min(b[k - 2], 1 - sum(a) + a[k - 2]));
    for (int i = k - 3; i >= 0; i--) {
      double a0 = std::max(a[i]/(1 - sum(x(span(i + 1, k - 2)))), 1 - (sum(b) - sum(b(span(i, k - 2))))/(1 - sum(x(span(i + 1, k - 2)))));
      double b0 = std::min(b[i]/(1 - sum(x(span(i + 1, k - 2)))), 1 - (sum(a) - sum(a(span(i, k - 2))))/(1 - sum(x(span(i + 1, k - 2)))));
      x[i] = (1 - sum(x(span(i + 1, k - 2))))*qtbeta1(u[i], eta[i], sum(eta) - sum(eta(span(i, k - 2))), a0, b0);
    }
    x[k - 1] = 1 - sum(x);
    return x;
  }')
rtdirichlet <- function(n, eta, a, b) {
  return(t(replicate(n, as.vector(rtdirichlet1(eta, a, b)))))
}

sourceCpp(code = '
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
#include <RcppEigen.h>
  
// [[Rcpp::depends(RcppArmadillo, RcppDist, RcppEigen)]]
  
using namespace Rcpp;
using namespace arma;
  
// [[Rcpp::export]]
Rcpp::List L2A2(arma::mat L, arma::mat A, IntegerVector s, IntegerVector r) {

// Read data information
int m = L.n_rows;
int n = L.n_cols;

 // Rearranged length proportion and angle proportions L2 and A2
NumericMatrix L2(m, n);
NumericMatrix A2(m, n);

for (int j = 0; j < m; j++){
  for (int i = 0; i < n; i++){
    for (int ii = 0; ii < n; ii++){
      int temp = ((s(j) + i*(1-2*r(j))) + n) % n;  // in case negative value, -1%4=-1, -1%4=3 in R
      if (temp == ii){
        L2(j, i) = L(j, ii);
        A2(j, i) = A(j, ii);
      }
    }
  }
}

return Rcpp::List::create(Rcpp::Named("L2") = L2, Rcpp::Named("A2") = A2);
}
')

v <- function(a, l, r) {
  n <- length(a)
  
  V <- matrix(0, 2, n + 1)
  V[1, 2] <- l[1]
  theta <- rep(0, n + 1)
  if (r == 1) {
    for (j in 1:(n - 1)) {
      index <- j + 1
      theta[index] <- theta[index - 1] + (pi - a[index])
      V[, index + 1] <- c(l[index]*cos(theta[index]), l[index]*sin(theta[index])) + V[, index]
    } 
  }
  if (r == 0) {
    for (j in 1:(n - 1)) {
      index <- j + 1
      theta[index] <- theta[index - 1] - (pi - a[index])
      V[, index + 1] <- c(l[index]*cos(theta[index]), l[index]*sin(theta[index])) + V[, index]
    }
  }
  return(V)
}

#' Dilate a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Scales the size of a polygonal chain to be greater or smaller, centered
#' about its centroid.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#' @param factor Positive floating-point number to scale the chain by.
#' A factor > 1 dilates the chain, while a factor > 0 and < 1 shrinks the chain.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the vertices of
#' the dilated chain.
dilate <- function(chain, factor = 1) {
  if (factor <= 0) {
    stop("Dilation factor must be greater than 0.\n")
  } else {
    # Calculate the centroid vector
    centroid <- t(matrix(rowSums(t(chain)) / nrow(chain))[, rep(1, each = nrow(chain))])
    chain <- (chain - centroid) * factor + centroid
  }
  return(chain)
}

#' Reflect a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Reflects a polygonal chain vertically or horizontally.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#' @param direction String containing direction in which to reflect the chain.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the vertices of
#' the reflected chain.
#'
reflect <- function(chain, direction = c("horizontal", "vertical")) {
  centroid <- t(matrix(rowSums(t(chain)) / nrow(chain))[, rep(1, each = nrow(chain))])
  if (direction == "horizontal") {
    # Reflect across the y-axis
    chain <- (chain - centroid) %*%
      matrix(c(-1, 0, 0, 1), ncol = 2, byrow = TRUE) + centroid
    
    # For the sf function
    colnames(chain) <- c("x", "y")
    
  } else if (direction == "vertical") {
    # Reflect across the x-axis
    chain <- (chain - centroid) %*%
      matrix(c(1, 0, 0, -1), ncol = 2, byrow = TRUE) + centroid
    
    # For the sf function
    colnames(chain) <- c("x", "y")
    
  } else {
    stop("Invalid or no direction provided.\n")
  }
  return(chain)
}

#' Rotate a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Rotates a polygonal chain by a specified angle about its centroid.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#' @param angle Rotation angle in degrees.
#' @param clockwise Rotate the chain clockwise if true, counterclockwise if false.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the vertices of
#' the rotated chain.
#'
rotate <- function(chain, angle, clockwise = TRUE) {
  # Convert argument angle to radians, as R's trigonometric functions use radians
  angle <- angle * (pi / 180)
  
  # Eliminate repeated row for now
  chain <- chain[-nrow(chain), ]
  
  # Calculate centroid vector
  centroid <- matrix(rowSums(t(chain)) / nrow(chain))[, rep(1, each = nrow(chain))]
  if (clockwise) {
    # Calculate counterclockwise rotation matrix
    rotation <- matrix(c(cos(angle), -sin(angle),
                         sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    
    # Rotate chain in place
    chain <- rotation %*% (t(chain) - centroid) + centroid
    
  } else {
    # Calculate clockwise rotation matrix
    rotation <- matrix(c(cos(angle), sin(angle),
                         -sin(angle), cos(angle)), ncol = 2, byrow = TRUE)
    
    # Rotate chain in place
    chain <- rotation %*% (t(chain) - centroid) + centroid
  }
  
  # Add repeated row back
  chain <- t(chain)
  chain <- rbind(chain, chain[1, ])
  
  # For the sf function
  colnames(chain) <- c("x", "y")
  
  # Transpose chain to original form
  return(chain)
}

#' Translate a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Moves every point of a polygonal chain by the same distance in a given
#' direction.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the vertices
#' of the polygonal chain.
#' @param x Horizontal distance to translate the chain by.
#' @param y Vertical distance to translate the chain by.
#'
#' @return A (k + 1) x 2 matrix containing the x-y coordinates of the vertices of
#' the translated chain.
translate <- function(chain, x = 0, y = 0) {
  # Horizontal translation
  chain[, "x"] <- chain[, "x"] + x
  # Vertical translation
  chain[, "y"] <- chain[, "y"] + y
  return(chain)
}

# Function to randomly transform a polygonal chain
random_transform <- function(chain) {
  # Randomly decide whether to apply each transformation
  apply_reflect <- sample(c(TRUE, FALSE), 1)
  apply_rotate <- sample(c(TRUE, FALSE), 1)
  apply_translate <- sample(c(TRUE, FALSE), 1)
  apply_dilate <- sample(TRUE, 1)
  
  # Apply transformations with random parameters
  if (apply_reflect) {
    direction <- sample(c("horizontal", "vertical"), 1)
    chain <- reflect(chain, direction)
  }
  if (apply_rotate) {
    angle <- runif(1, 0, 360) # Random angle between 0 and 360 degrees
    clockwise <- sample(c(TRUE, FALSE), 1)
    chain <- rotate(chain, angle, clockwise)
  }
  if (apply_translate) {
    x_shift <- runif(1, -.5, .5)  
    y_shift <- runif(1, -.5, .5)  
    chain <- translate(chain, x_shift, y_shift)
  }
  if (apply_dilate) {
    factor <- runif(1, 1, 1.5)  # Random scaling factor between 0.5 and 2
    chain <- dilate(chain, factor)
  }
  
  return(chain)
}

# Provided function to determine the new vertex order based on s and r
get_vertex_labels <- function(s, r) {
  # s: starting vertex (0, 1, or 2) in 0-indexing
  # r: orientation flag; 0 means keep order, 1 means reverse order of the other vertices
  ref_vertices <- c(0, 1, 2)
  start_pos <- which(ref_vertices == s)
  
  if (r == 0) {  # Keep clockwise ordering: simply rotate starting from s
    result <- c()
    for (i in 0:2) {
      pos <- ((start_pos - 1 + i) %% 3) + 1
      result[i + 1] <- ref_vertices[pos]
    }
  } else {  # Reverse order for the other vertices (counterclockwise)
    result <- c()
    result[1] <- s
    next_pos <- start_pos %% 3 + 1
    last_pos <- next_pos %% 3 + 1
    result[2] <- ref_vertices[last_pos]
    result[3] <- ref_vertices[next_pos]
  }
  
  return(result)
}

# -------------------------------------------------------------------------

m_z <- rep(50, 3)
m <- sum(m_z)
z <- c(rep(0, m_z[1]), rep(1, m_z[2]), rep(2, m_z[3]))
K <- 3
n <- 3
a_lower <- rep(0, n)
b_upper <- rep(1, n)

w <- c(1, 5, 10, 50, 100)

A_tensor <- array(0, dim = c(m, n, length(w)))
L_tensor <- array(0, dim = c(m, n, length(w)))
for (i in seq_along(w)) {
  eta_1 <- rep(1, n)*w[i]
  eta_2 <- c(3, 9, 30)*w[i]
  eta_3 <- c(15, 10, 5)*w[i]
  
  a_1 <- rtdirichlet(n = m_z[1], eta = eta_1, a = a_lower, b = b_upper)
  a_2 <- rtdirichlet(n = m_z[2], eta = eta_2, a = a_lower, b = b_upper)
  a_3 <- rtdirichlet(n = m_z[3], eta = eta_3, a = a_lower, b = b_upper)
  A <- rbind(a_1, a_2, a_3)
  A_tensor[, , i] <- A
  
  A_rad <- A*(n - 2)*pi
  
  L <- matrix(nrow = m, ncol = n)
  for (j in 1:m) {
    angles <- A_rad[j, ]
    side_lengths <- sin(angles) / sum(sin(angles))
    side_lengths <- c(side_lengths[3], side_lengths[-3]) 
    L[j, ] <- side_lengths/sum(side_lengths) # Sanity check for rounding error
  }
  L_tensor[, , i] <- L
  
  print(i)
}

save(A_tensor, L_tensor, z, K, m, n, file = "triangle_sims.RData")

# -------------------------------------------------------------------------

s <- 0:(n - 1)
r <- 0:1

A2_tensor <- array(0, dim = c(m, n, length(w)))
L2_tensor <- array(0, dim = c(m, n, length(w)))
for (i in seq_along(w)) {
  s <- c(0, sample(s, m - 1, replace = TRUE))
  r <- c(0, sample(r, m - 1, replace = TRUE))
  
  L2_A2 <- L2A2(L_tensor[, , i], A_tensor[, , i], s, r)
  
  A2_tensor[, , i] <- L2_A2$A2
  L2_tensor[, , i] <- L2_A2$L2
  
  print(i)
}

save(A2_tensor, L2_tensor, z, K, m, n, s, r, file = "triangle_sims_sr.RData")

# -------------------------------------------------------------------------

#Specify
ind <- 1

A <- A_tensor[, , ind]
L <- L_tensor[, , ind]

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

save(transformed_shapes, file = "transformed_shapes.RData")

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
