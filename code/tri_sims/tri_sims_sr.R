source("./code/tri_sims/funcs.R")

# -------------------------------------------------------------------------

m_z <- rep(50, 3)
m <- sum(m_z)
n <- 3

a_lower <- rep(0, n)
b_upper <- rep(1, n)

w <- c(1, 10, 50, 100, 500)

A_tensor <- array(0, dim = c(m, n, length(w)))
L_tensor <- array(0, dim = c(m, n, length(w)))
for (i in seq_along(w)) {
  eta_1 <- rep(10, n)*w[i]
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
    L[j, ] <- side_lengths/sum(side_lengths) 
  }
  L_tensor[, , i] <- L
  
  print(i)
}

z <- c(rep(0, m_z[1]), rep(1, m_z[2]), rep(2, m_z[3]))
K <- 3

# -------------------------------------------------------------------------

s_space <- 0:(n - 1)
r_space <- 0:1

A2_tensor <- array(0, dim = c(m, n, length(w)))
L2_tensor <- array(0, dim = c(m, n, length(w)))
s_matrix <- matrix(0, nrow = m, ncol = length(w))
r_matrix <- matrix(0, nrow = m, ncol = length(w))
for (i in seq_along(w)) {
  s_iter <- c(0, sample(s_space, m - 1, replace = TRUE))
  r_iter <- c(0, sample(r_space, m - 1, replace = TRUE))
  
  # Store the tracking vectors into the i-th column
  s_matrix[, i] <- s_iter
  r_matrix[, i] <- r_iter
  
  # Using the cleaner L2A2_try logic from earlier (assuming you renamed it L2A2)
  L2_A2 <- L2A2(L_tensor[, , i], A_tensor[, , i], s_iter, r_iter)
  
  A2_tensor[, , i] <- L2_A2$A2
  L2_tensor[, , i] <- L2_A2$L2
  
  print(i)
}
