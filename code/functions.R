#' Check whether a polygonal chain is closed
#' 
#' @author Kevin Jin
#'
#' @description
#' Determines if the given chain of coordinates is a polygonal chain or not
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#' @param reconstruct A boolean indicating whether the calling function is
#' the mapping function, which requires a rounding in the closed chain check.
#'
#' @return A logical value indicating whether the given chain is a polygonal
#' chain.
is_closed <- function(chain, reconstruct) {
  # Order of vertices must remain the same and might be disturbed after jitter.
  if (reconstruct) {
    # Also, you must allow some room for error because the reconstruct() 
    # functionand side length/angle proportion extraction processes might 
    # introduce error.
    if (identical(chain[1, ], round(chain[nrow(chain), ], digits = 10))) {
      is_closed <- TRUE
    } else {
      is_closed <- FALSE
    }
  } else {
    if (identical(chain[1, ], chain[nrow(chain), ])) {
      is_closed <- TRUE
    } else {
      is_closed <- FALSE
    }
  }
  return(is_closed)
}

#' Check whether a vertex produces a reflex angle in a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' Given the index of a vertex within a closed polygonal chain and the chain 
#' itself, close the chain without the vertex by drawing a third line and 
#' determine whether the given vertex is inside or outside of the resultant 
#' new polygonal chain.
#'
#' @param index A numeric containing the index of the vertex within the chain.
#' 
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#'
#' @return A logical value indicating whether the given vertex is in the
#' polygonal chain.
is_reflex <- function(index, chain) {
  require(sf)
  
  # Get offending vertex from the index
  point <- chain[index, ]
  
  # Cast the offending vertex to sf point object
  point_sf <- st_as_sf(data.frame(t(point)), coords = c("x", "y"))
  
  # Create a new polygonal chain without the offending vertex
  new_chain <- chain[-index, ]
  new_chain <- rbind(new_chain, new_chain[1, ]) # Duplicate first vertex
  
  # Cast the new polygonal chain to sf polygon object
  chain_sf <- st_sfc(st_polygon(list(new_chain)))
  
  # Determine whether offending vertex is inside the new polygonal chain
  is_reflex <- st_within(point_sf, chain_sf, sparse = FALSE)[, 1]
  
  return(is_reflex)
}

three_point_angle <- function(points) {
  pointA <- points[1, ]
  pointB <- points[2, ]
  pointC <- points[3, ]
  
  x1x2s <- (pointA[1] - pointB[1]) ^ 2
  x1x3s <- (pointA[1] - pointC[1]) ^ 2
  x2x3s <- (pointB[1] - pointC[1]) ^ 2
  
  y1y2s <- (pointA[2] - pointB[2]) ^ 2
  y1y3s <- (pointA[2] - pointC[2]) ^ 2
  y2y3s <- (pointB[2] - pointC[2]) ^ 2
  
  inner_ang <- (x1x2s + y1y2s + x2x3s + y2y3s - x1x3s - y1y3s) / (2 * sqrt(x1x2s + y1y2s) * sqrt(x2x3s + y2y3s))
  angle <- acos(inner_ang)
  
  if (is.na(angle)) {
    inner_ang <- min(max(inner_ang, -1), 1)
    angle <- ifelse(inner_ang == -1, pi, 0)
  }
  
  return(angle * (180 / pi))
}

get_interior_angles <- function(chain) {
  if (is_closed(chain, reconstruct = FALSE)) {
    # Number of vertices/rows is (k + 1) because the chain is closed
    n <- nrow(chain)
    
    # Loop over entire chain, checking each vertex for reflex angle
    reflex <- c()
    for (i in 1:n) {
      reflex <- c(reflex, is_reflex(i, chain))
    }
    
    # Create empty vector of length k to store interior angles
    angle <- rep(NA, n - 1)
    
    # Add penultimate vertex to first position for looping purposes
    a_chain <- rbind(chain[n - 1, ], chain)
    a_reflex <- append(reflex[n - 1], reflex)
    
    # Loop over entire chain, calculating the interior angles
    for (i in 2:n) {
      if (a_reflex[i]) { # Vertex produces a reflex angle; take the complement
        angle[i - 1] <- (360 - three_point_angle(a_chain[(i - 1):(i + 1), ]))
      } else { # Vertex does not produce a reflex angle
        angle[i - 1] <- three_point_angle(a_chain[(i - 1):(i + 1), ])
      }
    }
    
    # Normalize interior angles by the total to get relative interior angle
    angle <- compositional(angle, sum(angle))
    
  } else {
    stop("Argument is not a closed polygonal chain.")
  }
  return(angle)
}

#' Calculate the relative side lengths of a closed polygonal chain
#' 
#' @author Kevin Jin
#'
#' @description
#' If the given chain of coordinates is a closed polygonal chain, 
#' return a vector of its side lengths relative to the perimeter.
#'
#' @param chain A (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the polygonal chain.
#'
#' @return A vector of length k containing the side lengths of the polygonal 
#' chain.
get_side_lengths <- function(chain) {
  if (is_closed(chain, reconstruct = FALSE)) {
    # Eliminate repeated row for now
    chain <- chain[-nrow(chain), ]
    
    # Store side lengths in matrix
    side_lengths <- matrix(nrow = nrow(chain), ncol = 1)
    
    # Calculate the perimeter of the whole chain via Euclidean distance
    perimeter <- 0
    for (i in seq_len(nrow(chain))) {
      j <- i + 1
      if (i == (nrow(chain))) {
        j <- 1 # Loop back to first vertex once end of chain is reached
      }
      side_length <- sqrt(sum((chain[i, ] - chain[j, ]) ^ 2))
      side_lengths[i, ] <- side_length
      perimeter <- perimeter + side_length
    }
    
    # Normalize side lengths by the perimeter to get relative length
    side_lengths <- compositional(side_lengths, perimeter)
    
  } else {
    stop("Argument is not a closed polygonal chain.")
  }
  return(c(side_lengths))
}

#' Normalize a numerical vector by its total to produce its compositional data
#' 
#' @author Kevin Jin
#'
#' @description
#' Given a numerical vector, return a vector of its relative compositional data.
#'
#' @param data A numeric vector of length k containing data.
#' @param sum A numeric containing the sum of the vector.
#'
#' @return A vector of length k containing compositional data.
compositional <- function(data, sum) {
  normalized <- data / sum
  return(normalized)
}

#' Reconstruct all possible closed unit polygonal chains given relative interior 
#' angles and relative side lengths
#' 
#' @author Bryn Brakefield, Qiwei Li, Kevin Jin
#'
#' @description
#' Given two numerical vectors of compositional data containing relative 
#' interior angles and relative side lengths, reconstruct all possible 
#' unit polygonal chains and return the closed ones.
#'
#' @param angles A numeric vector of length n containing interior angle proportions.
#' @param side_lengths A numeric vector of length n containing side length proportions.
#'
#' @return A list or (k + 1) x 2 matrix containing the x-y coordinates of the 
#' vertices of the unit polygonal chain.
reconstruct <- function(angles, side_lengths) {
  # Get number of vertices
  n <- length(angles)
  
  # Convert compositional angle proportions back to radians
  angles <- angles * (n - 2) * pi
  
  # If the shape is a triangle, then ignore the provided side length
  # proportions and manually recalculate using the law of sines
  if (n == 3) {
    side_lengths <- sin(angles) / sum(sin(angles))
    side_lengths <- c(side_lengths[3], side_lengths[-3]) # Change ordering to be correct
  }
  
  pc0 <- expand.grid(replicate(n - 1, 0:1, simplify = FALSE))
  pc <- vector(mode = "list", length = 2^(n - 1))
  for (i in 1:(2 ^ (n - 1))) {
    V <- matrix(0, 2, n + 1)
    V[1, 2] <- side_lengths[1]
    theta <- rep(0, n + 1)
    for (j in 1:(n - 1)) {
      index <- j + 1
      if (pc0[i, j] == 1) {
        theta[index] <- theta[index - 1] + (pi - angles[index])
        V[, index + 1] <- c(side_lengths[index] * cos(theta[index]), 
                            side_lengths[index] * sin(theta[index])) + V[, index]
      }
      if (pc0[i, j] == 0) {
        theta[index] <- theta[index - 1] - (pi - angles[index])
        V[, index + 1] <- c(side_lengths[index] * cos(theta[index]), 
                            side_lengths[index] * sin(theta[index])) + V[, index]
      }
    }
    pc[[i]] <- V
  }
  
  # Create new list of closed chains (check is rounded to 10 digits)
  pc_closed <- list()
  for (i in 1:length(pc)) {
    if (is_closed(t(pc[[i]]), reconstruct = TRUE)) {
      pc_closed[[length(pc_closed) + 1]] <- list(t(pc[[1]]))
    } else {
      next
    }
  }
  
  return(pc_closed)
}

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

# -------------------------------------------------------------------------

# Function to calculate MAP estimate after burn-in
calc_MAP <- function(group_iter) {
  iters <- nrow(group_iter)
  post_burnin_iter <- group_iter[(floor(iters/2) + 1):iters, , drop = FALSE]  # Keep only post-burn-in samples
  
  # Compute the most frequent cluster assignment for each element
  apply(post_burnin_iter, 2, function(x) { 
    as.numeric(names(which.max(table(x))))
  })
}

# Function to calculate the Pairwise Probability Matrix (PPM) after burn-in
calc_PPM <- function(group_iter) {
  iters <- nrow(group_iter)
  post_burnin_iter <- group_iter[(floor(iters/2) + 1):iters, , drop = FALSE]  # Keep only post-burn-in samples
  U <- nrow(post_burnin_iter)  # Number of post-burn-in iterations
  n <- ncol(post_burnin_iter)  # Number of elements
  
  # Compute PPM as an n × n matrix
  PPM <- Reduce("+", lapply(1:U, function(u) {
    outer(post_burnin_iter[u, ], post_burnin_iter[u, ], FUN = "==")  # Compare elements within an iteration
  })) / U  # Average over all iterations
  
  return(PPM)
}

# Function to find the clustering that minimizes the sum of squared deviations (PPM Estimate)
calc_PPM_estimate <- function(group_iter, PPM) {
  iters <- nrow(group_iter)
  post_burnin_iter <- group_iter[(floor(iters/2) + 1):iters, , drop = FALSE]  # Keep only post-burn-in samples
  U <- nrow(post_burnin_iter)  # Number of post-burn-in iterations
  
  squared_deviations <- sapply(1:U, function(u) {
    assoc_matrix <- outer(post_burnin_iter[u, ], post_burnin_iter[u, ], FUN = "==")
    sum((assoc_matrix - PPM)^2)
  })
  
  best_iteration <- which.min(squared_deviations)
  return(post_burnin_iter[best_iteration, ])  # Return the best clustering from post-burn-in iterations
}

## Label Switching Function
switch_label <- function(z_store, K_final, Theta_store, Lambda_store){
  iter <- dim(z_store)[1]  # Number of valid iterations
  
  for (i in 1:iter){
    z_tmp <- z_store[i, ] + 1  # Convert to 1-based indexing
    Theta_tmp <- Theta_store[[i]]
    Lambda_tmp <- Lambda_store[[i]]
    
    order_tmp <- order(Theta_tmp[, 1])  # Order clusters based on the first column
    
    # Apply ordering to Theta and Lambda first
    Theta_store[[i]] <- Theta_tmp[order_tmp, ]
    Lambda_store[[i]] <- Lambda_tmp[order_tmp, ]
    
    # Now, relabel clusters in z_store after applying ordering
    label_map <- setNames(seq_along(order_tmp), order_tmp)  # Old label → New label mapping
    z_store[i, ] <- sapply(z_tmp, function(x) label_map[as.character(x)])
  }
  
  return(list("z_store" = z_store, "Theta_store" = Theta_store, "Lambda_store" = Lambda_store))
}

mean_k <- function(theta) {
  K <- dim(theta)[1]
  n <- dim(theta)[2]
  
  mean_theta_k <- array(dim = c(1, n, K))
  for (i in 1:K) {
    theta_sum_k <- sum(theta[i, ])
    mean_theta_k[, , i] <- theta[i, ]/theta_sum_k
  }
  return(mean_theta_k)
}

var_cov_k <- function(theta) {
  K <- dim(theta)[1]
  n <- dim(theta)[2]
  
  cov_mat_theta_k <- array(dim = c(n, n, K))
  for (i in 1:K) {
    theta_sum_k <- sum(theta[i, ])
    
    var_theta_k <- (theta[i, ] * (theta_sum_k - theta[i, ])) / ((theta_sum_k^2) * (theta_sum_k + 1))
    
    theta_k_outer <- outer(theta[i, ], theta[i, ])
    
    cov_mat_theta_k[, , i] <- -theta_k_outer / ((theta_sum_k^2) * (theta_sum_k + 1))
    
    diag(cov_mat_theta_k[, , i]) <- var_theta_k
  }
  return(cov_mat_theta_k)
}
