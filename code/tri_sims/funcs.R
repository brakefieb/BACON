library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(RcppEigen)

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
    angle <- runif(1, 0, 360) 
    clockwise <- sample(c(TRUE, FALSE), 1)
    chain <- rotate(chain, angle, clockwise)
  }
  if (apply_translate) {
    x_shift <- runif(1, -.5, .5)  
    y_shift <- runif(1, -.5, .5)  
    chain <- translate(chain, x_shift, y_shift)
  }
  if (apply_dilate) {
    factor <- runif(1, 1, 1.5)  
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
