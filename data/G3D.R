#https://ieeexplore.ieee.org/document/6239175

# -------------------------------------------------------------------------

library(dplyr)
library(stringr)
library(R.matlab)

# -------------------------------------------------------------------------

base_path <- "E:/G3D"
classes <- c("KickLeft", "KickRight", "PunchLeft", "PunchRight")
class_labels <- c("KickLeft" = 0, "KickRight" = 1, "PunchLeft" = 2, "PunchRight" = 3)

angles_matrix <- list()
sides_matrix <- list()
z <- c()
for (class_name in classes) {
  class_path <- file.path(base_path, class_name)
  mat_files <- list.files(class_path, pattern = "\\.mat$", full.names = TRUE, recursive = TRUE)
  
  for (mat_file in mat_files) {
    data <- readMat(mat_file)
    
    chain <- data$boundaryPoints
    chain <- rbind(chain, chain[1, ]) 
    colnames(chain) <- c("x", "y")
    
    angles <- get_interior_angles(chain)
    side_lengths <- get_side_lengths(chain)
    
    angles_matrix <- append(angles_matrix, list(angles))
    sides_matrix <- append(sides_matrix, list(side_lengths))
    z <- c(z, class_labels[class_name])
  }
}

A <- do.call(rbind, angles_matrix)
L <- do.call(rbind, sides_matrix)
z <- unname(z)

# -------------------------------------------------------------------------

save_path <- "E:/G3D/"
save(A, L, z, file = file.path(save_path, "G3D.RData"))
