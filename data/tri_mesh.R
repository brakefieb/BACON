#https://math.ntnu.edu.tw/~yueh/courses/MATLAB_MeshDataStructure.html

# -------------------------------------------------------------------------

library(R.matlab)

# -------------------------------------------------------------------------

mat_file <- "E:/mesh/tri_mesh/Triangle2D.mat"
data <- readMat(mat_file)
Triangle2D <- data$Triangle2D 

num_triangles <- dim(Triangle2D)[1]
A <- matrix(0, nrow=num_triangles, ncol=3) 
L <- matrix(0, nrow=num_triangles, ncol=3) 
for (i in 1:num_triangles) {
  chain <- Triangle2D[i,,] 
  chain <- matrix(chain, ncol=2, byrow=FALSE)
  colnames(chain) <- c("x", "y")
  
  angles <- get_interior_angles(chain)
  side_lengths <- get_side_lengths(chain)

  A[i, ] <- angles     
  L[i, ] <- side_lengths  
}

# -------------------------------------------------------------------------

save(A, L, file = "E:/mesh/tri_mesh/tri_mesh.RData")
