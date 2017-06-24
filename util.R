#==============================================================================
# Utility functions
#==============================================================================
inner_product <- function(v0, v1) {sum(v0 * v1)}

euclidean_norm <- function(v0) {sqrt(inner_product(v0, v0))}

a_project_onto_b <- function(a, b) {
  unit_b <- b / euclidean_norm(b)
  inner_product(a, unit_b)
}

center <- function(X) {
  p <- ncol(X)
  for (i in 1:p) {
    X[,i] <- X[,i] - mean(X[,i])
  }
  X
}

get_pc_variance <- function(W, X) {
  pc_var <- apply((X %*% W)^2, 2, mean)
  total_var <- mean(apply(X, 1, euclidean_norm)^2)
  list(pc_var = pc_var, total_var = total_var)
}
