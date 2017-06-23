#==============================================================================
# References:
# 1. (2008) Principal Component Analysis based on L1-norm Maximization
#    http://ajou.ac.kr/~nojunk/papers/L1PCA_TPAMI.pdf
# 2. (2013) A Fast Implementation of PCA-L1 Using Gram-Schmidt Orthogonalization
#    https://www.jstage.jst.go.jp/article/transinf/E96.D/3/E96.D_559/_article
#==============================================================================
# Find the principal components under the L1-norm.
#' @param X data matrix
#' @param m number of principal components desired
#' @param w initial vector
#' @examples 
#' X <- matrix(rnorm(1000), 100, 10)
#' rotation <- PCA_L1_max(X, 3)
#' Z <- X %*% t(rotation)
#' # Toy example from the paper
#' X <- cbind(c(-6:-2, 10, 0:4), c(-5:5))
#' w <- PCA_L1_max(X, 1)  # expect w_L1 = [0.8, 0.6]
PCA_L1_max <- function(X, m, w, turbo = F, remove_mean = T, tol = 1e-6) {
  p <- ncol(X)
  if (remove_mean) X <- center(X)
  if (missing(w)) w <- rnorm(p)
  if (turbo) basis <- prcomp(X)$rotation
  if (m > p) m <- p
  W <- matrix(0, m, p)
  W[1,] <- PCA_L1_max_single_vec(X, w, tol)
  if (m == 1) return(W)
  for (j in 2:m) {
    for (i in 1:nrow(X)) {
      xi <- X[i,]
      wj_1 <- W[j-1, ]
      X[i, ] <- xi - wj_1 * inner_product(wj_1, xi)
    }
    w <- rnorm(p)
    if (turbo) w <- turbo_start(basis[,j], W, j)
    W[j,] <- PCA_L1_max_single_vec(X, w, tol)
  }
  W
}

PCA_L1_max_single_vec <- function(X, w, tol = 1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  tX <- t(X)
  if (missing(w)) w <- rnorm(p)
  w <- w / euclidean_norm(w)
  
  p_i <- numeric(n)
  converge <- FALSE
  while (!converge) {
    old_w <- w
    polarity <- as.numeric(X %*% w)
    p_i <- 1 - as.numeric(polarity < 0) * 2
    w <- apply(tX * rep(p_i, each = ncol(X)), 1, sum)
    w <- w / euclidean_norm(w)
    # convergence check
    if (sum(abs(old_w - w)) < tol) {
      if (any(polarity == 0)) {
        w <- w + rnorm(p)
        w <- w / euclidean_norm(w)
      } else {
        converge <- T
      }
    }
  }
  w
}

get_pc_variance <- function(W, X) {
  pc_var <- apply((X %*% t(W))^2, 2, mean)
  total_var <- mean(apply(X, 1, euclidean_norm)^2)
  list(pc_var = pc_var, total_var = total_var)
}

#==============================================================================
# A Fast implementation of PCA-L1 using Gram-Schmidt orthogonalisation
# See reference at the top
#==============================================================================
turbo_start <- function(w, W, j) {
  for (k in 1:(j-1)) {
    w <- w - a_project_onto_b(w, W[k,]) * W[k,]
  } 
  w
}

#==============================================================================
# Utility functions
#==============================================================================
inner_product <- function(v0, v1) {
  sum(v0 * v1)
}

euclidean_norm <- function(v0) {
  sqrt(inner_product(v0, v0))
}

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
