# This file requires "util.R"
#==============================================================================
# Reference
# (2007) Fast principal component analysis using fixed-point algorithm
# http://www.sciencedirect.com/science/article/pii/S0167865507000438
#==============================================================================
# Fast PCA algorithm using iterative method
#' @param X data matrix, n data points x p features
#' @param h number of principal components desired
fastPCA <- function(X, h, tol = 1e-10) {
  XTX <- t(X) %*% X
  p <- ncol(X)
  Q <- matrix(0, nrow = ncol(X), ncol = h)
  for (i in 1:h) {
    eig_vec <- rnorm(p)
    err <- Inf
    # count <- 0   # verifying the claim of "few iterations" in the paper. TRUE.
    while (err > tol) {
      new_eig_vec <- drop(XTX %*% eig_vec)
      if (i > 1) {
        # do Gram-Schmidt step
        for (j in 1:(i-1)) {
          new_eig_vec <- new_eig_vec - a_project_onto_b(new_eig_vec, Q[,j]) * Q[,j]
        }
      }
      new_eig_vec <- new_eig_vec / euclidean_norm(new_eig_vec)
      err <- abs(inner_product(eig_vec, new_eig_vec) - 1)
      eig_vec <- new_eig_vec
      # count <- count + 1
    }
    # print(count)
    Q[,i] <- eig_vec
  }
  Q
}

#' Compute PCA rotation matrix using SVD decomposition.
#' [For comparison purpose]
svd_pca <- function(X, m) {svd(t(X), nu = m, nv = 0)$u}

#==============================================================================
# Main
#==============================================================================
source("util.R")
n <- 10000 # The paper suggested 100, but it doesn't show improvement. See 
# previous commit for details
p <- 2000 # 100 to 1000, then 2000 to 4000
h <- 10
X <- matrix(rnorm(n*p), n, p)
library(microbenchmark)
microbenchmark(
  ans <- fastPCA(X, h, tol = 0.01),
  ans2 <- svd_pca(X, h),
  times = 10
)

#==============================================================================
# Results and Summary
# Unit: seconds
#                             expr    min     lq   mean median     uq    max neval cld
# ans <- fastPCA(X, h, tol = 0.01) 108.20 108.80 110.55 110.62 111.12 113.92    10  a 
# ans2 <- svd_pca(X, h)            384.04 395.15 397.01 398.45 402.68 403.52    10   b
#------------------------------------------------------------------------------
# Sidenote: When tolerance is set to 0.01 as suggested in the paper, the eigenvectors 
# differ a lot from the 'true' ones. They still form a orthonomal basis, so the MSE 
# computed in the paper is legitimate, but they are off in terms of the explained 
# variance and the ordering of the eigenvectors. This goes away when the tolerance 
# is set to around 1e-10, but it leads to a significant increase of computation time.
