# Reference: Principal component analysis using QR decomposition (2012)
#==============================================================================
# Helper functions
#==============================================================================
#' Compute PCA rotation matrix using QR decomposition.
#' @param X n x p data matrix, n is the number of data, p is the number of features.
#' @param m desired number of principal components.
qr_pca <- function(X, m) {
  QR_fac <- qr(t(X))
  Q <- qr.Q(QR_fac)
  R <- qr.R(QR_fac)
  U <- svd(R, nu = m)$u
  Q %*% U
}
#' Compute PCA rotation matrix using SVD decomposition.
svd_pca <- function(X, m) {
  svd(t(X), nu = m, nv = 0)$u
}

#' Create a low-rank matrix.
#' @param n number of data
#' @param p number of features
#' @param r rank of the matrix
#' @examples 
#' X <- create_low_rank_matrix(100, 30, 10)
#' Matrix::rankMatrix(X)
create_low_rank_matrix <- function(n, p, r) {
  X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  if (r > min(n,p)) return(X)
  svd_X <- svd(X, nu = r, nv = r)
  svd_X$u %*% diag(svd_X$d[1:r]) %*% t(svd_X$v)
}

#==============================================================================
# Main
#==============================================================================
n <- 3000
p <- 10000
r <- 1000
m <- 10

# Check "qr_pca" gives the right solution as the direct "svd" approach.
# We take absolute values for the comparison because the factorisation is 
# unique up to a multiplicative constant of +/-1
X <- create_low_rank_matrix(n, p, r)
rotation_1 <- qr_pca(X, m = m)
rotation_2 <- svd_pca(X, m = m)
all.equal(abs(rotation_1), abs(rotation_2)) 

# Performance comparison
library(microbenchmark)
microbenchmark(
  r1 <- qr_pca(X, m), 
  r2 <- svd_pca(X, m = m),
  times = 100
)

#==============================================================================
# Results and Summary
# Unit: seconds
#         expr               min      lq    mean  median      uq     max  neval  cld
# r1 <- qr_pca(X, m)      829.34  840.30  848.40  850.79  855.51  861.90     10    a
# r2 <- svd_pca(X, m = m) 884.10  909.14  924.14  930.50  940.90  945.32     10    b
#------------------------------------------------------------------------------
#' Summary: The improvement of using the QR methods really comes from the fact 
#' that the data matrix X has a rank much smaller than the dimension (both n 
#' and p), so the QR step could "throw away" a big chunk of the matrix, easing 
#' the burdle of the subsequent SVD step. The result doesn't hold when the rank
#' of the matrix is close to the matrix dimension.
#' Opinion: The result seems somewhat artifical. Basically, it's like 
#' preconditioning the matrix before attempting the SVD, but I imagine there 
#' are some situations in practice that could benefit from this. 
#' Take-away: PCA factorisation can be speeded up if the matrix has low rank.
