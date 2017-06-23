#' Gram Schmidt
#' @examples
#' check: 
#' A <- rbind(c(2,-2,18), c(2,1,0), c(1,2,0))
#' gram_schmidt(A)
#' list(Q = qr.Q(qr(A)), R = qr.R(qr(A)))
gram_schmidt <- function(X) {
  X <- as.matrix(X)
  n <- ncol(X)
  p <- nrow(X)
  
  Q <- matrix(0, n, p)
  R <- matrix(0, p, p)
  
  x1 <- X[,1]
  Q[,1] <- x1 / euclidean_norm(x1)
  R[1,1] <- a_project_onto_b(x1, Q[,1])
  for (k in 2:ncol(X)) {
    xk <- X[,k]
    qk <- xk
    for (i in 1:(k-1)) {
      # modified
      R[i, k] <- a_project_onto_b(qk, Q[,i])
      qk <- qk - R[i, k] * Q[,i]
    }
    Q[,k] <- qk / euclidean_norm(qk)
    R[k,k] <- a_project_onto_b(xk, Q[,k])
  }
  list(Q = Q, R = R)
}

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
