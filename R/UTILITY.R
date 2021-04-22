##
## Matrix processing utility functions
##

#' GENERALIZED MATRIX POWER FUNCTION FOR SYMMETRIC MATRICES (by Daniel Kapla)
#' @keywords internal

matpow <- function(A, pow, symmetric = FALSE, tol = 1e-7) {
  if (nrow(A) != ncol(A)) {
    stop("Expected a square matix, but 'A' is ", nrow(A), " by ", ncol(A))
  }
  # Case study for negative, zero or positive power.
  if (pow > 0) {
    if (pow == 1) { return(A) }
    # Perform SVD and return power as A^pow = U diag(d^pow) V'.
    svdA <- La.svd(A)
    return(svdA$u %*% ((svdA$d^pow) * svdA$vt))
  } else if (pow == 0) {
    return(diag(nrow(A)))
  } else {
    # make QR decomposition.
    qrA <- qr(A, tol = tol)
    # Check rank of A.
    if (qrA$rank == nrow(A)) {
      # Full rank, calc inverse the classic way using A's QR decomposition
      return(matpow(solve(qrA), abs(pow), tol = tol))
    } else {
      # For singular matrices use the SVD decomposition for the power
      svdA <- svd(A)
      # Get (numerically) positive singular values.
      positives <- svdA$d > tol * svdA$d[1]
      # Apply the negative power to positive singular values and augment
      # the rest with zero.
      d <- c(svdA$d[positives]^pow, rep(0, sum(!positives)))
      # The pseudo inverse as A^pow = V diag(d^pow) U' for pow < 0.
      if (symmetric) {
        return(svdA$v %*% (d * t(svdA$v)))
      } else  {
        return(svdA$v %*% (d * t(svdA$u)))
      }
    }
  }
}

#' Approximates kronecker product decomposition. (Daniel Kapla)
#'
#' Approximates the matrices `A` and `B` such that
#'      C = A %x% B
#' with `%x%` the kronecker product of the matrixes `A` and `B`
#' of dimensions `dimA` and `dimB` respectively.
#'
#' @param C desired kronecker product result.
#' @param dimA length 2 vector of dimensions of \code{A}.
#' @param dimB length 2 vector of dimensions of \code{B}.
#'
#' @return list with attributes `A` and `B`.
#'
#' @examples
#' A <- matrix(seq(14), 7, 2)
#' B <- matrix(c(TRUE, FALSE), 3, 4)
#' C <- kronecker(A, B) # the same as 'C <- A %x% B'
#' approx.kronecker(C, dim(A), dim(B))
#'
#' @seealso C.F. Van Loan / Journal of Computational and Applied Mathematics
#'          123 (2000) 85-100 (pp. 93-95)
#'
#' @keywords internal
approx.kronecker <- function(C, dimA, dimB) {

  dim(C) <- c(dimB[1L], dimA[1L], dimB[2L], dimA[2L])
  R <- aperm(C, c(2L, 4L, 1L, 3L))
  dim(R) <- c(prod(dimA), prod(dimB))

  svdR <- svd(R, 1L, 1L)

  return(list(
    A = array(sqrt(svdR$d[1]) * svdR$u, dimA),
    B = array(sqrt(svdR$d[1]) * svdR$v, dimB)
  ))
}

#' Check whether a matrix is symmetric
#' @keywords internal
is.symmetric <- function(A) {
  stopifnot(is.matrix(A),
            dim(A)[1] == dim(A)[2])

  all.equal(A, t(A))
}

# -----------------------------------------------------------------------------

##
## Functions for random orthogonal matrices
##

#' POLAR PROJECTION OF A MATRIX ONTO ITS ORTHOGONAL PART
#' @keywords internal
polar_proj <- function(X) return(X %*% matpow(crossprod(X), -0.5))

#' draw random matrix from Stiefel manifold
#' @keywords internal
rstiefel <- function(p, d) {
  if (d > p) stop("The matrix needs to have more rows than columns.")

  polar_proj(matrix(rnorm(p * d), p, d))
}

#' #' ORTHOGONAL PROJECTION ONTO SPAN(A)
#' @keywords internal
orth_proj <- function(A) {
  tcrossprod(A %*% matpow(crossprod(A), -1), A)
}


# -----------------------------------------------------------------------------

##
## Different subspace distances
##

#' NORMALIZED SUBSPACE DISTANCE (FROBENIUS OF PROJECTION DIFFERENCE)
#' @keywords internal
subsp_dist <- function(X, Y) {
  #  if (all.equal(dim(X) != dim(Y)) stop("X and Y need to have the same dimensions.")
  if (is.vector(X)) X <- matrix(X, length(X), 1)
  if (is.vector(Y)) Y <- matrix(Y, length(Y), 1)

  return(norm(orth_proj(X) - orth_proj(Y), 'F') / sqrt(2 * ncol(X)))
}

#' UNNORMALIZED GRASSMANN DISTANCE (BASED ON PRINCIPAL ANGLES)
#' @keywords internal
grassmann_dist <- function(A, B, tol = max(dim(A))* sqrt(.Machine$double.eps),
                           normalized = FALSE) {

  if (is.vector(A)) A <- as.matrix(A)
  if (is.vector(B)) B <- as.matrix(B)

  if (nrow(A) != nrow(B)) stop("The inputs A and B are not compatible")

  # Create orthogonal bases from A and B
  orthA <- qr.Q(qr(A, tol))
  orthB <- qr.Q(qr(B, tol))

  # calculate the angles
  r <- min(ncol(A), ncol(B))
  angles <- svd(crossprod(orthA, orthB), nu=0, nv=0)$d

  # Return the grassmann distance
  return(ifelse(normalized, (sum(acos(pmin(1, angles))^2)^0.5) / sqrt(r * (pi/2)^2),
                sum(acos(pmin(1, angles))^2)^0.5))
}


# -----------------------------------------------------------------------------

##
## Time-constant reduced rank regression
##

#' RRR with additional full-rank predictors
#' taken from Reinsel and Velu (1998), Ch. 3, Thm. 3.1
#' u corresponds to their z
#'
#' Model: y = mu + ABX + Du + error
#' @export

rr.regression <- function(X, y, u = NULL, rank, Gamma_type = "identity") {

  if (nrow(X) != nrow(y)) {
    stop("X, y (and u) need to have the same dimensions.")
  }

  if (!is.null(u)) {
    if (nrow(X) != nrow(u)) stop("X, y and u need to have the same dimensions.")
  }

  if (is.null(u)) {

    Gamma_sq <- switch(
      Gamma_type,
      identity = diag(ncol(y)),
      cov_y = matpow(cov(y), -1/2),
      OLS = matpow(cov(y - X %*% MASS::ginv(X) %*% y), -1/2)
    )

    Sigma_xx <- cov(X)
    Sigma_yy <- cov(y)
    Sigma_yx <- cov(y, X)

    V <- svd(Gamma_sq %*% tcrossprod(Sigma_yx %*% matpow(Sigma_xx, -1),
               Sigma_yx) %*% Gamma_sq)$u[, 1:rank, drop = F]

    A <- matpow(Gamma_sq, -1) %*% V
    B <- t(V) %*% Gamma_sq %*% Sigma_yx %*% matpow(Sigma_xx, -1)
    mu <- colMeans(y) - A %*% B %*% colMeans(X)
    C <- A %*% B
    D <- NULL

    return(list(A = A, B = B, C = C, D = D, mu = mu))
  }

  Sigma_xx <- cov(X)
  Sigma_yy <- cov(y)
  Sigma_uu <- cov(u)
  Sigma_yx <- cov(y, X)
  Sigma_xu <- cov(X, u)
  Sigma_yu <- cov(y, u)

  Sigma_yx.u <- Sigma_yx - tcrossprod(Sigma_yu  %*% matpow(Sigma_uu, -1), Sigma_xu)
  Sigma_xx.u <- Sigma_xx - tcrossprod(Sigma_xu %*% matpow(Sigma_uu, -1), Sigma_xu)

  Gamma_sq <- switch(
    Gamma_type,
    identity = diag(ncol(y)),
    cov_y = matpow(Sigma_yy - tcrossprod(Sigma_yu %*% Sigma_uu, Sigma_yu), -1),
    OLS = matpow(Sigma_yy - cbind(Sigma_yx, Sigma_yu) %*%
      matpow(rbind(cbind(Sigma_xx, Sigma_xu),
                   cbind(t(Sigma_xu), Sigma_uu)), -1) %*% t(cbind(Sigma_yx, Sigma_yu)), -1)
  )

  V <- svd(
    Gamma_sq %*% tcrossprod(Sigma_yx.u %*% matpow(Sigma_xx.u, -1), Sigma_yx.u) %*% Gamma_sq, nv = 0
  )$u[, 1:rank, drop = F]

  A <- matpow(Gamma_sq, -1)%*% V
  B <- crossprod(V, Gamma_sq %*% Sigma_yx.u %*% matpow(Sigma_xx.u, -1))
  C <- A %*% B

  D <- Sigma_yu %*% matpow(Sigma_uu, -1) - C %*% Sigma_xu %*% matpow(Sigma_uu, -1)

  mu <- colMeans(y) - C %*% colMeans(X) - D %*% colMeans(u)


  list(A = A, B = B, C = C, D = D, mu = mu)
}



