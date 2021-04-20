## ############################################################################
##
## HELP FUNCTIONS FOR PARAMETER UPDATES
##
## ############################################################################


#' HELP 1A: update beta during the EM algorithm
#' @keywords internal

X_update_beta <- function(d, Omega, Gamma, X, y, u, kf, p, q, t) {

  # Update beta by calculating the roots of the first derivative as calculated
  # in the pdf, where our previous beta serves as a starting value

  s <- kf$states$smoothed

  Omega_inv <- matpow(Omega, -1, symmetric = TRUE)

  # Calculate the matrices \bar{P}_t
  P_bar <- apply(kf$covariances$`P_t^T`[-1, , ], 1, function(P) {
    P_t_bar <- matrix(NA, d, d)

    # symmetric, so we only need the upper triangle
    for (i in 1:d) {
      for (j in i:d) {
        P_t_bar[i, j] <- P_t_bar[j, i] <- {
          sum(diag((Omega_inv %*% P[((i - 1) * p + 1):(i * p), ((j - 1) * p + 1):(j * p)])))
        }
      }
    }

    P_t_bar
  })

  if (d == 1) P_bar <- matrix(P_bar, nrow = 1)

  pt1 <- matrix(
    rowSums(
      sapply(1:t, function(i) {
        kronecker(tcrossprod(X[i, ]), crossprod(matrix(s[i + 1, ], p, d), Omega_inv %*% matrix(s[i + 1, ], p, d))) +
          kronecker(tcrossprod(X[i, ]), matrix(P_bar[, i], d, d))
      })),
    q * d, q * d)

  y <- y - { if(!is.null(u)) t(Gamma %*% t(u)) else 0 }

  pt2 <- rowSums(sapply(1:t, function(i) {
    t(matrix(s[i + 1, ], p, d)) %*% Omega_inv %*% y[i, ] %*% t(X[i, ])
  }))

  vec_beta_transpose <- matpow(pt1, -1) %*% pt2

  return(matrix(vec_beta_transpose, q, d, byrow = TRUE))
}



#' HELP 1B: update alpha during the EM algorithm
#' @keywords internal

X_update_alpha <- function(d, Omega, Gamma = NULL,
                           X, y, u = NULL,
                           s, kf, p, q, t) {

  # Omega_inv <- matpow(Omega, -1)
  if (!is.null(u)) y <- y - t(Gamma %*% t(u))

  pt1 <- matrix(rowSums(sapply(1:t, function(i) y[i, ] %*% t(X[i, ]) %*%
                                 matrix(s[i + 1, ], q, d, byrow = T))), p, d)
  if (d == 1) {
    pt2 <- matrix(sum(sapply(1:t, function(i) {
      tcrossprod(matrix(s[i + 1, ], d, q) %*% X[i, ]) +
        kronecker(t(X[i, ]), diag(d)) %*% kf$covariances$`P_t^T`[i + 1, , ] %*%
        kronecker(X[i, ], diag(d))
    })), d, d)
  } else if (d > 1) {
    pt2 <- matrix(rowSums(sapply(1:t, function(i) {
      tcrossprod(matrix(s[i + 1, ], d, q) %*% X[i, ]) +
        kronecker(t(X[i, ]), diag(d)) %*% kf$covariances$`P_t^T`[i + 1, , ] %*%
        kronecker(X[i, ], diag(d))
    })), d, d)
  }

  return(pt1 %*% matpow(pt2, -1))
}


#' HELP 2: update Omega during the EM algorithm, optionally including the external variables
#' @keywords internal

X_update_Omega <- function(Z, Gamma, Omega_diagonal, y, u, kf, p, t) {

  s <- kf$states$smoothed

  if (Omega_diagonal) {
    # Update Omega
    C <- diag(matrix(rowSums(sapply(1:t, function(i) tcrossprod(y[i, ] - Z[i, , ] %*% s[i + 1, ] -
                                                   if(!is.null(u)) { Gamma %*% u[i, ] } else { 0 } ) +
                       Z[i, , ] %*% kf$covariances$`P_t^T`[i + 1, , ] %*% t(Z[i, , ]))), p, p))

    Omega_n <- diag(C / t)

  } else if (!Omega_diagonal) {

    Omega_n <- 1/t * matrix(rowSums(
      sapply(1:t, function(i) tcrossprod(y[i, ] - Z[i, , ] %*% s[i + 1, ] -
                                           if(!is.null(u)) { Gamma %*% u[i, ] } else { 0 } ) +
        Z[i, , ] %*% kf$covariances$`P_t^T`[i + 1, , ] %*% t(Z[i, , ]))), p, p)
  }

  if (det(Omega_n) <= 0) {
    warning("Omega is not positive definite and has been modified. \n")
    while (det(Omega_n <= 0)) {
      Omega_n <- Omega_n + min(1/(2*p), mean(abs(diag(Omega_n)))) * diag(p)
    }
  }


  if (!isTRUE(all.equal(Omega_n, t(Omega_n)))) {
    warning("Omega is not symmetric and has been symmetrized! \n")
    Omega_n <- 1/2 * (Omega_n + t(Omega_n))
  }

  Omega_n
}


#' HELP 3: Calculate the current value of the Likelihood and (if required) the
#'         matrix C for parameter estimation
#' @keywords internal
X_calc_Q <- function(kf, p, d, q, t, return_C = FALSE, model = "A") {

  s <- kf$states$smoothed

  Z <- kf$data$Z
  y <- kf$data$y
  X <- kf$data$X
  u <- kf$data$u

  # Calculate the matrices C and E needed for the M-step and evaluation of
  ## the Likelihood
  C <- rowSums(sapply(2:(t+1), function(i) {
    tcrossprod(s[i, ]) + kf$covariances$`P_t^T`[i, ,] +
      tcrossprod(s[i - 1, ]) + kf$covariances$`P_t^T`[i - 1, ,] -
      (tcrossprod(s[i, ], s[i-1, ]) + kf$covariances$`P_t-1t-2^T`[i - 1, , ]) -
      t(tcrossprod(s[i, ], s[i-1, ]) + kf$covariances$`P_t-1t-2^T`[i - 1, , ])
  }))


  if (model == "A") {
    C <- matrix(C, p * d, p * d)
    pref <- t * p * log(det(kf$parameters$Sigma))
    kron_mat <- kronecker(matpow(kf$parameters$Sigma, -1, symmetric = TRUE), diag(p))
  } else if (model == "B") {
    C <- matrix(C, q * d, q * d)
    pref <- t * q * log(det(kf$parameters$Sigma))
    kron_mat <- kronecker(diag(q), matpow(kf$parameters$Sigma, -1, symmetric = TRUE))
  }

  if(!is.null(u)) y <- t(t(y) - kf$parameters$Gamma %*% t(u))

  E <- matrix(rowSums(
    sapply(2:(t+1), function(i) {
      tcrossprod(y[i - 1, ]) -
        y[i - 1, ]%*% t(s[i, ]) %*% t(Z[i - 1, , ]) -
        Z[i - 1, , ] %*% s[i, ] %*% y[i - 1, ] +
        Z[i - 1, , ] %*% (tcrossprod(s[i, ]) + kf$covariances$`P_t^T`[i, , ]) %*% t(Z[i - 1, , ])
    })), p, p)



  ## EVALUATE LIKELIHOOD
  Q <- log(det(kf$covariances$`P_t^T`[1, , ])) + ifelse(model == "A", p, q) * d +
    pref + sum(diag(kron_mat %*% C)) +
    t * log(det(kf$parameters$Omega)) + sum(diag(matpow(kf$parameters$Omega, -1) %*% E))

  if (return_C) return(list(Q = Q, C = C))

  Q
}

#' HEPL 3.1:
#' Calculate the true data likelihood based on the latent states
#' @keywords internal

X_eval_lik <- function(kf, Omega, beta = NULL, alpha = NULL, Gamma = NULL, d, y, X) {

  #s <- kf$states$smoothed
  s <- kf$states$onestepahead
  t <- nrow(X)
  p <- ncol(y)
  q <- ncol(X)
  u <- kf$data$u

  if (!is.null(u)) y <- t(t(y) - Gamma %*% t(u))


  if (!is.null(beta)) {
    return(
      sum(
        sapply(1:t, function(i) {
          Sigma_t <- kf$data$Z[i, , ] %*% kf$covariances$`P_t^t-1`[i, , ] %*% t(kf$data$Z[i, , ]) + Omega

          -0.5 * log(det(Sigma_t)) -
                 0.5 * crossprod(y[i, ] - matrix(s[i, ], p, d) %*% t(beta) %*% X[i, ],
                                 matpow(Sigma_t, -1, symmetric = TRUE) %*%
                                   (y[i, ] - matrix(s[i, ], p, d) %*% t(beta) %*% X[i, ]))
    })))

  } else if (!is.null(alpha)) {
    return(
      sum(
        sapply(1:t, function(i) {

          Sigma_t <- kf$data$Z[i, , ] %*% kf$covariances$`P_t^t-1`[i, , ] %*% t(kf$data$Z[i, , ]) + Omega

          - 0.5 * log(det(Sigma_t)) -
                 0.5 * crossprod(y[i, ] - alpha %*% matrix(s[i, ], d, q) %*% X[i, ],
                                 matpow(Sigma_t, -1, symmetric = TRUE) %*%
                                   (y[i, ] - alpha %*% matrix(s[i, ], d, q) %*% X[i, ]))
        }
        )))
  }
}


#' HELP 4: Update Sigma_c
#' @keywords internal

X_update_Sigma <- function(C, p, dim, q, t, model = "A") {

  Sigma <- NULL

  if (dim == 1) {
    if (model == "A") {
      Sigma <- 1 / (t * p) * matrix(sum(svd(C)$d), dim, dim)}
    else if (model == "B") {
      Sigma <- 1 / (t * q) * matrix(sum(svd(C)$d), dim, dim)}
  }
  if (dim > 1)  {
    if (model == "A") {
      with(svd(C), {
        Sigma <<- matrix(1 / (t*p) *
                           rowSums(sapply(seq_along(d),
                                          function(i) d[i] * crossprod(matrix(u[, i], p, dim)))), dim, dim)
      })
    } else if (model == "B") {
      with(svd(C), {
        Sigma <<- matrix(1 / (t*q) *
                           rowSums(sapply(seq_along(d),
                                          function(i) d[i] * tcrossprod(matrix(u[, i], dim, q)))), dim, dim)
      })
    }
  }

  Sigma
}

#' HELP 5: Update Gamma (optional)
#' @keywords internal

X_update_Gamma <- function(kf, Z) {

  s <- kf$states$smoothed
  k <- ncol(kf$data$u)
  p <- ncol(kf$data$y)
  t <- nrow(kf$data$y)

  pt1 <- matrix(rowSums(
    sapply(1:t, function(i) {
    tcrossprod(kf$data$y[i, ] - Z[i, , ] %*% s[i + 1, ], kf$data$u[i, ])
  })), p, k)

  pt2 <- matrix(rowSums(
    matrix(sapply(1:t, function(i) tcrossprod(kf$data$u[i, ])), k)), k, k)

  pt1 %*% matpow(pt2, -1)

}



## ############################################################################
##
## ACTUAL PARAMETER UPDATING
##
## ############################################################################

#' Perform one EM-step
#'
#' Function that updates the parameters from the current output of the KF-function
#' i.e. one updating step in the EM algorithm
#' @param kf a Kalman filter object, output from eval_tvRRR() / filter_modelA() / filter_modelB() function
#' @param p dimension of the target
#' @param d latent dimension
#' @param t length of the time series
#' @param q dimension of the predictors
#' @param tol tolerance for convergence in iterative parameter updating
#' @param Omega_diagonal logical, indicates whether Omega is assumed to be diagonal
#'
#' @return
#' List of updated parameters:
#'  \item{Omega}{}
#'  \item{beta}{(when fitting model A)}
#'  \item{alpha}{(when fitting model B)}
#'  \item{Sigma}{}
#'  \item{Gamma}{(when considering external variables)}
#'  \item{Q}{value of the expected likelihood Q(Theta|Theta^j)}
#'
#' @keywords internal

update_pars <- function(kf, p, d, t, q, # kf-Object and the dimensions of the model
                        tol = 1e-3, maxit = 50,
                        Omega_diagonal = FALSE,
                        model = "A"
                        # details on the estimation algorithms
) {

  ## Put the stuff we need often into parameters:

  s <- kf$states$smoothed

  Z <- kf$data$Z
  y <- kf$data$y
  X <- kf$data$X
  u <- kf$data$u

  # E-step: Evaluate the likelihood
  tmp <- X_calc_Q(kf = kf, p = p, d = d, t = t, q = q, return_C = TRUE, model = model)
  Q <- tmp$Q
  C <- tmp$C


  ## Update Sigma using the SVD based formulae we derived

  Sigma <- X_update_Sigma(C = C, p = p, dim = d, q = q, t = t, model = model)

  iter <- 1

  ## UPDATING FOR MODEL A

  if (model == "A") {

    beta_n  <- kf$parameters$beta # "starting value for the steps"
    Omega_n <- kf$parameters$Omega
    Gamma_n <- kf$parameters$Gamma

    repeat {
      beta <- beta_n; Omega <- Omega_n
      if (!is.null(Gamma)) Gamma <- Gamma_n

      # Start by updating beta given the current set of parameters

      # Calculate new transition matrices Z
      Z <- sapply(1:t, function(i) kronecker(crossprod(X[i, ], beta), diag(p)))
      dim(Z) <- c(p, p*d, t); Z <- aperm(Z, c(3, 1, 2)); dimnames(Z)[[1]] <- paste0("t=", 1:t)

      if (!is.null(Gamma)) Gamma_n <- X_update_Gamma(kf = kf, Z = Z)

      Omega_n <- X_update_Omega(Z = Z, Gamma = Gamma_n,
                                Omega_diagonal = Omega_diagonal,
                                y = y, u = u, kf = kf, p = p, t = t)
      beta_n <- X_update_beta(d = d, Omega = Omega_n, Gamma = Gamma_n,
                              X = X, y = y, u = u,
                              kf = kf, p = p, t = t, q = q)

      # Conditions for stopping
      if (is.null(Gamma)) {
        if (norm(Omega_n - Omega, 'm') < tol |
            subsp_dist(beta_n, beta) < tol |
            iter > maxit) break
      } else if (!is.null(Gamma)) {
        if (norm(Omega_n - Omega, 'm') < tol |
            subsp_dist(beta_n, beta) < tol |
            iter > maxit |
            norm(Gamma_n - Gamma, 'm') < tol) break
      }

      iter <- iter + 1
    }


  } else if (model == "B") {

    ## PARAMETER UPDATES FOR MODEL B

    alpha_n <- kf$parameters$alpha # "starting value for the algorithm"
    Omega_n <- kf$parameters$Omega

    Gamma_n <- kf$parameters$Gamma

    if (!is.null(u)) {
      repeat {
        alpha <- alpha_n; Gamma <- Gamma_n
        alpha_n <- X_update_alpha(d = d, Gamma = Gamma_n, Omega = kf$parameters$Omega,
                                  X = X,
                                  y = y, u = u, s = s, kf = kf, p = p, q = q, t = t)

        Z <- sapply(1:t, function(i) kronecker(t(X[i, ]), alpha_n))
        dim(Z) <- c(p, q * d, t); Z <- aperm(Z, c(3, 1, 2)); dimnames(Z)[[1]] <- paste0("t=", 1:t)

        Gamma_n <- X_update_Gamma(kf = kf, Z = Z)

        if (subsp_dist(alpha_n, alpha) < tol |
            norm(Gamma_n - Gamma, 'm') < tol) break
      }
    } else {
      alpha_n <- X_update_alpha(d = d, Omega = kf$parameters$Omega,
                                X = X,
                                y = y, s = s, kf = kf, p = p, q = q, t = t)
      Z <- sapply(1:t, function(i) kronecker(t(X[i, ]), alpha_n))
      dim(Z) <- c(p, q * d, t); Z <- aperm(Z, c(3, 1, 2)); dimnames(Z)[[1]] <- paste0("t=", 1:t)
    }

    Omega_n <- X_update_Omega(Z = Z, Gamma = Gamma, Omega_diagonal = Omega_diagonal, y = y, u = u,
                              kf = kf, p = p, t = t)

  }


  list(Sigma = Sigma[, , drop = FALSE],
       Q = Q, # Q from previous run of the filter... not the new Q
       Omega = Omega_n,
       beta  = if (model == "A") beta_n else NULL,
       alpha = if (model == "B") alpha_n else NULL,
       Gamma = Gamma_n)
}
