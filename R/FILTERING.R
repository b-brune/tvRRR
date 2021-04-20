#' Filtering:
#'
#' Run the Kalman filter for a fixed set of parameters and prespecified model
#'
#' @param X the predictors
#' @param y the target variable
#' @param u additional predictors that come in as full rank
#' @param model specifies whether we're filtering in model "A" or model "B"
#' @param beta starting value for the algorithm, or time-constant parameter matrix (q x d)
#' @param alpha starting value for the algorithm, or the time-constant parameter matrix (p x d)
#' @param Gamma (optional), the time constant full rank (t x k) matrix
#' @param P_00 starting covariance for the algorithm (p*d x p*d matrix for A, (q*d x q*d) matrix for B)
#' @param Sigma column covariance of the states alpha_t (symmetric d x d matrix)
#' @param Omega error covariance in the measurement equation (symmetric p x p matrix)
#' @param d latent dimension (min. 1, no default)
#'
#' @returns An object of class \code{tvRRR} (LINK TO DOCUMENTATION OF tvRRR CLASS ->
#' DAS IST EIGENTLICH GELOGEN, mal gucken wie wir das lösen können)
#'
#' @export


eval_tvRRR <- function(X, y, u = NULL, d, model = "A", alpha, beta, ...) {

  if (model == "A") kf <- filter_modelA(X = X, y = y, u = u, alpha_00 = alpha, beta = beta, d = d, ...)

  if (model == "B") kf <- filter_modelB(X = X, y = y, u = u, alpha = alpha, beta_00 = t(beta), d = d, ...)

  class(kf) <- "tvRRR"

  return(kf)
}


## ############################################################################
##
## Filtering functions for models A and B
##
## ############################################################################


#' Function that runs the Kalman filter
#'
#' Input:
#' @param y the target variable (t x p matrix)
#' @param X the predictors (t x q matrix)
#' @param u (optional) additional predictors that do not necessarily vary in
#'          time (t x k matrix)
#' @param beta the time-constant part of the matrix (either NULL or a starting value)
#' @param alpha_00 starting value for the algorithm (either NULL or a starting value)
#' @param gamma (optional), the time constant full rank (t x k) matrix
#' @param P_00 starting covariance for the algorithm (p*d x p*d matrix)
#' @param Sigma column covariance of the states alpha_t (symmetric d x d matrix)
#' @param Omega error covariance in the measurement equation (symmetric p x p matrix)
#' @param d latent dimension (min. 1, no default)
#'
#' Output:
#' A named list of lists with elements
#' - states: filtered (the filtered states) -- one state matrix per row (t + 1 x p * d)
#'           smoothed (the smoothed states) -- one state matrix per row (t + 1 x p * d)
#' - covariances: the filtered and smoothed covariances and lag-1 covariances
#'     `P_t^t` filtered covariances -- array of dimensions (t+1, p*d, p*d)
#'     `P_t^t-1`predicted covariances -- (t, p*d, p*d)
#'     `P_t^T` smoothed covariances -- (t+1, p*d, p*d)
#'     `P_t-1t-2^T` smoothed lag-1 covariances -- (t, p*d, p*d)
#' - data: the data handed over to the algorithms
#'     `X` predictors -- (t, q)
#'     `y` responses -- (t, p)
#'     `u` additional predictors -- (t, k)
#'     `Z` transition matrices (X_t'beta (x) I_p) -- (t, p, p*d)
#' - parameters used during filtering: Sigma (d, d), Omega (p, p), beta (q, d)
#'
#' Symmetry of the matrices is enforced using the hint in
#' https://math.stackexchange.com/questions/2335831/covariance-matrix-p-for-an-extended-kalman-filter-not-symmetric
#' @export

filter_modelA <- function(y, X, u = NULL,
                          beta, alpha_00, Gamma = NULL,
                          P_00, Sigma, Omega, d,
                          return_covariances = TRUE,
                          Gamma_rrr = "identity"){

  t <- nrow(y)
  p <- ncol(y)

  ## -> directly work with vectorized alphas, i.e. store them in a matrix
  ## initialize matrices for storage of state estimated
  alpha_pred <- matrix(NA, t, p * d,     dimnames = list(paste0("t=", 1:t)))
  alpha_upd  <- matrix(NA, t + 1, p * d, dimnames = list(paste0("t=", 0:t)))
  alpha_s    <- matrix(NA, t + 1, p * d, dimnames = list(paste0("t=", 0:t)))

  ## initialize arrays for covariance storage and kalman gain and J-matrix
  P_tt  <- array(NA, c(t + 1, p * d, p * d), dimnames = list(paste0("t=", 0:t), NULL, NULL))
  P_tt1 <- array(NA, c(t, p * d, p * d),     dimnames = list(paste0("t=", 1:t), NULL, NULL))
  K_t   <- array(NA, c(t, p*d, p),           dimnames = list(paste0("t=", 1:t), NULL, NULL))
  P_tT  <- array(NA, c(t + 1, p * d, p * d), dimnames = list(paste0("t=", 0:t), NULL, NULL))
  J_t   <- array(NA, c(t, p * d, p * d),     dimnames = list(paste0("t=", 0:(t-1)), NULL, NULL))
  P_cov <- array(NA, c(t, p * d, p * d),     dimnames = list(paste0("t=", 1:t, ",t-1=", 0:(t-1)), NULL, NULL))


  ## initial conditions
  alpha_upd[1, ] <- alpha_00
  P_tt[1, , ]    <- P_00

  ## Make transition matrices from beta and x
  Z <- sapply(1:t, function(i) kronecker(crossprod(X[i, ], beta), diag(p)))
  dim(Z) <- c(p, p*d, t); Z <- aperm(Z, c(3, 1, 2)); dimnames(Z)[[1]] <- paste0("t=", 1:t)

  ## FILTERING

  state_cov <- kronecker(Sigma, diag(p))

  for (i in 1:t) {

    Zi <- Z[i, , ]

    # Prediction step
    alpha_pred[i, ] <- alpha_upd[i, ]
    P_tt1[i, , ]    <- P_tt[i, , ] + state_cov

    P_tt1i <- P_tt1[i, , ]

    # Calculate the Kalman gain
    K_t[i, , ] <- P_tt1i %*%
      crossprod(Zi, matpow(tcrossprod(Zi %*% P_tt1i, Zi) + Omega, -1, symmetric = TRUE))

    K_ti <- K_t[i, , ]

    # Updating step
    alpha_upd[i + 1, ] <- alpha_pred[i, ] +
      K_ti %*% (y[i, ] - Zi %*% alpha_pred[i, ] -
                  if (is.null(u)) { 0 } else { Gamma %*% u[i, ] })

    P_tt[i + 1, , ]    <- (diag(p * d) - K_ti %*% Zi) %*% P_tt1i
  }

  ## SMOOTHING
  alpha_s[t + 1, ] <- alpha_upd[t + 1, ]
  P_tT[t + 1, , ]  <- P_tt[t + 1, , ]

  for (i in t:1) {

    P_tti <- P_tt[i, , ]
    P_tt1i <- P_tt1[i, , ]

    J_t[i, , ] <- P_tti %*% matpow(P_tt1i, -1, symmetric = TRUE)

    J_ti <- J_t[i, , ]

    alpha_s[i, ] <- alpha_upd[i, ] + J_ti %*% (alpha_s[i + 1, ] - alpha_pred[i, ])

    P_tT[i, , ] <- P_tti + J_ti %*% (P_tT[i + 1, , ] - P_tt1i) %*% t(J_ti)

  }

  ## COVARIANCE SMOOTHING
  P_cov[t, , ] <- (diag(p * d) - K_t[t, , ] %*% Z[t, , ]) %*% P_tt[t, , ]

  for (i in t:2) {
    P_cov[i - 1, , ] <- P_tt[i, , ] %*% t(J_t[i - 1, , ]) +
      J_t[i, , ] %*% (P_cov[i, , ] - P_tt[i, , ]) %*% t(J_t[i - 1, , ])
  }


  return(list(states = list(
    onestepahead = alpha_pred,
    filtered = alpha_upd,
    smoothed = alpha_s),
    covariances = if (return_covariances) {
      list(`P_t^t` = P_tt,
           `P_t^t-1` = P_tt1,
           `P_t^T` = P_tT,
           `P_t-1t-2^T` = P_cov)
    } else NULL,
    prediction_covariance = P_tT[t + 1 , , ],
    data = list(X = X, y = y, u = u, Z = Z),
    parameters = list(Sigma = Sigma[, , drop = FALSE],
                      Omega = Omega,
                      beta = beta,
                      Gamma = Gamma)))
}




#' Function that runs the Kalman filter for model (B)
#'
#' It processes beta', not beta!!
#' @export
filter_modelB <- function(X, y, u = NULL,
                          P_00, Sigma, Omega,
                          beta_00, alpha, Gamma = NULL,
                          d,
                          return_covariances = TRUE,
                          Gamma_rrr = "identity") {

  t <- nrow(y)
  p <- ncol(y)
  q <- ncol(X)

  # INITIALIZATION OF STORAGE
  # -------
  beta_pred <- matrix(NA, t, q * d,   dimnames = list(paste0("t=", 1:t)))
  beta_upd  <- matrix(NA, t+1, q * d, dimnames = list(paste0("t=", 0:t)))
  beta_s    <- matrix(NA, t+1, q * d, dimnames = list(paste0("t=", 0:t)))

  ## initialize arrays for covariance storage and kalman gain and J-matrix
  P_tt  <- array(NA, c(t + 1, q * d, q * d), dimnames = list(paste0("t=", 0:t), NULL, NULL))
  P_tt1 <- array(NA, c(t, q * d, q * d),     dimnames = list(paste0("t=", 1:t), NULL, NULL))
  K_t   <- array(NA, c(t, q * d, p),         dimnames = list(paste0("t=", 1:t), NULL, NULL))
  P_tT  <- array(NA, c(t + 1, q * d, q * d), dimnames = list(paste0("t=", 0:t), NULL, NULL))
  J_t   <- array(NA, c(t, q * d, q * d),     dimnames = list(paste0("t=", 0:(t-1)), NULL, NULL))
  P_cov <- array(NA, c(t, q * d, q * d),     dimnames = list(paste0("t=", 1:t, ",t-1=", 0:(t-1)), NULL, NULL))
  # ----

  # initial conditions
  beta_upd[1, ] <- beta_00
  P_tt[1, , ]   <- P_00

  # make transition matrices
  Z <- sapply(1:t, function(i) kronecker(t(X[i, ]), alpha))
  dim(Z) <- c(p, q*d, t); Z <- aperm(Z, c(3, 1, 2)); dimnames(Z)[[1]] <- paste0("t=", 1:t)

  state_cov <- kronecker(diag(q), Sigma)
  # FILTERING:

  for (i in 1:t) {

    Zi <- Z[i, , ]

    # Prediction step
    beta_pred[i, ] <- beta_upd[i, ]
    P_tt1[i, , ] <- P_tt[i, , ] + state_cov

    P_tt1i <- P_tt1[i, , ]

    # Kalman gain
    K_t[i, , ] <- P_tt1i %*%
      crossprod(Zi, matpow(tcrossprod(Zi %*% P_tt1i, Zi) + Omega, -1, symmetric = TRUE))
    K_ti <- K_t[i, , ]

    # Updating step
    beta_upd[i + 1, ] <- beta_pred[i, ] +
      K_ti %*% (y[i, ] - Zi %*% beta_pred[i, ] -
                        if (is.null(u)) { 0 } else { Gamma %*% u[i, ] })

    P_tt[i + 1, , ] <- (diag(q * d) - K_ti %*% Zi) %*% P_tt1i
  }


  # SMOOTHING:

  beta_s[t + 1, ] <- beta_upd[t + 1, ]
  P_tT[t + 1, , ] <- P_tt[t + 1, , ]

  for (i in t:1) {
    P_tti <- P_tt[i, , ]
    P_tt1i <- P_tt1[i, , ]

    J_t[i, , ] <- P_tti %*% matpow(P_tt1i, -1, symmetric = TRUE)
    beta_s[i, ] <- beta_upd[i, ] + J_t[i, , ] %*% (beta_s[i + 1, ] - beta_pred[i, ])

    P_tT[i, , ] <- P_tti + J_t[i, , ] %*% (P_tT[i + 1, , ] - P_tt1i) %*% t(J_t[i, , ])
  }

  # COVARIANCE SMOOTHING:
  P_cov[t, , ] <- (diag(q * d) - K_t[t, , ] %*% Z[t, , ]) %*% P_tt[t, , ]

  for (i in t:2) {
    P_cov[i - 1, , ] <- P_tt[i, , ] %*% t(J_t[i - 1, , ]) +
      J_t[i, , ] %*% (P_cov[i, , ] - P_tt[i, , ]) %*% t(J_t[i - 1, , ])
  }

  return(list(states = list(
    onestepahead = beta_pred,
    filtered = beta_upd,
    smoothed = beta_s),
    covariances = if (return_covariances) {
      list(`P_t^t` = P_tt,
           `P_t^t-1` = P_tt1,
           `P_t^T` = P_tT,
           `P_t-1t-2^T` = P_cov)
    } else NULL,
    prediction_covariance = P_tT[t + 1 , , ],
    data = list(X = X, y = y, u = u, Z = Z),
    parameters = list(Sigma = Sigma[, , drop = FALSE],
                      Omega = Omega,
                      alpha = matrix(alpha, p, d),
                      Gamma = Gamma)))
}
