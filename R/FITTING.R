## ############################################################################
##
## tvRRR model selection based on BIC
##
## ############################################################################

#' Fit the tvRRR model
#'
#' The function fits the tvRRR model to an observed dataset using the EM-algorithm.
#' Rank selection is performed automatically based on the BIC criterion.
#' The model with the lowest BIC is selected.
#'
#' @param X predictors (t x q-dimensional numeric matrix)
#' @param y target (t x p-dimensional numeric matrix)
#' @param u additional predictors (t x k-dimensional numeric matrix)
#' @param model determines the model to be fit, either `"A"` or `"B"`.
#' @param d_max maximum latent dimension, needs to be an integer and can at most be min(p, q)
#' @param d (maximum) latent dimension, needs to be an integer and can at most be min(p, q)
#' @param select_rank logical, indicates whether the rank should be selected by BIC or whether
#'                    a model of rank `d` (or `d_max`) should be fit
#' @param Sigma_init optional, modifies the initial state covariance defined as \code{Sigma_init * diag(d)}
#' @param ...  additional parameters such as starting values handed over to the EM-algorithms
#'             for models A and B, see details.
#'
#' @details
#' The \code{tvRRR()} function fits the tvRRR models of class A and B.
#' The model
#' is a reduced rank regression model, where the coefficient matrix can be decomposed
#' in two ways: \eqn{C_t = \alpha_t \beta'} (Model A) or \eqn{C_t = \alpha \beta_t'} (Model B). For
#' details on the model fitting algorithm see Brune et al. (2021+).
#'
#' Precisely, model A is
#' \deqn{y_t = \alpha_t\beta'x_t + \Gamma u_t + \epsilon_t}
#'
#' and model B is given by
#'
#' \deqn{y_t = \alpha \beta'x_t + \Gamma u_t + \epsilon_t.}
#'
#' Both models are fitted by an EM-algorithm. The actual model fitting happens in
#' `fit_modelA` and `fit_modelB`. The function `tvRRR` calls these functions and,
#' if `select_rank = TRUE` performs model selection with BIC. Here, `d_max` denotes
#' the maximum dimension that should be fitted. If `d_max = NULL`, `d` is interpreted
#' as `d_max`. `tvRRR` determines starting values automatically, but also allows for
#' manually handing over appropriate starting values. For details, also see \code{\link{fit_modelA}()}and
#' \code{\link{fit_modelB}()}.
#'
#' Arguments that can be handed over to guide the behavior of the algorithm are
#'
#' \itemize{
#' \item `return_covariances` logical, indicates whether the state covariances
#'                          should be returned (might be necessary for evaluation
#'                          of the likelihood).
#' \item `silent` logical, indicates whether progress should be printed during
#'              model fitting
#' \item `maxit` maximum number of iterations for the EM algorithm
#' \item `tol_finish` tolerance for stopping the EM algorithm
#' \item `tol_EMstep` tolerance for iterative estimation during EM step
#' }
#'
#' Arguments regarding the starting values for the EM-algorithm are, if no initial
#' values are handed over:
#'
#' \itemize{
#' \item `initialize` either \code{"RRR"} or \code{"random"}, applies if no starting
#'                    values are handed over, then initialization is either carried
#'                    out randomly or from time-constant reduced rank regression.
#' \item `Gamma_rrr` if initialization is carried out by RRR, this specifies the
#'                   normalization of \eqn{\alpha} and \eqn{\beta}
#' }
#'
#'
#' If we want to hand over starting values explicitly:
#'
#' \itemize{
#' \item `alpha_00` initial value for the Kalman filter in model A
#' \item `alpha` starting value for the time-constant matrix \eqn{\alpha} in model B
#' \item `beta_00` initial value for the Kalman filter in model B. **ATTENTION:**
#'                 this has to be the transpose of \eqn{\beta}, i.e. \eqn{\beta'}
#' \item `beta` starting value for the time-constant matrix \eqn{\beta}
#' \item `Gamma` starting value for the fixed full-rank coefficient matrix
#' \item `P_00` starting state covariance (default `1000 * diag(p * d)` for
#'              model A and `1000 * diag(q * d)` for model B)
#' \item `Sigma` column covariance of states (default 0.1 * diag(d))
#' \item `Omega` error covariance (defaults to residual covariance from the starting values)
#' }
#'
#' Furthermore:
#'
#' `Omega_diagonal` logical, indicates whether Omega is assumed to be a
#'                  diagonal matrix (advisable if p is large)
#'
#' @return
#'  An object of class \code{tvRRR}, that is a named list of lists with elements
#' \describe{
#'  \item{states}{the filtered states, a named list with elements
#'  \itemize{\item filtered (the filtered states) -- one state matrix per row (t + 1 x p * d)
#'           \item smoothed (the smoothed states) -- one state matrix per row (t + 1 x p * d)
#'           \item one-step-ahead (one-step ahead predictions of the states -- one state matrix per row (t + 1 x p * d))
#'           }
#'           }
#'  \item{covariances}{the filtered and smoothed covariances and lag-1 covariances,
#'  if \code{return_covariances = TRUE}, a named list with elements
#'  \itemize{\item `P_t^t` filtered covariances -- array of dimensions (t+1, p*d, p*d)
#'     \item `P_t^t-1`predicted covariances -- (t, p*d, p*d)
#'     \item `P_t^T` smoothed covariances -- (t+1, p*d, p*d)
#'     \item `P_t-1t-2^T` smoothed lag-1 covariances -- (t, p*d, p*d),
#'     }
#'     else \code{NULL}}
#' \item{prediction_covariance}{contains `P_t^T[t+1, , ]` which is necessary for one-step
#'                              ahead prediction of \code{tvRRR} object.}
#' \item{data}{the data handed over to the algorithms, a named list with elements
#' \itemize{
#'     \item `X` predictors -- (t, q)
#'     \item `y` responses -- (t, p)
#'     \item `u` additional predictors -- (t, k)
#'     \item `Z` transition matrices (X_t'beta (x) I_p) -- (t, p, p*d)
#'     }
#'     }
#' \item{parameters}{the parameters used for filtering:
#' \itemize{
#'     \item Sigma -- (d, d)
#'     \item Omega -- (p, p)
#'     \item beta -- (q, d) (for model A)
#'     \item alpha -- (p, d) (for model B)
#' }}
#' \item{likelihoods}{list containing information on the data likelihood during
#' the fitting procedure, i.e.
#' \itemize{
#' \item \code{Q} the negative expected log likelihood obtained from EM
#' \item \code{logLik} actual data log Likelihood for each iteration
#' }
#' }
#' \item{convergence_information}{Message containing parameter stages at convergence}
#' \item{BIC}{Value of the information criterion for the selected model}
#' \item{rank}{The rank of the model}
#' }
#'
#' @seealso \code{\link[tvRRR]{fit_modelA}}, \code{\link[tvRRR]{fit_modelB}}
#'
#'
#'
#' @export

tvRRR <- function(X, y, u = NULL, model = "A",
                  d, select_rank = TRUE, d_max,
                  Sigma_init = 0.1, ...) {

  # p <- ncol(y)
  q <- ncol(X)
  # t <- nrow(X)

  if (!select_rank) {

    if (missing(d)) stop("Provide the latent dimension d or set `select_rank = TRUE`
                         for automatic determination.")

    fit <- fit_tvRRR(X = X, y = y, u = u, model = model, d = d, Sigma = Sigma_init * diag(d), ...)
    fit$rank <- d

    return(fit)
  }


  else if (select_rank) {


    if (missing(d_max) & !missing(d)) d_max <- d
    if (missing(d_max) & missing(d)) d_max <- q

    # if (length(criterion) > 1) criterion <- criterion[1]

    stopifnot(
      "d_max needs to be numeric" = { is.numeric(d_max) },
      "d_max needs to be an integer number" = { d_max %% 1 == 0 },
      "d_max must not be larger than min(ncol(X), ncol(y))" = { d_max <= min(ncol(X), ncol(y)) }
    )

    crit <- rep(NA, d_max)

    models <- vector("list", d_max)

    for (d_try in 1:d_max) {
      models[[d_try]] <- fit_tvRRR(X = X, y = y, u = u, d = d_try,
                                   model = model, Sigma = Sigma_init * diag(d_try), ...)
      crit[d_try] <- BIC_tvRRR(kf = models[[d_try]], d = d_try, model = model)

    }

    # What do I want to return?

    fit <- models[[which.min(crit)]]

    fit$BIC <- min(crit)
    fit$other_BICs <- crit
    names(fit$other_BICs) <- 1:d_max

    fit$convergence_information <- paste0(
      "Model selected based on BIC", "\n\n",
      "Tried d = 1 to d = ", d_max, "\n",
      "Rank selected: d = ", which.min(crit), " with BIC ", round(min(crit), 6), "\n\n",
      "Model diagnostics / Convergence information \n\n",
      fit$convergence_information)

    fit$rank <- which.min(crit)

    return(fit)
  }
}

## ############################################################################
##
## Model fitting function:
##
## ############################################################################

#' @rdname tvRRR
#' @export

fit_tvRRR <- function(X, y, u = NULL, d, model = "A",
                      ...) {

  # Check whether all arguments are numeric and whether the data dimensions
  # match

  if (is.null(u)) {
    stopifnot(
      "y needs to be a numeric matrix" = { length(dim(y)) == 2 & is.numeric(y) },
      "X needs to be a numeric matrix" = { length(dim(X)) == 2 & is.numeric(X) },
      "X and y need to have the same number of observations" = {nrow(X) == nrow(y)})
  } else if (!is.null(u)) {
    stopifnot(
      "y needs to be a numeric matrix" = { length(dim(y)) == 2 & is.numeric(y) },
      "X needs to be a numeric matrix" = { length(dim(X)) == 2 & is.numeric(X) },
      "u needs to be a numeric matrix" = { length(dim(u)) == 2 & is.numeric(u) },
      "X, y and u need to have the same number of observations" =
        { nrow(X) == nrow(y) & nrow(X) == nrow(u) & nrow(y) == nrow(u) }
    )
  }

  if (missing(d)) stop("The latent dimension d needs to be provided.")

  stopifnot(
    "d needs to be an integer number" = { d %% 1 == 0 },
    "d must not be larger than min(ncol(X), ncol(y))" = { d <= min(ncol(X), ncol(y)) }
  )




  if (model == "A") {

    fit <- fit_modelA(X = X, y = y, u = u, d = d, ...)

  } else if (model == "B") {

    fit <- fit_modelB(X = X, y = y, u = u, d = d, ...)

  }

  class(fit) <- c("tvRRR", class(fit))

  return(fit)
}



## ############################################################################
##
## Actual model fitting for models A and B
##
## ############################################################################


#' Runs the expectation maximization algorithm for Model A
#'
#' @param X predictors (t x q-dimensional)
#' @param y target (t x p-dimensional)
#' @param u additional predictors (t x q-dimensional)
#' @param d latent dimension
#' @param alpha_00 starting value for the algorithm, default NULL (RRR)
#' @param beta starting value for beta, default NULL (RRR)
#' @param Gamma starting value for the fixed full-rank coefficient matrix,
#'              default NULL (RRR)
#' @param P_00 starting state covariance (default 1000 * diag(p x d))
#' @param Sigma column covariance of states (default 0.01 * diag(d))
#' @param Omega error covariance (defaults to residual covariance from RRR,
#'              or respective starting values)
#' @param Omega_diagonal logical, indicates whether Omega is assumed to be a
#'                       diagonal matrix (advisable if p is large)
#' @param initialize either \code{"RRR"} or \code{"random"}, applies if no starting
#'                   values are handed over
#' @param Gamma_rrr type of normalization for the starting values obtained from RRR
#' @param return_covariances logical, indicates whether the state covariances
#'                           should be returned (might be necessary for evaluation
#'                           of the likelihood).
#' @param silent logical, indicates whether progress should be printed during
#'               model fitting
#' @param maxit maximum number of iterations for the EM algorithm
#' @param tol_finish tolerance for stopping the EM algorithm
#' @param tol_EMstep tolerance for iterative estimation during EM step
#'
#' Output:
#' @return
#' An object of class \code{tvRRR}, that is a named list of lists with elements
#' \item{states}{ The estimated states (i.e. coefficient matrices)
#'           \itemize{
#'           \item filtered (the filtered states)
#'           \item smoothed (the smoothed states)
#'           }
#'           }
#' \item{covariances}{The filtered and smoothed covariances and lag-1 covariances
#' (if \code{return_covariances = TRUE})
#' \itemize{
#'     \item \code{P_t^t} filtered covariances
#'     \item \code{P_t^t-1} predicted covariances
#'     \item \code{P_t^T} smoothed covariances
#'     \item \code{P_t-1t-2^T} smoothed lag-1 covariances
#'     }
#'     }
#' \item{data}{the data handed over to the algorithms
#' \itemize{
#'     \item \code{X} predictors
#'     \item \code{y} responses
#'     \item \code{u} additional predictors
#'     \item \code{Z} transition matrices (X_t'beta (x) I_p)
#'     }
#'     }
#' \item{parameters}{The parameters that have been fitted during the algorithm,
#' that is\itemize{
#' \item \code{Sigma} the column covariance of the states
#' \item \code{Omega} the error covariance
#' \item \code{beta} (for model A)
#' \item \code{alpha} (for model B)
#' }
#' }
#' \item{likelihoods}{list containing Q and data loglikelihood for each iteration}
#' \item{convergence_information}{Message containing parameter stages at convergence}
#'
#'
#' @export

# alpha_00 <- beta <- u <- Gamma <- Omega <- Sigma <- NULL;tol_EMstep <- tol_finish <- 1e-4; P_00 <- 1000 * diag(d*p); maxit <- 100; Omega_diagonal <- F; Gamma_rrr <- "OLS"


fit_modelA <- function(X, y, u = NULL,
                       d,
                       beta = NULL, # model parameters
                       alpha_00 = NULL,
                       Gamma = NULL,
                       P_00 = 1000 * diag(ncol(y) * d),  # starting values
                       Sigma = NULL, # column covariance
                       Omega = NULL, # measurement error covariance
                       Omega_diagonal = FALSE,
                       maxit = 100,
                       silent = FALSE,
                       tol_finish = 1e-3,
                       tol_EMstep = 1e-3,
                       return_covariances = FALSE,
                       Gamma_rrr = "identity",
                       initialize = "RRR") {

  p <- ncol(y)
  t <- nrow(X)
  q <- ncol(X)
  if (!is.null(u)) k <- ncol(u)

  if (!all(dim(P_00) == c(p * d, p * d)) |
      !isTRUE(is.symmetric(P_00)) |
      !(qr(P_00)$rank == p * d))  {
    stop("Problems with P_00.")
  }

  stopifnot("maxit needs to be at least 2" = { maxit > 1 })

  # ---------------------------------------------------------------------------
  # Starting values for the algorithm:
  # ---------------------------------------------------------------------------

  if ( is.null(beta) & is.null(alpha_00) ) {

    if (initialize == "RRR") {
    with(rr.regression(X = X, y = y, u = u, rank = d, Gamma_type = Gamma_rrr), {
      alpha_00 <<- c(A)
      beta <<- t(B)
      if (!is.null(u) & is.null(Gamma)) Gamma <<- D
    })
    } else if (initialize == "random") {
      alpha_00 <- rstiefel(p, d)
      beta <- rstiefel(q, d)
    }

  } else if (is.null(beta) & !is.null(alpha_00)) {

    beta <- rstiefel(q, d)
    if (!is.null(u) & is.null(Gamma)) Gamma <- rstiefel(p, k)

  } else if (!is.null(beta) & is.null(alpha_00)) {

    alpha_00 <- rstiefel(p, d)
    if (!is.null(u) & is.null(Gamma)) Gamma <- rstiefel(p, k)

  }

  if (all(alpha_00 == 0)) {
    stop("The algorithm cannot be initialized with zeros.
          If you did not hand over starting values there may be
          something wrong with your data.")
  }

  if (is.null(Sigma)) Sigma <- 0.1 * diag(d)

  if (is.null(Omega)) {
    if (Omega_diagonal) {
      Omega <- diag(diag(stats::cov(y - t(tcrossprod(matrix(alpha_00, p, d), beta) %*% t(X) +
                                     if (!is.null(u)) { Gamma %*% t(u) } else { 0 } ))))
    } else Omega <- stats::cov(y - t(tcrossprod(matrix(alpha_00, p, d), beta) %*% t(X) +
                                if (!is.null(u)) { Gamma %*% t(u) } else { 0 } ))
  }


  #  --------------------------------------------------------------------------
  # Model fitting: EM Algorithm
  # ---------------------------------------------------------------------------

  # Perform an initial run of the filter:
  kf <- filter_modelA(y = y, X = X, u = u,
                      beta = beta, alpha_00 = alpha_00, Gamma = Gamma,
                      P_00 = P_00, Sigma = Sigma,
                      Omega = Omega, d = d)

  Q <- rep(NA, maxit)
  loglik <- rep(NA, maxit)


  # Perform a first update -- first M-step
  update <- update_pars(kf, p = p, d = d, t = t, q = q,
                        Omega_diagonal = Omega_diagonal,
                        model = "A")

  # Initial likelihood and Q
  loglik[1] <- X_eval_lik(kf = kf, beta = beta, Omega = Omega, Gamma = Gamma,
                          d = d, y = y, X = X)
  Q[1] <- update$Q

  # Run the EM algorithm
  iter <- 1

  while (iter < maxit) {

    # Store the previous parameters:
    old_beta <- update$beta
    old_Omega <- update$Omega
    old_Sigma <- update$Sigma
    old_Gamma <- update$Gamma

    # Run the filter with the updated parameters:
    kf <- filter_modelA(y = y, X = X, u = u, Sigma = update$Sigma,
                        Omega = update$Omega, Gamma = update$Gamma,
                        alpha_00 = kf$states$smoothed[1, ],
                        P_00 = kf$covariance$`P_t^T`[1, , ],
                        d = d, beta = update$beta)

    # Update the parameters:
    update <- update_pars(kf, p = p, d = d, t = t, q = q,
                          Omega_diagonal = Omega_diagonal,
                          tol = tol_EMstep)

    iter <- iter + 1

    # Store the likelihoods from previous iteration:
    Q[iter]      <- update$Q
    loglik[iter] <- X_eval_lik(kf = kf, beta = old_beta,
                               Omega = old_Omega,
                               Gamma = old_Gamma,
                               d = d, y = y, X = X)

    # Check for problems / convergence:
    if (!is.finite(Q[iter]) | is.nan(Q[iter])) {
      cat("Q is ", Q[iter], ". Algorithm terminated.")
      break
    }
    if ((norm(old_Omega - update$Omega, "m") < tol_finish &
         subsp_dist(old_beta, update$beta) < tol_finish &
         norm(old_Sigma - update$Sigma, "m") < tol_finish &
         abs((loglik[iter] - loglik[iter - 1]) / loglik[iter - 1]) < tol_finish) |
        loglik[iter] < loglik[iter - 1]
    ) {

      break
    }

  }

  # One last run of the filter with the final parameter estimates:
  kf <- filter_modelA(y = y, X = X, u = u, Sigma = update$Sigma,
                      Gamma = update$Gamma,
                      Omega = update$Omega, alpha_00 = kf$states$smoothed[1, ],
                      P_00 = kf$covariance$`P_t^T`[1, , ],
                      d = d, beta = update$beta,
                      return_covariances = TRUE)

  # Evaluate the likelihoods:
  Q <- c(Q[1:iter],
         X_calc_Q(kf = kf, p = p, d = d, q = q, t = t, return_C = F, model = "A"))

  loglik <- c(loglik[1:iter],
              X_eval_lik(kf = kf, Omega = kf$parameters$Omega,
                         beta = kf$parameters$beta, Gamma = kf$parameters$Gamma,
                         d = d, y = y, X = X))

  if (!return_covariances) kf$covariances <- NULL

  iter <- iter + 1

  conv_message <- paste0("Algorithm stopped after ", iter,
                         " iterations. \n Relative likelihood difference: ",
                         round((loglik[iter] - loglik[iter - 1]) / abs(loglik[iter - 1]), 6),
                         "\n Difference in Omega: ", round(norm(old_Omega - update$Omega, "m"), 6),
                         "\n Difference in beta: ",  round(subsp_dist(old_beta, update$beta), 6),
                         "\n Difference in Sigma: ", round(norm(old_Sigma - update$Sigma, "m"), 6))


  if (!silent) {
    message(conv_message)
  }

  return(c(kf, iter = iter, likelihoods = list(Q = Q, loglik = loglik),
         convergence_information = conv_message))
}


#' Runs the expectation maximization algorithm for Model B
#' Input:
#' @param X predictors (t x q-dimensional)
#' @param y target (t x p-dimensional)
#' @param u additional predictors (t x q-dimensional)
#' @param d latent dimension
#' @param alpha starting value for the algorithm, default NULL (RRR)
#' @param beta_00 starting value for beta', default NULL (RRR)
#' @param Gamma starting value for the fixed full-rank coefficient matrix,
#'              default NULL (RRR)
#' @param P_00 starting state covariance (default 1000 * diag(p x d))
#' @param Sigma column covariance of states (default 0.01 * diag(d))
#' @param Omega error covariance (defaults to residual covariance from RRR,
#'              or respective starting values)
#' @param Omega_diagonal logical, indicates whether Omega is assumed to be a
#'                       diagonal matrix (advisable if p is large)
#' @param initialize either \code{"RRR"} or \code{"random"}, applies if no starting
#'                   values are handed over
#' @param Gamma_rrr type of normalization for the starting values obtained from RRR
#' @param return_covariances logical, indicates whether the state covariances
#'                           should be returned (might be necessary for evaluation
#'                           of the likelihood).
#' @param silent logical, indicates whether progress should be printed during
#'               model fitting
#' @param maxit maximum number of iterations for the EM algorithm
#' @param tol_finish tolerance for stopping the EM algorithm
#' @param tol_EMstep tolerance for iterative estimation during EM step
#'
#' Output:
#' @return
#' A named list of lists with elements
#' - states: filtered (the filtered states)
#'           smoothed (the smoothed states)
#' - covariances: the filtered and smoothed covariances and lag-1 covariances
#' (if \code{return_covariances = TRUE})
#'     `P_t^t` filtered covariances
#'     `P_t^t-1`predicted covariances
#'     `P_t^T` smoothed covariances
#'     `P_t-1t-2^T` smoothed lag-1 covariances
#' - data: the data handed over to the algorithms
#'     `X` predictors
#'     `y` responses
#'     `Z` transition matrices (X_t'beta (x) I_p)
#' - parameters used during filtering: Sigma, Omega, beta
#' - likelihoods: list containing Q and data loglikelihood for each iteration
#' - convergence_information: Message containing parameter stages at convergence
#'
#' @export

fit_modelB <- function(X, y, u = NULL, d,
                       alpha = NULL,
                       beta_00 = NULL,
                       Gamma = NULL,
                       P_00 = 1000 * diag(ncol(X) * d),
                       Sigma = NULL,
                       Omega = NULL,
                       Omega_diagonal = FALSE,
                       maxit = 100,
                       silent = FALSE,
                       tol_finish = 1e-3,
                       tol_EMstep = 1e-3,
                       return_covariances = FALSE,
                       initialize = "RRR",
                       Gamma_rrr = "identity") {

  p <- ncol(y)
  t <- nrow(y)
  q <- ncol(X)
  if (!is.null(u)) k <- ncol(u)

  if (!all(dim(P_00) == c(q * d, q * d)) |
      !isTRUE(is.symmetric(P_00)) |
      !(qr(P_00)$rank == q * d))  {
    stop("Problems with P_00.")
  }


  stopifnot("maxit needs to be at least 2" = { maxit > 1 })

  # ---------------------------------------------------------------------------
  # Starting values for the parameters:
  # ---------------------------------------------------------------------------

  if ( is.null(beta_00) & is.null(alpha) ) {

    if (initialize == "RRR") {
      with(rr.regression(X = X, y = y, u = u, rank = d, Gamma_type = Gamma_rrr), {
        alpha <<- A
        beta_00 <<- B
        if (!is.null(u) & is.null(Gamma)) Gamma <<- D
      })
    } else if (initialize == "random") {
      alpha <- rstiefel(p, d)
      beta_00 <- t(rstiefel(q, d))
      if (!is.null(u) & is.null(Gamma)) Gamma <- rstiefel(p, k)
    }

  } else if (is.null(beta_00) & !is.null(alpha)) {

    beta_00 <- t(rstiefel(q, d))
    if (!is.null(u) & is.null(Gamma)) Gamma <- rstiefel(p, k)

  } else if (!is.null(beta_00) & is.null(alpha)) {

    alpha <- rstiefel(p, d)
    if (!is.null(u) & is.null(Gamma)) Gamma <- rstiefel(p, k)

  }

  if (all(beta_00 == 0)) {
    stop("The algorithm cannot be initialized with zeros.
          If you did not hand over starting values there may be
          something wrong with your data.")
  }


  if (is.null(Sigma)) Sigma <- 0.1 * diag(d)

  if (is.null(Omega)) {
    if (Omega_diagonal) {
      Omega <- diag(diag(stats::cov(y - t(alpha %*% beta_00 %*% t(X) +
                                     if (!is.null(u)) { Gamma %*% t(u) } else { 0 } ))))
    } else Omega <- stats::cov(y - t(alpha %*% beta_00 %*% t(X) +
                                if (!is.null(u)) { Gamma %*% t(u) } else { 0 }))
  }

  # ---------------------------------------------------------------------------
  # Model fitting: Run the EM-algorithm
  # ---------------------------------------------------------------------------

  # Initial run of the filter
  kf <- filter_modelB(X = X, y = y, u = u,
                      beta_00 = beta_00,
                      alpha = alpha, Gamma = Gamma,
                      P_00 = P_00, Sigma = Sigma,
                      Omega = Omega, d = d)

  Q <- rep(NA, maxit)
  loglik <- rep(NA, maxit)

  update <- update_pars(kf, p = p, d = d, t = t, q = q,
                        Omega_diagonal = Omega_diagonal,
                        model = "B")

  # Likelihoods at start of algorithm:
  loglik[1] <- X_eval_lik(kf = kf, Omega = Omega, alpha = alpha,
                          d = d, X = X, y = y)
  Q[1] <- update$Q

  # Run the EM-algorithm:
  iter <- 1

  while (iter < maxit) {

    old_alpha <- update$alpha; old_Omega <- update$Omega; old_Gamma <- update$Gamma
    old_Sigma <- update$Sigma

    # E-step:
    kf <- filter_modelB(y = y, X = X, u = u, Gamma = update$Gamma, Sigma = update$Sigma,
                        Omega = update$Omega, beta_00 = kf$states$smoothed[1, ],
                        P_00 = kf$covariance$`P_t^T`[1, , ],
                        d = d, alpha = update$alpha)

    # M-step:
    update <- update_pars(kf, p = p, d = d, t = t, q = q,
                          Omega_diagonal = Omega_diagonal, # calc_Q = TRUE,
                          tol = tol_EMstep, model = "B")

    iter <- iter + 1

    Q[iter] <- update$Q
    loglik[iter] <- X_eval_lik(kf = kf, alpha = old_alpha, Omega = old_Omega,
                               Gamma = old_Gamma,
                               d = d, y = y, X = X)

    # Check for convergence / errors:
    if (!is.finite(Q[iter]) | is.nan(Q[iter])) {
      cat("Q is ", Q[iter], ". Algorithm terminated.")
      break
    }

    # Check conditions for stopping
    if ((norm(old_Omega - update$Omega, "m") < tol_finish &
         subsp_dist(old_alpha, update$alpha) < tol_finish &
         norm(old_Sigma - update$Sigma, "m") < tol_finish &
         abs((loglik[iter] - loglik[iter - 1]) / loglik[iter - 1]) < tol_finish) |
        loglik[iter] < loglik[iter - 1]
    ) {

      break
    }

  }

  # One last run of the filter with the final set of parameters:
  kf <- filter_modelB(y = y, X = X, u = u, Gamma = update$Gamma, Sigma = update$Sigma,
                      Omega = update$Omega, beta_00 = kf$states$smoothed[1, ],
                      P_00 = kf$covariance$`P_t^T`[1, , ],
                      d = d, alpha = update$alpha,
                      return_covariances = TRUE)

  # Evaluate the likelihoods for this set of parameters:
  Q <- c(Q[1:iter], X_calc_Q(kf = kf, p = p, d = d, q = q, t = t, return_C = F, model = "B"))
  loglik <- c(loglik[1:iter], X_eval_lik(kf = kf, Omega = kf$parameters$Omega,
                                        alpha = kf$parameters$alpha, Gamma = kf$parameters$Gamma,
                                        d = d, y = y, X = X))

  iter <- iter + 1

  conv_message <-  paste0("Algorithm terminated after ", iter,
                          " iterations. \n Relative likelihood difference: ",
                          round((loglik[iter] - loglik[iter - 1]) / abs(loglik[iter - 1]), 6),
                          "\n Difference in Omega: ", round(norm(old_Omega - update$Omega, "m"), 6),
                          "\n Difference in beta: ",  round(subsp_dist(old_alpha, update$alpha), 6),
                          "\n Difference in Sigma: ", round(norm(old_Sigma - update$Sigma, "m"), 6))


  if (!silent) {
    message(conv_message)
  }


  if (!return_covariances) kf$covariances <- NULL

  return(c(kf, iter = iter, likelihoods = list(Q = Q, loglik = loglik),
           convergence_information = conv_message))
}


