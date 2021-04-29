
## #############################################################################
##
## Fitted- and print methods for objects of class tvRRR
##
## #############################################################################

#' Calculate fitted values from tvRRR object
#'
#' @param object an object of class \code{tvRRR}
#' @param type either \code{filtered} or \code{smoothed}, defaults to filtered
#' @param model \code{"A"} or \code{"B"}, can be handed over, but otherwise is determined
#' automatically from the \code{tvRRR} object
#' @param ... ignored
#'
#' @export

fitted.tvRRR <- function(object, type = "filtered", model, ...) {

  kf <- object

  if (missing(model)) model <- ifelse(is.null(kf$parameters$alpha), "A", "B")

  t <- nrow(kf$states$filtered) - 1
  p <- ncol(kf$data$y)
  q <- ncol(kf$data$X)
  d <- ifelse(model == "A", ncol(kf$states$filtered) / p, ncol(kf$states$filtered) / q)
  u <- kf$data$u


  if (type == "filtered") {
    fitted <- switch(
      model,
      A = t(sapply(1:t, function(i) matrix(kf$states$filtered[i + 1, ], p, d) %*% t(kf$parameters$beta) %*% kf$data$X[i, ] +
                     if (!is.null(kf$data$u)) kf$parameters$Gamma %*% u[i, ] else 0)),
      B = t(sapply(1:t, function(i) kf$parameters$alpha %*%  matrix(kf$states$filtered[i + 1, ], d, q) %*% kf$data$X[i, ] +
                     if (!is.null(kf$data$u)) kf$parameters$Gamma %*% u[i, ] else 0))
    )
  } else if (type == "smoothed") {
    fitted <- switch(
      model,
      A = t(sapply(1:t, function(i) matrix(kf$states$smoothed[i + 1, ], p, d) %*% t(kf$parameters$beta) %*% kf$data$X[i, ] +
                     if (!is.null(kf$data$u)) kf$parameters$Gamma %*% u[i, ] else 0)),
      B = t(sapply(1:t, function(i) kf$parameters$alpha %*%  matrix(kf$states$smoothed[i + 1, ], d, q) %*% kf$data$X[i, ]+
                     if (!is.null(kf$data$u)) kf$parameters$Gamma %*% u[i, ] else 0))
    )
  }

  return(fitted)
}


#' Print Method for tvRRR object
#'
#' @param x an object of class \code{tvRRR}
#' @param ... ignored
#'
#' @returns A short summary of the \code{tvRRR} object.
#'
#' @export
print.tvRRR <- function(x, ...) {

  kf <- x

  d <- ifelse(is.null(kf$parameters$alpha), dim(kf$parameters$beta)[2], dim(kf$parameters$alpha)[2])
  model <- ifelse(is.null(kf$parameters$alpha), "A", "B")

  cat("Time-varying reduced rank regression model of type", model,
      "with latent rank d =", d, "\n\n")

  cat("Convergence information: \n\n")

  cat(kf$convergence_information, "\n\n")

  cat("Estimated state covariance matrix: \n\n")

  print(kf$parameters$Sigma)

  cat("\n",
      "Proportion of variance explained per time series: \n")

  vars <- colMeans((fitted.tvRRR(kf) - kf$data$y)^2) / apply(kf$data$y, 2, stats::var)

  if (!is.null(colnames(kf$data$y))) names(vars) <- colnames(kf$data$y)

  print(t(vars))
}

#' Predict method for tvRRR object
#' @param object object of class tvRRR as returned by \code{\link[tvRRR]{tvRRR}()} function
#' @param newdata can be specified differently with differing behavior: \itemize{
#'                \item if \code{NULL}, the \code{\link[tvRRR]{fitted.tvRRR}()} method is called and
#'                fitted values for the initial data are calculated
#'                \item a named list of length 1 or 2 containing predictors X and u for predictions
#'                      using the last state fitted
#'                \item a named list of length 2 or 3 containing targets y, predictors X and, optionally
#'                      additional predictors u; here the estimated states are updated in every step
#'                      of the algorithm
#'                }
#' @param silent specifies whether a message should be printed regarding the type of prediction
#' @param ... ignored
#'
#'
#' @returns The predicted target variables in the specified setup (as a matrix)
#'
#' @export


predict.tvRRR <- function(object, newdata = NULL, silent = FALSE, ...) {

  kf <- object

  if (is.null(newdata)) return(fitted.tvRRR(kf))

  # Unterscheide zwischen Dimension 1 und höherer Dimension für die predictions

  model <- ifelse(is.null(kf$parameters$alpha), "A", "B")

  if (!is.null(newdata)) {

    if ( !(all( names(newdata) %in% c("y", "X", "u") )) ) stop("newdata needs to be a named list.")

    X <- newdata$X
    y <- newdata$y
    u <- newdata$u


    if (is.null(y)) {

      t <- nrow(kf$data$X)
      q <- ncol(X)

      if (model == "A") {
        yhat <- t(matrix(kf$states$filtered[t + 1, ], p, d) %*% t(kf$parameters$beta) %*% t(X) +
                    if (!is.null(u)) kf$parameters$Gamma %*% t(u) else 0
        )
      } else if (model == "B") {
        yhat <- t(kf$parameters$alpha %*% matrix(kf$states$filtered[t + 1, ], d, q) %*% t(X) +
                    if (!is.null(u)) kf$parameters$Gamma %*% t(u) else 0
        )
      }
      if (!silent) message("Prediction without updating.")

      return(yhat)

    }

    if (!is.null(y)) {

        p <- ncol(y)
        q <- ncol(X)
        t <- nrow(kf$data$y)
        # if(!is.null(u)) k <- ncol(u)

        d <- ifelse(model == "A", dim(kf$parameters$beta)[2], dim(kf$parameters$alpha)[2])

      alpha <- {
        if (model == "A") {
          matrix(kf$states$smoothed[t + 1, ], p, d) ## BRAUCHE ICH HIER t ODER t+1???
        } else if (model == "B") {
          kf$parameters$alpha
        }
      }

      beta <- {
        if (model == "B") {
          t(matrix(kf$states$smoothed[t + 1, ], d, q))
        } else if (model == "A") {
          kf$parameters$beta
        }
      }

      new <- eval_tvRRR(X = X, y = y, u = u, model = model,
                        alpha = alpha,
                        beta = beta,
                        Gamma = kf$parameters$Gamma,
                        Omega = kf$parameters$Omega,
                        Sigma = kf$parameters$Sigma,
                        P_00 = kf$prediction_covariance,
                        d = d)

      if (model == "A") {
        yhat <- t(sapply(1:nrow(y), function(i) {
          matrix(new$states$onestepahead[i, ], p, d) %*% t(new$parameters$beta) %*% X[i, ] +
            if (!is.null(u)) new$parameters$Gamma %*% u[i, ] else 0
        }))
      } else if (model == "B") {
        yhat <- t(sapply(1:nrow(y), function(i) {
          new$parameters$alpha %*% matrix(new$states$onestepahead[i, ], d, q) %*% X[i, ] +
            if (!is.null(u)) new$parameters$Gamma %*% u[i, ] else 0
        }))
      }

      if (!silent) message("Prediction with subsequent state updates, one-step-ahead.")

      return(yhat)
    }
  }
}


## ############################################################################
##
## BIC information criterion
##
## ############################################################################

#' @keywords internal
BIC_tvRRR <- function(kf, d, model = "A", ...) {

  t <- nrow(kf$data$X)

  if (model == "A") {
    K <- d * (d + 1) / 2 + ifelse(d > 1, prod(dim(kf$parameters$beta)), length(kf$parameters$beta))
  } else if (model == "B") {
    K <- d * (d + 1) / 2 + ifelse(d > 1, prod(dim(kf$parameters$alpha)), length(kf$parameters$alpha))
  }

  -2 * kf$likelihoods.loglik[kf$iter] + K * log(t)
}

#' @keywords internal
MSFE_tvRRR <- function(kf, newdata = NULL, model = "A") {
  yhat <- stats::predict(kf, newdata = newdata)
  mean((newdata$y - yhat)^2)
}

#' Corrected version (for small samples) taken from
#' Burnham, Anderson (2004): Multimodel Inference: Understanding AIC and BIC in model selection
#' Recommended if t / K < 40, but as it converges to zero with growing t it should be used anyway
#' for that correction is set as default
# #' @export
# AIC.tvRRR <- function(kf, d, model = "A", correct = TRUE, ...) {
#
#   t <- nrow(kf$data$X)
#
#   if (model == "A") {
#     K <- d * (d + 1) / 2 + ifelse(d > 1, prod(dim(kf$parameters$beta)), length(kf$parameters$beta))
#   } else if (model == "B") {
#     K <- d * (d + 1) / 2 + ifelse(d > 1, prod(dim(kf$parameters$alpha)), length(kf$parameters$alpha))
#   }
#
#   if (!correct) {
#     return(-2 * kf$likelihoods.logLik[kf$iter] + K * 2)
#   } else {
#     return(-2 * kf$likelihoods.loglik[kf$iter] + 2 * K + 2 * K * (K + 1) / (t - K - 1))
#   }
# }
