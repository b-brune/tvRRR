
#' DATA GENERATING FUNCTIONS:
#' This needs to be documented very carefully
#' @param type \code{= c("deterministic", "break", "VARbreak", "rw", "fixed")},
#'             see Details for a description of the different kinds
#' @param p dimension of the target vector
#' @param q dimension of the predictors
#' @param d latent dimension of the model
#' @param t length of the time series
#' @param X optional, the predictors used to generate the time series.
#'          Not used if \code{type = "VARbreak"}. If not handed over, it is drawn
#'          randomly from a q-variate normal distribution with mean zero and
#'          covariance identity
#' @param model which model do we want? Either "A" or "B"
#' @param forecast how many additional variables do we want to generate for
#'                 forecasting experiments
#' @param ... additional variables specific to the different types, see Details.
#'
#' @returns A named list with the desired dataset, containing components
#' \item{y}{the target variable (t x p) matrix}
#' \item{X}{the predictors (t x q) matrix}
#' \item{alpha}{if \code{model = "A"}, an array of dimensions (t, p, d),
#'              if \code{model = "B"}, a (p, d) matrix}
#' \item{beta}{if \code{model = "B"}, an array of dimensions (t, q, d),
#'              if \code{model = "A"}, a (q, d) matrix}
#' \item{Omega}{the error covariance (p, p) dimensional}
#' \item{Sigma}{if \code{type = "rw"} the column covariance of the states, (d, d) matrix}
#'
#' @details
#'
#' Optional parameters that may be handed over to the function are
#' \itemize{
#' \item `alpha_0` optional
#' \item `alpha_1` optional (only used if \code{model = "A"} and
#'                \code{type = c("break", "VARbreak", "deterministic")} )
#' \item `beta_0` optional
#' \item `beta_1` optional (only used if \code{model = "B"} and
#'                \code{type = c("break", "VARbreak", "deterministic")} )
#' \item `Omega` error covariance, defaults to the identity matrix
#' \item `shift` only used if \code{type = c("break", "VARbreak", "deterministic")}),
#'              guides how far the matrices are apart (see Details)
#' \item `Sigma` column covariance if \code{type = "rw"}
#' \item `Delta` row covariance if \code{type = "rw"}, defaults to the identity
#' \item `breakpoint` position of the structural break in the dataset for
#'                   \code{type = c("break", "VARbreak")}

#'
#' }
#'
#' @export

# The transition matrices are generated as from the following four mechanisms:
#
# Fixed: $\balpha_t \equiv \balpha \forall t=1,...,T$, or respectively $\bbeta_t \equiv \bbeta \forall t=1,...,T$; $\balpha,\bbeta$ are randomly drawn orthogonal matrices.
#
# Deterministic: deterministic transition between two matrices $\balpha_1$ and $\balpha_2$ with shift parameter $s$,
#  i.e.
# $\balpha_t = -s \frac{T-t}{T-1} \balpha_1 + s \frac{t-1}{T-1} \balpha_2$
# The shift controls how strong the parameter difference is.
#
#  Random walk: The matrices are generated as a random walk according to the actual model that we assume during the filtering procedure, i.e.
# $\balpha_t \sim \norm_{p,d}(\balpha_{t-1}; \id_p,\Sigma_c)$
# where $\Sigma_c\in\R^{d\times d}$ is drawn randomly from the Wishart distribution.
#
# Structural break: We generate
# $\balpha_t = \begin{cases}  -s \cdot \alpha_1, & t \leq T/2   \\
# s \cdot \alpha_2, & t > T/2
# \end{cases} $
# $s$ is specified as in the deterministic model; $\alpha_1$ and $\alpha_2$ are randomly drawn orthogonal matrices.
#
# The corresponding data set is then drawn according to the model equation
# $$ \by_t = \balpha_t\bbeta'\bx_t + \eps_t $$
# where $(\bx_t)_{t=1,...,T}$ are $q$-dimensional vectors of independent $\text{MA}(2)$ processes with parameters $(0.9, 0.5)$.
# The errors are drawn from a multivariate normal distribution with covariance matrix $\Omega$ which is either diagonal with $U[0.5, 1.5]$ entries, or randomly drawn from the Wishart distribution.
# In addition, we generate a \enquote{validation dataset} by generating observations $T+1,...,T+f$, $f=T/8$ according to above models.




dataset <- function(type, p, d, q, t, X, forecast, model = "A", ...) {

  if (!missing(forecast)) len <- t + forecast else { len <- t; forecast <- 0 }

  if (missing(X)) {
    X <- mvtnorm::rmvnorm(n = len, mean = rep(0, q), sigma = diag(q))
  } else if (!missing(X)) {
    t <- nrow(X)
    q <- ncol(X)
  }

  dat <- switch(
    type,
    deterministic = make_dataset_deterministic(X = X, p = p, d = d, model = model, forecast = forecast, ...),
    `break` = make_dataset_break(X = X, p = p, d = d, model = model, forecast = forecast, ...),
    VARbreak = make_VARdataset_break(t = t, p = p, d = d, model = model, forecast = forecast, ...),
    rw = make_dataset_rw(X = X, p = p, d = d, model = model, forecast = forecast, ...),
    fixed = make_dataset_fixed(X = X, p = p, d = d, model = model, forecast = forecast, ...),
    VARdeterm = make_VARdataset_determ(t = t, p = p, d = d, model = model, forecast = forecast, ...)
  )

  return(dat)
}


#' Dataset that has a deterministic transition
#' @keywords internal

make_dataset_deterministic <- function(X, alpha_0 = NULL, alpha_1 = NULL,
                                       shift = 1, beta_0 = NULL, beta_1 = NULL,
                                       Omega = NULL,
                                       p, d, model = "A", forecast = 0,
                                       ...) {

  # Initializiation of parameters
  if (is.null(alpha_0)) {
    if (is.null(p) | is.null(d)) {
      stop("Either alpha_0 or its dimensions p and d need to be provided.")
    } else if (!is.null(p) & !is.null(d)) {
      alpha_0 <- rstiefel(p, d)
    }
  }

  if (model == "A") {
    if (is.null(alpha_1)) {
      if (is.null(p) | is.null(d)) {
        stop("Either alpha_1 or its dimensions p and d need to be provided.")
      } else if (!is.null(p) & !is.null(d)) {
        alpha_1 <- rstiefel(p, d)
      }
    }
  }

  if (is.null(beta_0)) {
    if (is.null(d)) {
      stop("Either beta_0 or its dimension d need to be provided.")
    } else if (!is.null(d)) {
      beta_0 <- rstiefel(ncol(X), d)
    }
  }

  if (model == "B") {
    if (is.null(beta_1)) {
      if (is.null(d)) {
        stop("Either beta_1 or its dimension d needs to be provided.")
      } else if (!is.null(d)) {
        beta_1 <- rstiefel(ncol(X), d)
      }
    }
  }


  if (is.null(Omega)) {
    if (is.null(p)) {
      stop ("p must be provided.")
    } else {
      Omega <- diag(p)
    }
  }


  ## Observation 1:t used for model fitting, additionally we generate
  ## forecast observation that we want to forecast
  t <- nrow(X) - forecast
  len <- t + forecast

  if (model == "A") {

    alpha <- lapply(1:len, function(i) -shift * (t - i) /(t - 1) * alpha_0 +
                      shift * (i - 1) / (t - 1) * alpha_1)

    alpha <- aperm(array(unlist(alpha), dim = c(p, d, len)), perm = c(3, 1, 2))

    y <- t(sapply(1:len, function(i) alpha[i, , ] %*% t(beta_0) %*% X[i, ] +
                  t(mvtnorm::rmvnorm(1, sigma = Omega))))

    beta <- beta_0

  } else if (model == "B") {

    q <- ncol(X)

    beta <- lapply(1:len, function(i) -shift * (t - i) /(t - 1) * beta_0 +
                     shift * (i - 1) / (t - 1) * beta_1)

    beta <- aperm(array(unlist(beta), dim = c(q, d, len)), perm = c(3, 1, 2))

    y <- t(sapply(1:len, function(i) alpha_0 %*% t(beta[i, , ]) %*% X[i, ] +
                  t(mvtnorm::rmvnorm(1, sigma = Omega))))

    alpha <- alpha_0
  }

  return(list(y = y, X = X, alpha = alpha, beta = beta, Omega = Omega))
}

#' VAR-dataset with deterministic transition
#' @keywords internal

make_VARdataset_determ <- function(t = 100,
                                  alpha_0 = NULL, alpha_1 = NULL, beta_0 = NULL,
                                  beta_1 = NULL, Omega = NULL, p, d, model = "A",
                                  shift = 1,
                                  forecast = 0, lag = 1, breakpoint = floor(t/2)) {

  if (lag != 1) stop("Sorry, lags higher than one haven't been implemented yet.")

  if (is.null(alpha_0)) {
    if (is.null(p) | is.null(d)) {
      stop("Either alpha_0 or its dimensions p and d need to be provided.")
    } else if (!is.null(p) & !is.null(d)) {
      alpha_0 <- rstiefel(p, d)
      alpha_1 <- rstiefel(p, d)

    }
  }

  if (is.null(beta_0)) {
    if (is.null(d)) {
      stop("Either beta or its dimension d need to be provided.")
    } else if (!is.null(d)) {
      beta_0 <- rstiefel(lag * p, d)
      beta_1 <- rstiefel(lag * p, d)
    }
  }

  if (is.null(Omega)) {
    if (is.null(p)) {
      stop ("p must be provided.")
    } else {
      Omega <- diag(p)
    }
  }

  len <- t + forecast

  if (model == "A") {

    alpha <- lapply(1:len, function(i) -shift * (t - i) /(t - 1) * alpha_0 +
                      shift * (i - 1) / (t - 1) * alpha_1)

    alpha <- aperm(array(unlist(alpha), dim = c(p, d, len)), perm = c(3, 1, 2))

    y <- matrix(NA, len + 1, p)
    y[1, ] <- stats::rnorm(p)


    errors <- mvtnorm::rmvnorm(len, sigma = Omega)


    for (i in 2:(len + 1)) {
      y[i, ] <- alpha[i - 1, , ] %*% t(beta_0) %*% y[i - 1, ] + errors[i -1, ]
    }

    X <- y[1:len, ]
    y <- y[-1, ]

    return(list(y = y, X = X, errors = errors, alpha = alpha, beta = beta_0, Omega = Omega))
  } else if (model == "B") {

    q <- lag * p

    beta <- lapply(1:len, function(i) -shift * (t - i) /(t - 1) * beta_0 +
                     shift * (i - 1) / (t - 1) * beta_1)

    beta <- aperm(array(unlist(beta), dim = c(q, d, len)), perm = c(3, 1, 2))

    y <- matrix(NA, len + 1, p)
    y[1, ] <- stats::rnorm(p)

    beta <- lapply(1:len, function(i) if (i <= breakpoint) -shift * beta_0 else shift * beta_1)

    beta <- aperm(array(unlist(beta), dim = c(q, d, len)), perm = c(3, 1, 2))

    errors <- mvtnorm::rmvnorm(len, sigma = Omega)

    for (i in 2:(len + 1)) {
      y[i, ] <- alpha_0 %*% t(beta[i - 1, , ]) %*% y[i - 1, ] + errors[i - 1, ]
    }

    X <- y[1:len, ]
    y <- y[-1, ]


    return(list(y = y, X = X, errors = errors, alpha = alpha_0, beta = beta, Omega = Omega))
  }

}




#' VAR-dataset that contains a break
#' @keywords internal

make_VARdataset_break <- function(t = 100,
                                  alpha_0 = NULL, alpha_1 = NULL, beta_0 = NULL,
                                  beta_1 = NULL, Omega = NULL, p, d, model = "A",
                                  shift = 1,
                                  forecast = 0, lag = 1, breakpoint = floor(t/2)) {

  if (lag != 1) stop("Sorry, lags higher than one haven't been implemented yet.")

  if (model == "A") {
    if (is.null(alpha_0)) {
      if (is.null(p) | is.null(d)) {
        stop("Either alpha_0 or its dimensions p and d need to be provided.")
      } else if (!is.null(p) & !is.null(d)) {
        alpha_0 <- rstiefel(p, d)
      }
    }

    if (is.null(alpha_1)) {
      if (is.null(p) | is.null(d)) {
        stop("Either alpha_1 or its dimensions p and d need to be provided.")
      } else if (!is.null(p) & !is.null(d)) {
        alpha_1 <- rstiefel(p, d)
      }
    }

    if (is.null(beta_0)) {
      if (is.null(d)) {
        stop("Either beta or its dimension d need to be provided.")
      } else if (!is.null(d)) {
        beta_0 <- rstiefel(lag * p, d)
      }
    }

    if (is.null(Omega)) {
      if (is.null(p)) {
        stop ("p must be provided.")
      } else {
        Omega <- diag(p)
      }
    }

    len <- t + forecast

    alpha <- lapply(1:len, function(i) if (i <= breakpoint) -shift * alpha_0 else shift * alpha_1)
    alpha <- aperm(array(unlist(alpha), dim = c(p, d, len)), perm = c(3, 1, 2))

    y <- matrix(NA, len + 1, p)
    y[1, ] <- stats::rnorm(p)


    errors <- mvtnorm::rmvnorm(len, sigma = Omega)


    for (i in 2:(len + 1)) {
      y[i, ] <- alpha[i - 1, , ] %*% t(beta_0) %*% y[i - 1, ] + errors[i -1, ]
    }

    X <- y[1:len, ]
    y <- y[-1, ]

    return(list(y = y, X = X, errors = errors, alpha = alpha, beta = beta_0, Omega = Omega))
  }

  if (model == "B") {
    if (is.null(alpha_0)) {
      if (is.null(p) | is.null(d)) {
        stop("Either alpha_0 or its dimensions p and d need to be provided.")
      } else if (!is.null(p) & !is.null(d)) {
        alpha_0 <- rstiefel(p, d)
      }
    }

    if (is.null(beta_0)) {
      if (is.null(d)) {
        stop("Either beta or its dimension d need to be provided.")
      } else if (!is.null(d)) {
        beta_0 <- rstiefel(p * lag, d)
      }
    }

    if (is.null(beta_1)) {
      if (is.null(d)) {
        stop("Either beta or its dimension d need to be provided.")
      } else if (!is.null(d)) {
        beta_1 <- rstiefel(p * lag, d)
      }
    }

    if (is.null(Omega)) {
      if (is.null(p)) {
        stop ("p must be provided.")
      } else {
        Omega <- diag(p)
      }
    }

    len <- t + forecast
    q <- p * lag

    y <- matrix(NA, len + 1, p)
    y[1, ] <- stats::rnorm(p)

    beta <- lapply(1:len, function(i) if (i <= breakpoint) -shift * beta_0 else shift * beta_1)

    beta <- aperm(array(unlist(beta), dim = c(q, d, len)), perm = c(3, 1, 2))

    errors <- mvtnorm::rmvnorm(len, sigma = Omega)

    for (i in 2:(len + 1)) {
      y[i, ] <- alpha_0 %*% t(beta[i - 1, , ]) %*% y[i - 1, ] + errors[i - 1, ]
    }

    X <- y[1:len, ]
    y <- y[-1, ]


    return(list(y = y, X = X, errors = errors, alpha = alpha_0, beta = beta, Omega = Omega))
  }

}



#' Dataset that contains a structural break
#' @keywords internal
make_dataset_break <- function(X, alpha_0 = NULL, alpha_1 = NULL, shift = 1,
                               beta_0 = NULL, beta_1 = NULL, Omega = NULL,
                               p, d, model = "A", forecast = 0,
                               breakpoint = floor(nrow(X)/2), ...) {


  q <- ncol(X)

  if (model == "A") {
    if (is.null(alpha_0)) {
      if (is.null(p) | is.null(d)) {
        stop("Either alpha_0 or its dimensions p and d need to be provided.")
      } else if (!is.null(p) & !is.null(d)) {
        alpha_0 <- rstiefel(p, d)
      }
    }

    if (is.null(alpha_1)) {
      if (is.null(p) | is.null(d)) {
        stop("Either alpha_1 or its dimensions p and d need to be provided.")
      } else if (!is.null(p) & !is.null(d)) {
        alpha_1 <- rstiefel(p, d)
      }
    }

    if (is.null(beta_0)) {
      if (is.null(d)) {
        stop("Either beta or its dimension d need to be provided.")
      } else if (!is.null(d)) {
        beta_0 <- rstiefel(q, d)
      }
    }

    if (is.null(Omega)) {
      if (is.null(p)) {
        stop ("p must be provided.")
      } else {
        Omega <- diag(p)
      }
    }

    t <- nrow(X) - forecast
    len <- t + forecast

    alpha <- lapply(1:len, function(i) if (i <= breakpoint) -shift * alpha_0 else shift * alpha_1)

    alpha <- aperm(array(unlist(alpha), dim = c(p, d, len)), perm = c(3, 1, 2))

    y <- t(sapply(1:len, function(i) alpha[i, , ] %*% t(beta_0) %*% X[i, ] +
                  t(mvtnorm::rmvnorm(1, sigma = Omega))))

    return(list(y = y, X = X, alpha = alpha, beta = beta_0, Omega = Omega))
  }

  if (model == "B") {
    if (is.null(alpha_0)) {
      if (is.null(p) | is.null(d)) {
        stop("Either alpha_0 or its dimensions p and d need to be provided.")
      } else if (!is.null(p) & !is.null(d)) {
        alpha_0 <- rstiefel(p, d)
      }
    }

    if (is.null(beta_0)) {
      if (is.null(d)) {
        stop("Either beta or its dimension d need to be provided.")
      } else if (!is.null(d)) {
        beta_0 <- rstiefel(ncol(X), d)
      }
    }

    if (is.null(beta_1)) {
      if (is.null(d)) {
        stop("Either beta or its dimension d need to be provided.")
      } else if (!is.null(d)) {
        beta_1 <- rstiefel(ncol(X), d)
      }
    }

    if (is.null(Omega)) {
      if (is.null(p)) {
        stop ("p must be provided.")
      } else {
        Omega <- diag(p)
      }
    }

    t <- nrow(X) - forecast
    len <- t + forecast
    q <- ncol(X)

    beta <- lapply(1:len, function(i) if (i <= (t / 2)) -shift * beta_0 else shift * beta_1)

    beta <- aperm(array(unlist(beta), dim = c(q, d, len)), perm = c(3, 1, 2))

    y <- t(sapply(1:len, function(i) alpha_0 %*% t(beta[i, , ]) %*% X[i, ] +
                  t(mvtnorm::rmvnorm(1, sigma = Omega))))

    return(list(y = y, X = X, alpha = alpha_0, beta = beta, Omega = Omega))
  }
}

#' Dataset where the time-varying parameter matrices perform a random walk in
#' euclidean space where Cov(vec(alpha_t)) = kronecker(Sigma, Delta), Sigma
#' is dxd (column covariance) and Delta is pxp (row covariance)
#' @keywords internal

make_dataset_rw <- function(X, alpha_0 = NULL, beta_0 = NULL, Omega = NULL,
                            p, d, q = ncol(X), Sigma = NULL,
                            Delta = if (model == "A") diag(p) else diag(q),
                            model = "A", forecast = 0,
                            ...) {

  if (is.null(Omega)) {
    if (is.null(p)) {
      stop ("p must be provided.")
    } else {
      Omega <- diag(p)
    }
  }

  t <- nrow(X) - forecast
  len <- t + forecast

  q <- ncol(X)

  if (is.null(Sigma)) {
    Sigma <- 1 / t * drop(stats::rWishart(1, d, Sigma = 1 / q * diag(d)))
  }


  rw <- switch(model,
               A = mvtnorm::rmvnorm(len, sigma = kronecker(Sigma, Delta)),
               B = mvtnorm::rmvnorm(len, sigma = kronecker(Delta, Sigma)))

  eps <- mvtnorm::rmvnorm(len, sigma = Omega)


  if (model == "A") {
    if (is.null(alpha_0)) {
      if (is.null(p) | is.null(d)) {
        stop("Either alpha_0 or its dimensions p and d need to be provided.")
      } else if (!is.null(p) & !is.null(d)) {
        alpha_0 <- rstiefel(p, d)
      }
    }

    if (is.null(beta_0)) {
      if (is.null(d)) {
        stop("Either beta or its dimension d need to be provided.")
      } else if (!is.null(d)) {
        beta_0 <- rstiefel(ncol(X), d)
      }
    }

    alpha <- apply(rbind(c(alpha_0), rw), 2, cumsum)[-1, ]

    y <- t(sapply(1:len,
                function(i) matrix(alpha[i, ], p, d) %*% t(beta_0) %*%
                  X[i, ] + eps[i, ]))

    alpha <- t(alpha)
    dim(alpha) <- c(p, d, len)
    alpha <- aperm(alpha, c(3, 1, 2))
  }

  if (model == "B") {

    if (is.null(alpha_0)) {
      if (is.null(p) | is.null(d)) {
        stop("Either alpha or its dimensions p and d need to be provided.")
      } else if (!is.null(p) & !is.null(d)) {
        alpha_0 <- rstiefel(p, d)
      }
    }

    if (is.null(beta_0)) {
      if (is.null(d)) {
        stop("Either beta or its dimension d need to be provided.")
      } else if (!is.null(d)) {
        beta_0 <- rstiefel(ncol(X), d)
      }
    }

    betas <- apply(rbind(c(beta_0), rw), 2, cumsum)[-1, ]

    y <- t(sapply(1:len,
                function(i) {
                  alpha_0 %*% t(matrix(betas[i, ], q, d)) %*% X[i, ] + eps[i, ]
                  }))

    betas <- t(betas)
    dim(betas) <- c(q, d, len)
    betas <- aperm(betas, c(3, 1, 2))
  }

  return(list(y = y, X = X,
              alpha = if (model == "A") alpha else alpha_0,
              beta = if (model == "B") betas else beta_0, Omega = Omega,
              Sigma = Sigma, Delta = Delta))
}


#' Matrices are actually fixed
#' @keywords internal
make_dataset_fixed <- function(X, alpha = NULL, beta = NULL, Omega = NULL,
                               p, d, model = "A", forecast = 0) {

  if (is.null(alpha)) alpha <- rstiefel(p, d)
  if (is.null(beta)) beta <- rstiefel(ncol(X), d)
  if (is.null(Omega)) Omega <- diag(p)

  t <- nrow(X) - forecast
  len <- t + forecast

  q <- ncol(X)

  eps <- mvtnorm::rmvnorm(len, sigma = Omega)

  y <- t(sapply(1:len,
              function(i) alpha %*% t(beta) %*% X[i, ] + eps[i, ]))

  if (model == "A") alpha <- aperm(array(alpha, dim = c(p, d, len)), c(3, 1, 2))

  if (model == "B") beta <- aperm(array(beta, dim = c(q, d, len)), c(3, 1, 2))

  return(list(y = y, X = X, alpha = alpha, beta = beta, Omega = Omega))
}
