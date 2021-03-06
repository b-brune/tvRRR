% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FITTING.R
\name{fit_modelB}
\alias{fit_modelB}
\title{Runs the expectation maximization algorithm for Model B
Input:}
\usage{
fit_modelB(
  X,
  y,
  u = NULL,
  d,
  alpha = NULL,
  beta_00 = NULL,
  Gamma = NULL,
  P_00 = 1000 * diag(ncol(X) * d),
  Sigma = NULL,
  Omega = NULL,
  Omega_diagonal = FALSE,
  maxit = 100,
  silent = FALSE,
  tol_finish = 0.001,
  tol_EMstep = 0.001,
  return_covariances = FALSE,
  initialize = "RRR",
  Gamma_rrr = "identity"
)
}
\arguments{
\item{X}{predictors (t x q-dimensional)}

\item{y}{target (t x p-dimensional)}

\item{u}{additional predictors (t x q-dimensional)}

\item{d}{latent dimension}

\item{alpha}{starting value for the algorithm, default NULL (RRR)}

\item{beta_00}{starting value for beta', default NULL (RRR)}

\item{Gamma}{starting value for the fixed full-rank coefficient matrix,
default NULL (RRR)}

\item{P_00}{starting state covariance (default 1000 * diag(p x d))}

\item{Sigma}{column covariance of states (default 0.01 * diag(d))}

\item{Omega}{error covariance (defaults to residual covariance from RRR,
or respective starting values)}

\item{Omega_diagonal}{logical, indicates whether Omega is assumed to be a
diagonal matrix (advisable if p is large)}

\item{maxit}{maximum number of iterations for the EM algorithm}

\item{silent}{logical, indicates whether progress should be printed during
model fitting}

\item{tol_finish}{tolerance for stopping the EM algorithm}

\item{tol_EMstep}{tolerance for iterative estimation during EM step

Output:}

\item{return_covariances}{logical, indicates whether the state covariances
should be returned (might be necessary for evaluation
of the likelihood).}

\item{initialize}{either \code{"RRR"} or \code{"random"}, applies if no starting
values are handed over}

\item{Gamma_rrr}{type of normalization for the starting values obtained from RRR}
}
\value{
A named list of lists with elements
\itemize{
\item states: filtered (the filtered states)
smoothed (the smoothed states)
\item covariances: the filtered and smoothed covariances and lag-1 covariances
(if \code{return_covariances = TRUE})
\code{P_t^t} filtered covariances
\code{P_t^t-1}predicted covariances
\code{P_t^T} smoothed covariances
\verb{P_t-1t-2^T} smoothed lag-1 covariances
\item data: the data handed over to the algorithms
\code{X} predictors
\code{y} responses
\code{Z} transition matrices (X_t'beta (x) I_p)
\item parameters used during filtering: Sigma, Omega, beta
\item likelihoods: list containing Q and data loglikelihood for each iteration
\item convergence_information: Message containing parameter stages at convergence
}
}
\description{
Runs the expectation maximization algorithm for Model B
Input:
}
