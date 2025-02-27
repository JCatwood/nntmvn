% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmvn_sample.R
\name{rtmvn}
\alias{rtmvn}
\title{Draw one sample from a truncated multivariate normal (TMVN) distribution
using sequential nearest neighbor (SNN) method}
\usage{
rtmvn(
  cens_lb,
  cens_ub,
  m = 30,
  covmat = NULL,
  locs = NULL,
  cov_name = NULL,
  cov_parm = NULL,
  NN = NULL,
  ordering = 0,
  seed = NULL
)
}
\arguments{
\item{cens_lb}{lower bound vector for TMVN of length n}

\item{cens_ub}{upper bound vector for TMVN of length n}

\item{m}{positive integer for the number of nearest neighbors used}

\item{covmat}{n-by-n dense covariance matrix, either \code{covmat} or \code{locs},
\code{cov_name}, and \code{cov_parms} need to be provided}

\item{locs}{location matrix n X d}

\item{cov_name}{covariance function name from the \code{GpGp} package}

\item{cov_parm}{parameters for the covariance function from the \code{GpGp} package}

\item{NN}{n X m matrix for nearest neighbors. i-th row is the nearest neighbor indices of y_i. \code{NN[i, 1]} should be \code{i}}

\item{ordering}{\code{0} for do not reorder, \code{1} for variance descending order, \code{2} for maximin ordering}

\item{seed}{set seed for reproducibility}
}
\value{
a vector of length n representing the underlying GP responses
}
\description{
Draw one sample from a truncated multivariate normal (TMVN) distribution
using sequential nearest neighbor (SNN) method
}
\examples{
library(nntmvn)
library(TruncatedNormal)
set.seed(123)
x <- matrix(seq(from = 0, to = 1, length.out = 51), ncol = 1)
cov_name <- "matern15_isotropic"
cov_parm <- c(1.0, 0.1, 0.001) #'' variance, range, nugget
cov_func <- getFromNamespace(cov_name, "GpGp")
covmat <- cov_func(cov_parm, x)
lb <- rep(-Inf, nrow(x))
ub <- rep(-1, nrow(x))
m <- 30
samp_SNN <- matrix(NA, 3, nrow(x))
for (i in 1:3) {
  samp_SNN[i, ] <- nntmvn::rtmvn(lb, ub, m = m, covmat = covmat, locs = x, ordering = 0)
}
samp_TN <- TruncatedNormal::rtmvnorm(3, rep(0, nrow(x)), covmat, lb, ub)
qqplot(samp_SNN, samp_TN, xlim = range(samp_SNN, samp_TN), ylim = range(samp_SNN, samp_TN))
abline(a = 0, b = 1, lty = "dashed", col = "red")

}
