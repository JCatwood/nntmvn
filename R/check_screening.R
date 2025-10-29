#' Draw one sample of the underlying GP responses for a partially censored Gaussian
#' process using sequential nearest neighbor (SNN) method
#'
#' @importFrom TruncatedNormal rtmvnorm
#' @importFrom utils getFromNamespace
#' @importFrom RANN nn2
#' @importFrom R.utils withTimeout
#' @import GpGp
#' @import tidyr
#' @import dplyr
#' @param cens_lb lower bound vector for TMVN of length n
#' @param cens_ub upper bound vector for TMVN of length n
#' @param covmat n-by-n dense covariance matrix, either `covmat` or `locs`,
#' `cov_name`, and `cov_parm` need to be provided
#' @param locs location matrix n X d
#' @param cov_name covariance function name from the `GpGp` package
#' @param cov_parm parameters for the covariance function from the `GpGp` package
#' @param m_vec a vector of `m` values (int) to be tested
#' @param N the number of samples to generate for each test index and each `m` to evaluate 1st-order moment
#' @param ind_test a vector of indices indexing the locs where we check the 1st-order moment convergence, by default, 10 random indices are generated
#' @return a vector of length n representing the underlying GP responses
#' @export
#' @examples
#' library(GpGp)
#' library(nntmvn)
#' library(lhs)
#' library(ggplot2)
#' set.seed(123)
#' n <- 500
#' locs <- lhs::randomLHS(n, 2)
#' lb <- rep(-Inf, n)
#' ub <- rep(0, n)
#'
#' # using covariance matrix
#' covmat <- GpGp::matern15_isotropic(c(1, 0.1, 0.001), locs)
#' first_mmt <- check_1stmmt(lb, ub, covmat,
#'   m_vec = seq(from = 10, to = 30, by = 10)
#' )
#' plot(first_mmt)
#'
#' # using locs, cov_name, and cov_parm
#' cov_name <- "matern15_isotropic"
#' cov_parm <- c(1, 0.1, 0.001)
#' first_mmt <- check_1stmmt(lb, ub,
#'   locs = locs, cov_name = cov_name, cov_parm = cov_parm,
#'   m_vec = seq(from = 10, to = 30, by = 10)
#' )
#' plot(first_mmt) + theme(text = element_text(size = 14))
#'
check_1stmmt <- function(cens_lb, cens_ub, covmat = NULL,
                         locs = NULL, cov_name = NULL, cov_parm = NULL,
                         m_vec = seq(from = 10, to = 100, by = 10), N = 1000, ind_test = NULL) {
  if (any(cens_lb > cens_ub)) {
    stop("There exists entry in cens_lb greater than its counterpart in cens_ub\n")
  }
  n <- length(cens_lb)
  if (n <= max(m_vec)) {
    stop("The dimension of the input TMVN distribution should be greater than the
         maximum of m_vec ", max(m_vec))
  }
  if (!is.null(ind_test)) {
    if (min(ind_test) < 1 || max(ind_test) > n) {
      stop("The input ind_test is out of the range from 1 to n, n = ", n)
    }
  }

  if (is.null(covmat)) {
    if (any(c(is.null(locs), is.null(cov_name), is.null(cov_parm)))) {
      stop("locs, cov_name, and cov_parm must be provided when covmat is NULL")
    }
    covfunc <- getFromNamespace(cov_name, "GpGp")
    covmat <- covfunc(cov_parm, locs)
    NN <- RANN::nn2(locs, k = max(m_vec) + 1)[[1]]
  } else {
    NN <- corr_nn(covmat, max(m_vec))
  }

  if (is.null(ind_test)) {
    ind_test <- sample(1:n, min(n, 10), replace = FALSE)
  }

  first_mmt <- matrix(nrow = length(ind_test), ncol = length(m_vec))
  for (i in 1:length(ind_test)) {
    for (j in 1:length(m_vec)) {
      ind <- ind_test[i]
      m <- m_vec[j]
      NN_i <- NN[ind, 1:(m + 1)]
      a_i <- cens_lb[NN_i]
      b_i <- cens_ub[NN_i]
      covmat_i <- covmat[NN_i, NN_i]
      tryCatch(
        {
          samp <- withTimeout(
            expr = {
              TruncatedNormal::rtmvnorm(
                N, rep(0, m + 1),
                covmat_i, a_i, b_i
              )
            },
            timeout = 30, onTimeout = "error"
          )
        },
        TimeoutException = function(ex) {
          message("Timeout (30s) occurred at ind = ", ind, ", m = ", m)
          samp <- matrix(data = NA, nrow = sample_sz, ncol = length(a_i))
        }
      )
      first_mmt[i, j] <- mean(samp[, 1])
    }
  }
  rownames(first_mmt) <- paste0("Loc ", ind_test)
  colnames(first_mmt) <- paste0("m = ", m_vec)
  class(first_mmt) <- "nntmvn_1stmmt_check"
  return(first_mmt)
}

#' Plot function for the `nntmvn_1stmmt_check` class
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#' @return a ggplot object of class "gg" and "ggplot"
#' @export
plot.nntmvn_1stmmt_check <- function(x) {
  ind_test <- as.integer(sub("Loc ", "", rownames(x)))
  first_mmt_df <- as.data.frame(x[,])
  first_mmt_df <- first_mmt_df %>%
    mutate(loc_ind = ind_test) %>%
    pivot_longer(cols = starts_with("m = "), names_to = "m", values_to = "mean") %>%
    mutate(m = as.numeric(sub("m = ", "", m)))
  ggplot(first_mmt_df, mapping = aes(x = m, y = mean, group = loc_ind)) +
    geom_line()
}
