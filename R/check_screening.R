#' Check the convergence of 1st moment with m at selected indices of a TMVN distribution
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
#' @param ind_test a vector of indices indexing the locs where we check the 1st-order moment convergence, by default, 10 random indices are used
#' @return a matrix summarizing the 1st moments evaluated at increasing m at the selected indices
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
#' first_mmt <- tmvn_check_converge(lb, ub, covmat,
#'   m_vec = seq(from = 10, to = 50, by = 10)
#' )
#' plot(first_mmt)
#'
#' # using locs, cov_name, and cov_parm
#' cov_name <- "matern15_isotropic"
#' cov_parm <- c(1, 0.1, 0.001)
#' first_mmt <- tmvn_check_converge(lb, ub,
#'   locs = locs, cov_name = cov_name, cov_parm = cov_parm,
#'   m_vec = seq(from = 10, to = 50, by = 10)
#' )
#' plot(first_mmt) + theme(text = element_text(size = 14))
#'
tmvn_check_converge <- function(cens_lb, cens_ub, covmat = NULL,
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
  class(first_mmt) <- "nntmvn_1stmmt_pred"
  return(first_mmt)
}

#' Check the convergence of 1st moment with m at selected indices of a PTMVN distribution with zero mean
#' @param y responses before censoring, of length n
#' @param cens_lb `cens_lb` and `cens_ub` define the censoring region, of length n
#' @param cens_ub `cens_lb` and `cens_ub` define the censoring region, of length n
#' @param covmat n-by-n dense covariance matrix, either `covmat` or `locs`,
#' `cov_name`, and `cov_parm` need to be provided
#' @param locs location matrix n X d
#' @param cov_name covariance function name from the `GpGp` package
#' @param cov_parm parameters for the covariance function from the `GpGp` package
#' @param m_vec a vector of `m` values (int) to be tested
#' @param N the number of samples to generate for each test index and each `m` to evaluate 1st-order moment
#' @param ind_test a vector of indices indexing the locs where we check the 1st-order moment convergence, by default, 10 random indices are used. If some test loc is not censored, the function treats it as unobserved
#' @return a matrix summarizing the 1st moments evaluated at increasing m at the selected indices
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
#' covmat <- GpGp::matern15_isotropic(c(1, 0.1, 0.01), locs)
#' y <- as.vector(t(chol(covmat)) %*% rnorm(n))
#' check_obj <- ptmvn_check_converge(y, lb , ub, covmat,
#'                                   m_vec = seq(from = 10, to = 50, by = 10)
#' )
#' first_mmt <- check_obj$pred
#' plot(first_mmt)
#' pred_err <- check_obj$error
#' plot(pred_err)
#' @export
ptmvn_check_converge <- function(y, cens_lb, cens_ub, covmat = NULL,
                                 locs = NULL, cov_name = NULL, cov_parm = NULL,
                                 m_vec = seq(from = 10, to = 100, by = 10), N = 1000,
                                 ind_test = NULL) {
  mask_cens <- (y < cens_ub) & (y > cens_lb)
  ind_cens <- which(mask_cens)
  if (length(ind_cens) == 0) {
    stop("There is no censored responses")
  }

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
    ind_test <- sample(ind_cens, min(length(ind_cens), 10), replace = FALSE)
  }

  first_mmt <- matrix(nrow = length(ind_test), ncol = length(m_vec))
  pred_err <- matrix(nrow = length(ind_test), ncol = length(m_vec))
  for (i in 1:length(ind_test)) {
    for (j in 1:length(m_vec)) {
      ind <- ind_test[i]
      m <- m_vec[j]
      NN_i <- NN[ind, 1:(m + 1)]
      a_i <- cens_lb[NN_i]
      b_i <- cens_ub[NN_i]
      y_i <- y[NN_i]
      covmat_i <- covmat[NN_i, NN_i]
      mask_cens_i <- mask_cens[NN_i]
      ind_obs_i <- which(!mask_cens_i[2:(m + 1)]) + 1

      if (length(ind_obs_i) > 0) {
        schur <- solve(
          covmat_i[ind_obs_i, ind_obs_i, drop = F],
          covmat_i[ind_obs_i, -ind_obs_i, drop = F]
        )
        cond_mean <- as.vector(t(y_i[ind_obs_i]) %*% schur)
        cond_covmat <- covmat_i[-ind_obs_i, -ind_obs_i, drop = F] -
          covmat_i[-ind_obs_i, ind_obs_i, drop = F] %*% schur
        cond_covmat <- (cond_covmat + t(cond_covmat)) / 2
        a_tmp <- a_i[-ind_obs_i]
        b_tmp <- b_i[-ind_obs_i]
      } else {
        cond_mean <- numeric(m + 1)
        cond_covmat <- covmat_i
        a_tmp <- a_i
        b_tmp <- b_i
      }

      if (!mask_cens_i[1]) { # if the first entry is not censored
        a_tmp[1] <- -Inf
        b_tmp[1] <- Inf
      }

      tryCatch(
        {
          samp <- withTimeout(
            expr = {
              TruncatedNormal::rtmvnorm(
                N, cond_mean,
                cond_covmat, a_tmp, b_tmp
              )
            },
            timeout = 30, onTimeout = "error"
          )
          if (length(a_tmp) == 1) {
            samp <- matrix(samp, ncol = 1)
          }
        },
        TimeoutException = function(ex) {
          message("Timeout (30s) occurred at ind = ", ind, ", m = : ", m)
          samp <- matrix(data = NA, nrow = sample_sz, ncol = length(a_tmp))
        }
      )
      first_mmt[i, j] <- mean(samp[, 1])
      pred_err[i, j] <- first_mmt[i, j] - y[ind]
    }
  }
  rownames(first_mmt) <- paste0("Loc ", ind_test)
  colnames(first_mmt) <- paste0("m = ", m_vec)
  class(first_mmt) <- "nntmvn_1stmmt_pred"
  rownames(pred_err) <- paste0("Loc ", ind_test)
  colnames(pred_err) <- paste0("m = ", m_vec)
  class(pred_err) <- "nntmvn_1stmmt_error"
  return(list(pred = first_mmt, error = pred_err))
}

#' Plot function for the `nntmvn_1stmmt_pred` class
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#' @return a ggplot object of class "gg" and "ggplot"
#' @export
plot.nntmvn_1stmmt_pred <- function(x) {
  ind_test <- as.integer(sub("Loc ", "", rownames(x)))
  first_mmt_df <- as.data.frame(x[, ])
  first_mmt_df <- first_mmt_df %>%
    mutate(loc_ind = ind_test) %>%
    pivot_longer(cols = starts_with("m = "), names_to = "m", values_to = "expectation") %>%
    mutate(m = as.numeric(sub("m = ", "", m)))
  ggplot(first_mmt_df, mapping = aes(x = m, y = expectation, group = loc_ind)) +
    geom_line() + 
    geom_text(data = first_mmt_df %>% filter(m == min(m)), 
              mapping = aes(label = loc_ind, x = m, y = expectation), 
              position = position_jitter(width = 1, height = 0))
}

#' Plot function for the `nntmvn_1stmmt_error` class
#' @import tidyr
#' @import dplyr
#' @import ggplot2
#' @return a ggplot object of class "gg" and "ggplot"
#' @export
plot.nntmvn_1stmmt_error <- function(x) {
  ind_test <- as.integer(sub("Loc ", "", rownames(x)))
  first_mmt_df <- as.data.frame(x[, ])
  first_mmt_df <- first_mmt_df %>%
    mutate(loc_ind = ind_test) %>%
    pivot_longer(cols = starts_with("m = "), names_to = "m", values_to = "error") %>%
    mutate(m = as.numeric(sub("m = ", "", m)))
  ggplot(first_mmt_df, mapping = aes(x = m, y = error, group = loc_ind)) +
    geom_line() + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_text(data = first_mmt_df %>% filter(m == min(m)), 
              mapping = aes(label = loc_ind, x = m, y = error), 
              position = position_jitter(width = 1, height = 0))
}
