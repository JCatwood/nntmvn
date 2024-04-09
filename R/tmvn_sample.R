#' Simulate truncated multivariate normal (TMVN) as censored locations using
#'  nearest neighbors
#'
#' @importFrom TruncatedNormal rtmvnorm
#' @import GpGp
#' @param y uncensored responses of length n, where n is the number of locations
#' @param cens_lb lower bound vector for TMVN of length n
#' @param cens_ub upper bound vector for TMVN of length n
#' @param mask_cens mask for censored responses (also locations) of length n
#' @param NN n X m matrix for nearest neighbors. i-th row is the nearest neighbor indices of y_i. `NN[i, 1]` should be `i`
#' @param locs location matrix n X d
#' @param cov_name covariance function name from the `GpGp` package
#' @param cov_parm parameters for `cov_name`
#' @param covmat (optional) n-by-n dense covariance matrix
#' @return a vector of length n with imputed results for censored responses
#'
#' @export
seq_Vecc_samp_func <- function(y, cens_lb, cens_ub, mask_cens, NN, locs, cov_name,
                               cov_parm, covmat = NULL,
                               method = c("exact", "Vecchia"),
                               seed = NULL) {
  if (is.null(covmat)) {
    covfunc <- getFromNamespace(cov_name, "GpGp")
  } else {
    covfunc <- NULL
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  ind_cens <- which(mask_cens)
  method <- method[1]
  for (i in ind_cens) {
    NN_row <- NN[i, ]
    if (is.null(covfunc)) {
      covmat_sub <- covmat[NN_row, NN_row]
    } else {
      covmat_sub <- covfunc(cov_parm, locs[NN_row, , drop = FALSE])
    }
    mask_cens_sub <- mask_cens[NN_row]
    y_sub <- y[NN_row]
    cens_lb_sub <- cens_lb[NN_row]
    cens_ub_sub <- cens_ub[NN_row]
    n_cens_sub <- sum(mask_cens_sub)
    if (method == "exact") {
      if (n_cens_sub == length(NN_row)) {
        cond_covmat_sub_cens <- covmat_sub
        cond_mean_sub_cens <- rep(0, n_cens_sub)
      } else {
        tmp_mat <- covmat_sub[mask_cens_sub, !mask_cens_sub] %*%
          solve(covmat_sub[!mask_cens_sub, !mask_cens_sub])
        cond_covmat_sub_cens <- covmat_sub[mask_cens_sub, mask_cens_sub] -
          tmp_mat %*% covmat_sub[!mask_cens_sub, mask_cens_sub]
        cond_mean_sub_cens <- as.vector(tmp_mat %*% y_sub[!mask_cens_sub])
      }
      samp_cens_sub <- t(TruncatedNormal::rtmvnorm(
        1, cond_mean_sub_cens,
        cond_covmat_sub_cens,
        cens_lb_sub[mask_cens_sub], cens_ub_sub[mask_cens_sub]
      ))
    } else if (method == "Vecchia") {
      stop("Vecchia method is not implemented")
    } else {
      stop("Undefined method\n")
    }
    y[i] <- samp_cens_sub[1]
    mask_cens[i] <- FALSE
  }
  y
}
