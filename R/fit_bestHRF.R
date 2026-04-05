


#' Compute OLS Estimates
#'
#' @param y Numeric vector. Response (nT x 1).
#' @param X Numeric matrix. Design matrix (nT x p).
#' @return List with beta_hat (p x 1) and XtX_inv (p x p).
#' @keywords internal
compute_ols <- function(y, X) {
  XtX_inv <- solve(crossprod(X))
  beta_hat <- XtX_inv %*% crossprod(X, y)
  list(beta_hat = as.numeric(beta_hat), XtX_inv = XtX_inv)
}


#' Compute Residual Variance
#'
#' @param y Numeric vector. Response (nT x 1).
#' @param X Numeric matrix. Design matrix (nT x p).
#' @param beta_hat Numeric vector. OLS estimates (p x 1).
#' @param df Integer. Degrees of freedom (nT - p).
#' @return Scalar residual variance estimate (sigma_hat^2).
#' @keywords internal
compute_residual_variance <- function(y, X, beta_hat, df) {
  resid <- y - X %*% beta_hat
  sum(resid^2) / df
}


#' Build Task Design Matrix
#'
#' Wraps make_design() to build the task portion of the design matrix
#' for a specific HRF parameterization.
#'
#' @param EVs Event data for one subject.
#' @param nT Integer. Number of timepoints.
#' @param TR Numeric. Repetition time.
#' @param a1,b1,c HRF parameters.
#' @param onsets,offsets Logical. Include onset/offset regressors.
#' @return List with design (nT x n_task matrix) and task_names.
#' @keywords internal
build_task_design <- function(EVs, nT, TR, a1, b1, c, onsets, offsets) {
  a2 <- (16 / sqrt(6)) * sqrt(a1) * sqrt(b1)
  b2 <- b1

  taper_start <- get_taper_start(
    a1 = a1, b1 = b1, a2 = a2, b2 = b2,
    c = c, TR = TR, deriv = 0
  )

  result <- make_design(
    EVs = EVs, nTime = nT, TR = TR, dHRF = 0,
    onset = onsets, offset = offsets,
    taper_start = taper_start,
    a1 = a1, b1 = b1, c = c, a2 = a2, b2 = b2
  )

  list(
    design = result$design,
    task_names = colnames(result$design)
  )
}


#' Build Nuisance Block
#'
#' Combines user-provided nuisance regressors with DCT bases for
#' high-pass filtering. Called once per subject (shared across voxels).
#'
#' @param nuisance_mat Numeric matrix of nuisance regressors (or NULL).
#' @param nT Integer. Number of timepoints.
#' @param TR Numeric. Repetition time.
#' @param hpf Numeric. High-pass filter cutoff frequency.
#' @return Numeric matrix: [nuisance | DCT_bases] (or just DCT if nuisance is NULL).
#' @keywords internal
build_nuisance_block <- function(nuisance_mat, nT, TR, hpf) {
  n_dct <- round(fMRItools::dct_convert(T_ = nT, TR = TR, f = hpf))
  DCT <- fMRItools::dct_bases(nT, n_dct)

  if (!is.null(nuisance_mat)) {
    cbind(nuisance_mat, DCT)
  } else {
    DCT
  }
}


#' Build Full Design Matrix
#'
#' Assembles the complete design matrix: [task | intercept | nuisance_block].
#' Matches the column order used in GLM_multi.R.
#'
#' @param task_design Matrix. Task regressors (nT x n_task).
#' @param nuisance_block Matrix. Nuisance + DCT (nT x n_nuisance).
#' @param nT Integer. Number of timepoints.
#' @return List with X_full (nT x p) and n_task.
#' @keywords internal
build_full_design <- function(task_design, nuisance_block, nT) {
  X_full <- cbind(task_design, rep(1, nT), nuisance_block)
  list(
    X_full = X_full,
    n_task = ncol(task_design)
  )
}
