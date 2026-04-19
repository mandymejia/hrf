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


#' Fit GLM for a Single Voxel
#'
#' Runs OLS and extracts task-specific betas and their covariance.
#'
#' @param y_v Numeric vector. BOLD timeseries for one voxel (nT x 1).
#' @param X_full Numeric matrix. Full design matrix (nT x p).
#' @param n_task Integer. Number of task regressors (first n_task columns of X_full).
#' @return List with beta_task, sigma2, df, Cov_task.
#' @keywords internal
fit_voxel_glm <- function(y_v, X_full, n_task) {
  ols <- compute_ols(y_v, X_full)
  df <- length(y_v) - ncol(X_full)
  sigma2 <- compute_residual_variance(y_v, X_full, ols$beta_hat, df)

  beta_task <- ols$beta_hat[1:n_task]
  Cov_task <- sigma2 * ols$XtX_inv[1:n_task, 1:n_task, drop = FALSE]

  list(
    beta_task = beta_task,
    sigma2 = sigma2,
    df = df,
    Cov_task = Cov_task
  )
}


#' Compute Contrasts
#'
#' Computes contrast estimates, standard errors, t-statistics, and p-values
#' from task betas and their covariance matrix.
#'
#' @param beta_task Numeric vector. Task beta estimates (n_task x 1).
#' @param Cov_task Numeric matrix. Covariance of task betas (n_task x n_task).
#' @param df Integer. Degrees of freedom.
#' @param A Numeric matrix. Contrast matrix (n_contrasts x n_task).
#' @return List with est, SE, tstat, pval (each length n_contrasts).
#' @keywords internal
compute_contrasts <- function(beta_task, Cov_task, df, A) {
  est <- as.numeric(A %*% beta_task)
  var_contrasts <- diag(A %*% Cov_task %*% t(A))
  SE <- sqrt(var_contrasts)
  tstat <- est / SE
  pval <- 2 * pt(-abs(tstat), df = df)

  list(est = est, SE = SE, tstat = tstat, pval = pval)
}


#' Apply P-Value Adjustment
#'
#' Applies multiple comparison correction per contrast across voxels.
#'
#' @param pval_mat Numeric matrix. P-values (n_voxels x n_contrasts).
#' @param method Character. Adjustment method (default "BH").
#' @return Adjusted p-value matrix (same dims).
#' @keywords internal
apply_p_adjustment <- function(pval_mat, method = "BH") {
  pval_adj <- pval_mat
  for (k in seq_len(ncol(pval_mat))) {
    pval_adj[, k] <- p.adjust(pval_mat[, k], method = method)
  }
  pval_adj
}


#' Resolve HRF Map from Regularize Result
#'
#' Determines which HRF map to use based on whether regularize
#' produced subject-level results (seffects) or population only.
#'
#' @param regularize_result Output from regularize_allHRFs.
#' @param subject_idx Integer subject index (required for seffects mode).
#' @return data.frame with columns: voxel, a1, b1, c.
#' @keywords internal
resolve_hrf_map <- function(regularize_result, subject_idx = NULL) {
  has_seffects <- !is.null(regularize_result$subject_results)

  if (has_seffects) {
    if (is.null(subject_idx)) {
      stop("subject_idx is required when regularize was run with seffects=TRUE")
    }
    winning_id <- regularize_result$subject_results$winning_candidate_id[subject_idx]
    cm <- regularize_result$candidate_maps[[winning_id]]
    data.frame(voxel = cm$voxel, a1 = cm$a1, b1 = cm$b1, c = cm$c)
  } else {
    pa <- regularize_result$pop_avg
    data.frame(voxel = pa$voxel, a1 = pa$a1_snapped, b1 = pa$b1_snapped, c = pa$c_snapped)
  }
}


#' Fit GLM for a Group of Voxels Sharing One HRF
#'
#' All voxels in the group share the same design matrix. Fits OLS
#' and computes contrasts for each voxel.
#'
#' @param vox_idx Integer vector. Voxel indices in this group.
#' @param y Numeric matrix. BOLD data (nT x nV_total).
#' @param X_full Numeric matrix. Full design matrix (nT x p).
#' @param XtX_inv Numeric matrix. Precomputed (X'X)^-1 (p x p).
#' @param n_task Integer. Number of task regressors.
#' @param A Numeric matrix. Contrast matrix (n_contrasts x n_task).
#' @return List with beta, est, SE, tstat, pval matrices (rows = voxels in group).
#' @keywords internal
fit_hrf_group <- function(vox_idx, y, X_full, XtX_inv, n_task, A) {
  nT <- nrow(y)
  df <- nT - ncol(X_full)
  n_contrasts <- nrow(A)
  n_vox <- length(vox_idx)

  beta_out <- matrix(NA_real_, nrow = n_vox, ncol = n_task)
  est_out <- matrix(NA_real_, nrow = n_vox, ncol = n_contrasts)
  SE_out <- matrix(NA_real_, nrow = n_vox, ncol = n_contrasts)
  tstat_out <- matrix(NA_real_, nrow = n_vox, ncol = n_contrasts)
  pval_out <- matrix(NA_real_, nrow = n_vox, ncol = n_contrasts)

  for (j in seq_len(n_vox)) {
    v <- vox_idx[j]
    y_v <- y[, v]

    if (mean(abs(y_v)) < 1e-6 || var(y_v) < 1e-6) next

    beta_hat <- as.numeric(XtX_inv %*% crossprod(X_full, y_v))
    resid <- y_v - X_full %*% beta_hat
    sigma2 <- sum(resid^2) / df

    beta_task <- beta_hat[1:n_task]
    Cov_task <- sigma2 * XtX_inv[1:n_task, 1:n_task, drop = FALSE]
    contr <- compute_contrasts(beta_task, Cov_task, df, A)

    beta_out[j, ] <- beta_task
    est_out[j, ] <- contr$est
    SE_out[j, ] <- contr$SE
    tstat_out[j, ] <- contr$tstat
    pval_out[j, ] <- contr$pval
  }

  list(beta = beta_out, est = est_out, SE = SE_out,
       tstat = tstat_out, pval = pval_out, df = df)
}


#' Fit Best HRF GLM for a Single Subject
#'
#' Fits a voxel-wise GLM using personalized HRFs from regularize_allHRFs.
#' Each voxel gets its own design matrix based on the population or
#' subject-specific HRF map. Voxels sharing the same HRF are grouped
#' for efficiency.
#'
#' @param regularize_result Output from \code{regularize_allHRFs()}.
#' @param BOLD_file Character. Path to subject's CIFTI file.
#' @param EVs Event data for this subject.
#' @param nuisance_file Character. Path to nuisance regressor file (or NULL).
#' @param TR Numeric. Repetition time in seconds.
#' @param subject_idx Integer. Subject index (required when seffects=TRUE).
#' @param contrasts Numeric matrix. Contrast matrix A (n_contrasts x n_task).
#'   Default NULL uses identity (each beta tested individually).
#' @param brainstructures Character vector. Brain structures to model.
#' @param resamp_res Integer. Resampling resolution.
#' @param hpf Numeric. High-pass filter cutoff.
#' @param onsets,offsets Logical. Include onset/offset regressors.
#' @param p_adjust_method Character. P-value adjustment method.
#' @param verbose Integer. Verbosity level.
#'
#' @return A list with class \code{"bestHRF"} containing betas, contrasts,
#'   df, contrast_matrix, subject_idx, hrf_map.
#'
#' @export
fit_bestHRF <- function(regularize_result,
                        BOLD_file,
                        EVs,
                        nuisance_file = NULL,
                        TR,
                        subject_idx = NULL,
                        contrasts = NULL,
                        brainstructures = c("left", "right"),
                        resamp_res = 10000,
                        hpf = 0.01,
                        onsets = FALSE,
                        offsets = FALSE,
                        p_adjust_method = "BH",
                        verbose = 1) {

  # Resolve HRF map
  hrf_map <- resolve_hrf_map(regularize_result, subject_idx)
  if (verbose > 0) {
    mode <- if (!is.null(regularize_result$subject_results)) "subject" else "population"
    cat("fit_bestHRF:", mode, "mode,", nrow(hrf_map), "voxels\n")
  }

  # Load BOLD data
  if (verbose > 0) cat("Loading BOLD data...\n")
  bold_data <- load_bold_data(BOLD_file, brainstructures, resamp_res, smoothing = FALSE)
  y <- t(as.matrix(bold_data$BOLD_xii))  # nT x nV
  nT <- nrow(y)
  nV <- ncol(y)

  # Build nuisance block (shared across all voxels)
  nuisance_mat <- if (!is.null(nuisance_file)) {
    as.matrix(utils::read.table(nuisance_file, header = FALSE))
  } else {
    NULL
  }
  nuisance_block <- build_nuisance_block(nuisance_mat, nT, TR, hpf)

  # Determine n_task and set up contrasts
  first_td <- build_task_design(EVs, nT, TR, hrf_map$a1[1], hrf_map$b1[1], hrf_map$c[1], onsets, offsets)
  n_task <- ncol(first_td$design)
  task_names <- first_td$task_names
  A <- if (is.null(contrasts)) diag(n_task) else contrasts
  if (ncol(A) != n_task) stop("Contrast matrix has ", ncol(A), " columns but there are ", n_task, " task regressors")

  # Allocate output
  beta_mat <- matrix(NA_real_, nrow = nV, ncol = n_task)
  est_mat <- matrix(NA_real_, nrow = nV, ncol = nrow(A))
  SE_mat <- matrix(NA_real_, nrow = nV, ncol = nrow(A))
  tstat_mat <- matrix(NA_real_, nrow = nV, ncol = nrow(A))
  pval_mat <- matrix(NA_real_, nrow = nV, ncol = nrow(A))

  # Group voxels by shared HRF
  hrf_map$hrf_key <- paste(hrf_map$a1, hrf_map$b1, hrf_map$c, sep = "_")
  voxel_groups <- split(hrf_map$voxel, hrf_map$hrf_key)
  unique_hrfs <- unique(hrf_map[, c("a1", "b1", "c", "hrf_key")])

  if (verbose > 0) cat("Fitting", nrow(unique_hrfs), "unique HRF groups...\n")

  # Fit each group
  df_val <- NA_integer_
  for (g in seq_len(nrow(unique_hrfs))) {
    td <- build_task_design(EVs, nT, TR, unique_hrfs$a1[g], unique_hrfs$b1[g], unique_hrfs$c[g], onsets, offsets)
    fd <- build_full_design(td$design, nuisance_block, nT)
    XtX_inv <- solve(crossprod(fd$X_full))
    vox_idx <- voxel_groups[[unique_hrfs$hrf_key[g]]]

    grp <- fit_hrf_group(vox_idx, y, fd$X_full, XtX_inv, n_task, A)

    beta_mat[vox_idx, ] <- grp$beta
    est_mat[vox_idx, ] <- grp$est
    SE_mat[vox_idx, ] <- grp$SE
    tstat_mat[vox_idx, ] <- grp$tstat
    pval_mat[vox_idx, ] <- grp$pval
    df_val <- grp$df
  }

  # P-value adjustment
  valid <- which(!is.na(pval_mat[, 1]))
  pval_adj_mat <- matrix(NA_real_, nrow = nV, ncol = nrow(A))
  if (length(valid) > 0) {
    pval_adj_mat[valid, ] <- apply_p_adjustment(pval_mat[valid, , drop = FALSE], p_adjust_method)
  }

  if (verbose > 0) cat("fit_bestHRF complete.\n")

  xii <- bold_data$BOLD_xii

  result <- list(
    betas = ciftiTools::newdata_xifti(xii, beta_mat),
    contrasts = list(
      est = ciftiTools::newdata_xifti(xii, est_mat),
      SE = ciftiTools::newdata_xifti(xii, SE_mat),
      tstat = ciftiTools::newdata_xifti(xii, tstat_mat),
      pval = ciftiTools::newdata_xifti(xii, pval_mat),
      pval_adj = ciftiTools::newdata_xifti(xii, pval_adj_mat)
    ),
    df = df_val,
    contrast_matrix = A,
    subject_idx = subject_idx,
    hrf_map = hrf_map
  )
  class(result) <- "bestHRF"
  result
}
