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
#' @return Numeric matrix: nuisance columns + DCT bases (or just DCT if nuisance is NULL).
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
#' Assembles the complete design matrix: task + intercept + nuisance_block.
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
  pval <- 2 * stats::pt(-abs(tstat), df = df)

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
    pval_adj[, k] <- stats::p.adjust(pval_mat[, k], method = method)
  }
  pval_adj
}


#' Validate fit_bestHRF Inputs
#'
#' Checks consistency between `use`, the regularize result, and subject_idx.
#'
#' @param regularize_result Output from regularize_allHRFs.
#' @param use Either "personalized" or "adapted".
#' @param subject_idx Subject index (required when use = "personalized").
#' @keywords internal
validate_bestHRF_inputs <- function(regularize_result, use, subject_idx) {
  has_seffects <- !is.null(regularize_result$subject_results)

  if ("personalized" %in% use) {
    if (!has_seffects) {
      stop("use including 'personalized' requires regularize_allHRFs to have been run with seffects = TRUE. ",
           "Re-run with seffects = TRUE, or call fit_bestHRF with use = 'adapted' only.")
    }
    if (is.null(subject_idx)) {
      stop("subject_idx is required when use includes 'personalized'.")
    }
    n_subjects <- nrow(regularize_result$subject_results)
    if (subject_idx < 1 || subject_idx > n_subjects) {
      stop("subject_idx = ", subject_idx, " is out of range (1 to ", n_subjects, ").")
    }
  }
}


#' Resolve HRF Map from Regularize Result
#'
#' Returns the HRF map for a single mode: the subject's personalized map
#' (winning candidate) or the population-average ("adapted") map.
#'
#' @param regularize_result Output from regularize_allHRFs.
#' @param mode One of "personalized" or "adapted".
#' @param subject_idx Integer subject index (required when mode = "personalized").
#' @return data.frame with columns: voxel, a1, b1, c.
#' @keywords internal
resolve_hrf_map <- function(regularize_result, mode, subject_idx = NULL) {
  validate_bestHRF_inputs(regularize_result, mode, subject_idx)

  if (mode == "personalized") {
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


#' Fit All Voxels Given an HRF Map
#'
#' Iterates over unique HRFs in `hrf_map_df`, building a design matrix per
#' unique HRF and fitting all voxels that share it. Returns matrices of
#' betas, contrast estimates, SEs, t-stats, raw + adjusted p-values.
#'
#' @param hrf_map_df data.frame with columns voxel, a1, b1, c.
#' @param y Numeric matrix. BOLD data (nT x nV).
#' @param EVs Event data for one subject.
#' @param nT Integer. Number of timepoints.
#' @param TR Numeric. Repetition time.
#' @param nuisance_block Numeric matrix. Nuisance + DCT (nT x n_nuisance).
#' @param A Numeric matrix. Contrast matrix (n_contrasts x n_task).
#' @param n_task Integer. Number of task regressors.
#' @param nV Integer. Total number of voxels.
#' @param label Character. Used in verbose log line.
#' @param onsets,offsets Logical. Include onset/offset regressors.
#' @param p_adjust_method Character. P-value adjustment method.
#' @param verbose Integer. Verbosity level.
#' @return List of matrices: beta_mat, est_mat, SE_mat, tstat_mat, pval_mat,
#'   pval_adj_mat (each nV x n_task or nV x n_contrasts), and df scalar.
#' @keywords internal
fit_all_voxels <- function(hrf_map_df, y, EVs, nT, TR, nuisance_block,
                           A, n_task, nV, label,
                           onsets = FALSE, offsets = FALSE,
                           p_adjust_method = "BH", verbose = 1) {
  beta_mat  <- matrix(NA_real_, nrow = nV, ncol = n_task)
  est_mat   <- matrix(NA_real_, nrow = nV, ncol = nrow(A))
  SE_mat    <- matrix(NA_real_, nrow = nV, ncol = nrow(A))
  tstat_mat <- matrix(NA_real_, nrow = nV, ncol = nrow(A))
  pval_mat  <- matrix(NA_real_, nrow = nV, ncol = nrow(A))

  hrf_map_df$hrf_key <- paste(hrf_map_df$a1, hrf_map_df$b1, hrf_map_df$c, sep = "_")
  voxel_groups <- split(hrf_map_df$voxel, hrf_map_df$hrf_key)
  unique_hrfs  <- unique(hrf_map_df[, c("a1", "b1", "c", "hrf_key")])

  if (verbose > 0) cat("Fitting", label, "-", nrow(unique_hrfs), "unique HRF groups...\n")

  # df is invariant across HRF groups (same nT, same X_full ncol). Set once.
  df_val <- NA_integer_
  for (g in seq_len(nrow(unique_hrfs))) {
    td <- build_task_design(EVs, nT, TR, unique_hrfs$a1[g], unique_hrfs$b1[g], unique_hrfs$c[g], onsets, offsets)
    fd <- build_full_design(td$design, nuisance_block, nT)
    XtX_inv <- solve(crossprod(fd$X_full))
    vox_idx <- voxel_groups[[unique_hrfs$hrf_key[g]]]

    grp <- fit_hrf_group(vox_idx, y, fd$X_full, XtX_inv, n_task, A)

    beta_mat[vox_idx, ]  <- grp$beta
    est_mat[vox_idx, ]   <- grp$est
    SE_mat[vox_idx, ]    <- grp$SE
    tstat_mat[vox_idx, ] <- grp$tstat
    pval_mat[vox_idx, ]  <- grp$pval

    if (is.na(df_val)) df_val <- grp$df
    else if (df_val != grp$df) stop("df differs across HRF groups: ", df_val, " vs ", grp$df)
  }

  valid <- which(!is.na(pval_mat[, 1]))
  pval_adj_mat <- matrix(NA_real_, nrow = nV, ncol = nrow(A))
  if (length(valid) > 0) {
    pval_adj_mat[valid, ] <- apply_p_adjustment(pval_mat[valid, , drop = FALSE], p_adjust_method)
  }

  list(beta_mat = beta_mat, est_mat = est_mat, SE_mat = SE_mat,
       tstat_mat = tstat_mat, pval_mat = pval_mat, pval_adj_mat = pval_adj_mat,
       df = df_val)
}


#' Package Raw Fit Matrices into a List of xifti Objects
#'
#' Wraps each numeric matrix produced by \code{\link{fit_all_voxels}} into an
#' xifti object using \code{xii} as the template.
#'
#' @param raw List from \code{fit_all_voxels()}: must contain beta_mat, est_mat,
#'   SE_mat, tstat_mat, pval_mat, pval_adj_mat.
#' @param xii Template xifti object.
#' @return List with \code{$betas} and \code{$contrasts$\{est,SE,tstat,pval,pval_adj\}}.
#' @keywords internal
package_results <- function(raw, xii) {
  list(
    betas = ciftiTools::newdata_xifti(xii, raw$beta_mat),
    contrasts = list(
      est      = ciftiTools::newdata_xifti(xii, raw$est_mat),
      SE       = ciftiTools::newdata_xifti(xii, raw$SE_mat),
      tstat    = ciftiTools::newdata_xifti(xii, raw$tstat_mat),
      pval     = ciftiTools::newdata_xifti(xii, raw$pval_mat),
      pval_adj = ciftiTools::newdata_xifti(xii, raw$pval_adj_mat)
    )
  )
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
#' @param use Character vector. Subset of \code{c("personalized", "adapted")};
#'   pass either (or both) to request the corresponding HRF fit alongside the
#'   always-present working fit. Defaults to "personalized". Including
#'   "personalized" requires `regularize_result` to have been built with
#'   `seffects = TRUE` and a non-NULL `subject_idx`.
#' @param subject_idx Integer. Subject index (required when use includes
#'   "personalized").
#' @param contrasts Numeric matrix. Contrast matrix A (n_contrasts x n_task).
#'   Default NULL uses identity (each beta tested individually).
#' @param working_hrf List with elements \code{a1}, \code{b1}, \code{c}
#'   parameterizing the canonical "working" HRF that produces the
#'   \code{$working} section of the result. Defaults to the SPM canonical
#'   shape (\code{a1 = 6, b1 = 1, c = 1/6}).
#' @param brainstructures Character vector. Brain structures to model.
#' @param resamp_res Integer. Resampling resolution.
#' @param hpf Numeric. High-pass filter cutoff.
#' @param onsets,offsets Logical. Include onset/offset regressors.
#' @param p_adjust_method Character. P-value adjustment method.
#' @param verbose Integer. Verbosity level.
#'
#' @return A list with class \code{"bestHRF"} containing:
#'   \describe{
#'     \item{working}{Canonical HRF fit: betas and contrasts (as xifti). Always
#'       present.}
#'     \item{personalized}{Subject-specific best HRF fit (when use is
#'       "personalized" or "both").}
#'     \item{adapted}{Population-average HRF fit (when use is "adapted" or
#'       "both").}
#'     \item{df}{Degrees of freedom}
#'     \item{contrast_matrix}{The contrast matrix A used}
#'   }
#'
#' @export
fit_bestHRF <- function(regularize_result,
                        BOLD_file,
                        EVs,
                        nuisance_file = NULL,
                        TR,
                        use = c("personalized", "adapted"),
                        subject_idx = NULL,
                        contrasts = NULL,
                        working_hrf = list(a1 = 6, b1 = 1, c = 1/6),
                        brainstructures = c("left", "right"),
                        resamp_res = 10000,
                        hpf = 0.01,
                        onsets = FALSE,
                        offsets = FALSE,
                        p_adjust_method = "BH",
                        verbose = 1) {

  # Vector of one or both: c("personalized") | c("adapted") | c("personalized","adapted").
  use <- match.arg(use, several.ok = TRUE)
  modes <- use
  hrf_maps <- setNames(
    lapply(modes, function(m) resolve_hrf_map(regularize_result, m, subject_idx)),
    modes
  )
  if (verbose > 0) cat("fit_bestHRF: modes =", paste(modes, collapse = ", "),
                       "(", nrow(hrf_maps[[1]]), "voxels )\n")

  # Load BOLD data
  if (verbose > 0) cat("Loading BOLD data...\n")
  bold_data <- load_bold_data(BOLD_file, brainstructures, resamp_res, smoothing = TRUE, scale = TRUE)
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

  # Determine n_task and set up contrasts. Use working_hrf params (always
  # valid SPM canonical) since task structure depends on EVs, not HRF shape.
  first_td <- build_task_design(EVs, nT, TR,
                                working_hrf$a1, working_hrf$b1, working_hrf$c,
                                onsets, offsets)
  n_task <- ncol(first_td$design)
  task_names <- first_td$task_names
  A <- if (is.null(contrasts)) diag(n_task) else contrasts
  if (ncol(A) != n_task) stop("Contrast matrix has ", ncol(A), " columns but there are ", n_task, " task regressors")

  xii <- bold_data$BOLD_xii

  # Fit working HRF (same HRF for all voxels)
  working_map <- data.frame(voxel = 1:nV, a1 = working_hrf$a1, b1 = working_hrf$b1, c = working_hrf$c)
  raw_working <- fit_all_voxels(working_map, y, EVs, nT, TR, nuisance_block,
                                A, n_task, nV, "working HRF",
                                onsets, offsets, p_adjust_method, verbose)

  # Fit each requested mode (personalized and/or adapted)
  mode_sections <- list()
  for (m in modes) {
    raw_m <- fit_all_voxels(hrf_maps[[m]], y, EVs, nT, TR, nuisance_block,
                            A, n_task, nV, m,
                            onsets, offsets, p_adjust_method, verbose)
    sec <- package_results(raw_m, xii)
    sec$hrf_assignments <- hrf_maps[[m]]
    if (m == "personalized") sec$subject_idx <- subject_idx
    mode_sections[[m]] <- sec
  }

  if (verbose > 0) cat("fit_bestHRF complete.\n")

  working_section <- package_results(raw_working, xii)
  working_section$hrf_assignments <- working_map

  result <- c(list(working = working_section), mode_sections,
              list(df = raw_working$df, contrast_matrix = A))

  class(result) <- "bestHRF"
  result
}
