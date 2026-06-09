#' Build the Population HRF Template from fit_allHRFs Results
#'
#' Aggregates per-subject best-HRF parameters from \code{fit_allHRFs()} into a
#' population template: per-voxel mean time-to-peak / FWHM, unmasked + smoothed,
#' snapped back to the HRF grid. Per-subject candidate scoring lives in
#' \code{\link{fit_bestHRF}} now -- this function is aggregation only.
#'
#' @param workingHRF_results Results from \code{fit_workingHRF()}. Needed for
#'   the activation masks (which voxels are kept in pop_avg) and the xifti
#'   geometry template attached for downstream plotting.
#' @param allHRF_results Results from \code{fit_allHRFs()}.
#' @param verbose Integer verbosity level.
#'
#' @return A list with class \code{"regularizeHRFs"} containing:
#'   \describe{
#'     \item{pop_avg}{Population average HRF map per voxel.}
#'     \item{best_params_df}{Per-subject best HRF parameters (intermediate).}
#'     \item{winning_c}{Winning c value (scalar).}
#'     \item{c_votes}{Vote counts per c value.}
#'     \item{hrf_grid}{HRF grid with t2p/fwhm columns.}
#'     \item{mask_prop_NA}{Population activation mask.}
#'   }
#'
#' @export
regularize_allHRFs <- function(workingHRF_results,
                                allHRF_results,
                                verbose = 1) {

  # Step 0: Precompute metrics on HRF grid
  hrf_grid <- allHRF_results$hrf_grid
  if (verbose > 0) cat("\nPrecomputing t2p/fwhm for HRF grid...\n")
  metrics <- get_hrf_metrics(hrf_grid$a1, hrf_grid$b1, hrf_grid$c)
  hrf_grid$time_to_peak <- metrics$time_to_peak
  hrf_grid$FWHM <- metrics$FWHM

  # Step 1 (pre): Determine winning c
  if (verbose > 0) cat("\nDetermining winning c value...\n")
  winning_result <- determine_winning_c(allHRF_results, workingHRF_results, verbose = verbose)

  # Step 1-2: Extract best params per subject with t2p/fwhm
  if (verbose > 0) cat("\nExtracting best params and converting to t2p/fwhm...\n")
  best_params_df <- extract_best_params_per_subject(
    allHRF_results, workingHRF_results,
    winning_result$winning_c, hrf_grid, verbose = verbose
  )

  # Step 3: Average t2p/fwhm across subjects
  if (verbose > 0) cat("\nComputing population average...\n")
  mask_prop_NA <- workingHRF_results[["activation_masks"]][["mask_prop_NA"]]
  pop_mask_voxels <- which(!is.na(mask_prop_NA))

  pop_avg <- stats::aggregate(
    cbind(time_to_peak, FWHM) ~ voxel,
    data = best_params_df[best_params_df$mask & best_params_df$voxel %in% pop_mask_voxels, ],
    FUN = mean
  )
  names(pop_avg)[2:3] <- c("t2p_mean", "fwhm_mean")

  if (verbose > 0) cat("Population average computed for", nrow(pop_avg), "voxels\n")

  # Step 3.5: Unmask (fill + smooth) t2p_mean / fwhm_mean before snapping
  pop_avg <- unmask_pop_avg(pop_avg, workingHRF_results, mask_prop_NA, verbose = verbose)

  # Step 4: Snap back to a1/b1 grid
  if (verbose > 0) cat("\nSnapping t2p/fwhm back to grid...\n")
  snapped <- snap_to_grid_t2p_fwhm(pop_avg$t2p_mean, pop_avg$fwhm_mean, hrf_grid, winning_result$winning_c)
  pop_avg$a1_snapped <- snapped$a1
  pop_avg$b1_snapped <- snapped$b1
  pop_avg$c_snapped <- snapped$c

  result <- list(
    pop_avg        = pop_avg,
    best_params_df = best_params_df,
    winning_c      = winning_result$winning_c,
    c_votes        = winning_result$c_votes,
    hrf_grid       = hrf_grid,
    mask_prop_NA   = mask_prop_NA
  )

  # Attach a xifti geometry template (from any successful workingHRF subject)
  # so plot.regularizeHRFs methods don't need the user to pass one.
  attr(result, "xii") <- extract_xii_template(workingHRF_results)

  class(result) <- "regularizeHRFs"

  if (verbose > 0) cat("\nregularize_allHRFs complete.\n")
  return(result)
}


#' Compute HRF Metrics (Time-to-Peak and FWHM)
#'
#' Converts HRF parameters (a1, b1, c) to interpretable metrics:
#' time-to-peak and full-width-at-half-maximum.
#'
#' @param a1 Numeric vector of a1 parameter values.
#' @param b1 Numeric vector of b1 parameter values.
#' @param c Numeric vector of c parameter values.
#' @param tapered Logical. If TRUE (default), applies tapering for c > 0
#'   to handle the undershoot properly when computing FWHM.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{a1, b1, c}{Input parameter values}
#'     \item{time_to_peak}{Time to peak (= a1 - b1)}
#'     \item{FWHM}{Full-width-at-half-maximum of the HRF}
#'   }
#'
#' @keywords internal
get_hrf_metrics <- function(a1, b1, c, tapered = TRUE) {
  n <- length(a1)
  time_to_peak <- numeric(n)
  FWHM <- numeric(n)

  # Time grid: TR=2, sr_factor=100, max_time=20
  inds <- seq(1/100, 20, 1/100) * 2

  for (i in 1:n) {
    a1_i <- a1[i]
    b1_i <- b1[i]
    c_i <- c[i]

    # Compute a2 from a1, b1 (standard relationship)
    a2_i <- (16 / sqrt(6)) * sqrt(a1_i) * sqrt(b1_i)
    b2_i <- b1_i

    # Time-to-peak is simply a1 - b1
    time_to_peak[i] <- a1_i - b1_i

    # Compute HRF curve
    hrf_raw <- HRF_calc(t = inds, deriv = 0, a1 = a1_i, b1 = b1_i,
                        a2 = a2_i, b2 = b2_i, c = c_i)

    # Apply tapering if needed (for c > 0, the undershoot can cause issues)
    if (tapered && c_i > 0) {
      peak2_idx <- which.min(hrf_raw)
      peak2_time <- inds[peak2_idx]
      taper_start <- min(peak2_time, 25)

      # Only taper if HRF hasn't returned to baseline by t=30
      if (abs(hrf_raw[which.min(abs(inds - 30))]) > 0.01) {
        hrf_use <- HRF_calc(t = inds, deriv = 0, a1 = a1_i, b1 = b1_i,
                            a2 = a2_i, b2 = b2_i, c = c_i,
                            taper_start = taper_start, taper_end = 30, taper_power = 1)
      } else {
        hrf_use <- hrf_raw
      }
    } else {
      hrf_use <- hrf_raw
    }

    # Compute FWHM from the HRF curve
    peak_val <- max(hrf_use)
    peak_idx <- which.max(hrf_use)
    vals_left <- hrf_use[1:peak_idx]
    vals_right <- hrf_use[peak_idx:length(hrf_use)]

    x1 <- min(which(vals_left > 0.5 * peak_val))
    x2 <- max(which(vals_right > 0.5 * peak_val))
    time1 <- inds[x1]
    time2 <- inds[peak_idx + x2 - 1]
    FWHM[i] <- time2 - time1
  }

  return(data.frame(a1 = a1, b1 = b1, c = c, time_to_peak = time_to_peak, FWHM = FWHM))
}


#' Determine Winning C Value Across Subjects
#'
#' Compares RSS between c values (e.g., c=0 vs c=1/6) across all subjects
#' and voxels to determine which c value provides better fits overall.
#'
#' @param allHRF_results Results from \code{fit_allHRFs()}.
#' @param workingHRF_results Results from \code{fit_workingHRF()}.
#' @param verbose Integer verbosity level.
#'
#' @return A list with:
#'   \describe{
#'     \item{winning_c}{The c value with more wins (lower RSS more often)}
#'     \item{c_votes}{Named vector with vote counts per c value}
#'     \item{c_keys}{Character vector of c value keys}
#'   }
#'
#' @keywords internal
determine_winning_c <- function(allHRF_results, workingHRF_results, verbose = 1) {
  # Get masks
  mask_prop_NA <- workingHRF_results[["activation_masks"]][["mask_prop_NA"]]
  pop_mask_voxels <- which(!is.na(mask_prop_NA))
  subject_masks <- workingHRF_results[["activation_masks"]][["masks"]]

  n_subjects <- length(allHRF_results[["subject_results"]])
  n_voxels <- nrow(subject_masks)

  # Get c value keys from first subject
  c_keys <- names(allHRF_results[["subject_results"]][[1]][["glm_result"]][["mGLM0s"]][["cortexL"]][["best_per_c"]])

  if (verbose > 0) cat("C values found:", c_keys, "\n")

  # Initialize vote counts
  c_votes <- setNames(rep(0, length(c_keys)), c_keys)

  # Count votes across all subjects
  for (i in 1:n_subjects) {
    subj_result <- allHRF_results[["subject_results"]][[i]][["glm_result"]][["mGLM0s"]]
    subj_mask <- subject_masks[, i]

    # Get RSS for each c value (combine cortexL and cortexR)
    rss_by_c <- lapply(c_keys, function(k) {
      c(subj_result[["cortexL"]][["best_per_c"]][[k]]$rss,
        subj_result[["cortexR"]][["best_per_c"]][[k]]$rss)
    })
    names(rss_by_c) <- c_keys

    # Apply masks: both subject mask and population mask
    both_mask <- subj_mask & (seq_len(n_voxels) %in% pop_mask_voxels)

    # For masked voxels, find which c has lower RSS and count votes
    for (k in c_keys) {
      rss_k <- rss_by_c[[k]][both_mask]
      # Compare against other c values
      other_keys <- setdiff(c_keys, k)
      wins <- rep(TRUE, sum(both_mask))
      for (ok in other_keys) {
        wins <- wins & (rss_k < rss_by_c[[ok]][both_mask])
      }
      c_votes[k] <- c_votes[k] + sum(wins)
    }
  }

  # Determine winning c
  winning_c_key <- names(which.max(c_votes))
  winning_c <- as.numeric(winning_c_key)

  if (verbose > 0) {
    cat("C votes:", paste(names(c_votes), "=", c_votes, collapse = ", "), "\n")
    cat("Winning c:", winning_c, "\n")
  }

  return(list(
    winning_c = winning_c,
    c_votes = c_votes,
    c_keys = c_keys
  ))
}


#' Extract Best HRF Parameters Per Subject
#'
#' Extracts the best HRF parameters for each voxel and subject using the
#' winning c value. Combines cortexL and cortexR data and adds t2p/fwhm
#' metrics via lookup from the precomputed hrf_grid.
#'
#' @param allHRF_results Results from \code{fit_allHRFs()}.
#' @param workingHRF_results Results from \code{fit_workingHRF()}.
#' @param winning_c The winning c value from \code{determine_winning_c()}.
#' @param hrf_grid HRF grid with precomputed time_to_peak and FWHM columns.
#' @param verbose Integer verbosity level.
#'
#' @return A data.frame with columns:
#'   \describe{
#'     \item{voxel}{Voxel index (1 to n_voxels)}
#'     \item{subject}{Subject index}
#'     \item{a1, b1, c}{Best HRF parameters for winning c}
#'     \item{rss}{Residual sum of squares}
#'     \item{time_to_peak, FWHM}{HRF metrics looked up from grid}
#'     \item{mask}{Subject activation mask (TRUE/FALSE)}
#'   }
#'
#' @keywords internal
extract_best_params_per_subject <- function(allHRF_results, workingHRF_results,
                                             winning_c, hrf_grid, verbose = 1) {
  # Get subject masks

subject_masks <- workingHRF_results[["activation_masks"]][["masks"]]
  n_subjects <- length(allHRF_results[["subject_results"]])
  n_voxels <- nrow(subject_masks)

  # Convert winning_c to character key
  winning_c_key <- as.character(winning_c)

  # Create lookup for t2p/fwhm from hrf_grid
  # Key: "a1_b1_c" -> row index in hrf_grid
  hrf_grid$key <- paste(hrf_grid$a1, hrf_grid$b1, hrf_grid$c, sep = "_")

  if (verbose > 0) cat("Extracting best params for", n_subjects, "subjects...\n")

  best_params_list <- vector("list", n_subjects)

  for (i in seq_len(n_subjects)) {
    if (verbose > 1 && i %% 100 == 0) cat("  Subject", i, "/", n_subjects, "\n")

    subj_result <- allHRF_results[["subject_results"]][[i]][["glm_result"]][["mGLM0s"]]
    subj_mask <- subject_masks[, i]

    # Get best params for winning c (cortexL + cortexR)
    best_L <- subj_result[["cortexL"]][["best_per_c"]][[winning_c_key]]
    best_R <- subj_result[["cortexR"]][["best_per_c"]][[winning_c_key]]

    # Combine hemispheres
    a1_vec <- c(best_L$a1, best_R$a1)
    b1_vec <- c(best_L$b1, best_R$b1)
    c_vec <- c(best_L$c, best_R$c)
    rss_vec <- c(best_L$rss, best_R$rss)

    # Create lookup keys and get t2p/fwhm
    keys <- paste(a1_vec, b1_vec, c_vec, sep = "_")
    grid_idx <- match(keys, hrf_grid$key)

    best_params_list[[i]] <- data.frame(
      voxel = seq_len(n_voxels),
      subject = i,
      a1 = a1_vec,
      b1 = b1_vec,
      c = c_vec,
      rss = rss_vec,
      time_to_peak = hrf_grid$time_to_peak[grid_idx],
      FWHM = hrf_grid$FWHM[grid_idx],
      mask = subj_mask
    )
  }

  best_params_df <- do.call(rbind, best_params_list)

  if (verbose > 0) {
    cat("Built best_params_df:", nrow(best_params_df), "rows\n")
  }

  return(best_params_df)
}


#' Snap t2p/fwhm Values to Nearest Grid Point
#'
#' Converts population-averaged t2p and fwhm values back to the nearest
#' (a1, b1) grid point for a fixed c value. Uses normalized Euclidean
#' distance in (t2p, fwhm) space to find the closest match.
#'
#' @param t2p Numeric vector of time-to-peak values.
#' @param fwhm Numeric vector of FWHM values.
#' @param hrf_grid HRF grid with precomputed time_to_peak and FWHM columns.
#' @param fixed_c The fixed c value to snap to.
#'
#' @return A data.frame with columns: a1, b1, c (snapped to grid).
#'
#' @keywords internal
snap_to_grid_t2p_fwhm <- function(t2p, fwhm, hrf_grid, fixed_c) {
  # Filter grid to matching c value
  grid_subset <- hrf_grid[abs(hrf_grid$c - fixed_c) < 1e-6, ]

  # Normalization ranges
  t2p_range <- max(grid_subset$time_to_peak) - min(grid_subset$time_to_peak)
  fwhm_range <- max(grid_subset$FWHM) - min(grid_subset$FWHM)

  # Build matrices for vectorized distance computation
  # grid_t2p and grid_fwhm are length-G vectors (G = grid subset size)
  grid_t2p <- grid_subset$time_to_peak
  grid_fwhm <- grid_subset$FWHM

  # Compute distance matrix: (n_voxels x G)
  # Each row is a voxel, each column is a grid point
  t2p_dist <- outer(t2p, grid_t2p, function(a, b) ((a - b) / t2p_range)^2)
  fwhm_dist <- outer(fwhm, grid_fwhm, function(a, b) ((a - b) / fwhm_range)^2)
  dist_matrix <- sqrt(t2p_dist + fwhm_dist)

  # Find nearest grid point for each voxel
  nearest_idx <- apply(dist_matrix, 1, which.min)

  data.frame(
    a1 = grid_subset$a1[nearest_idx],
    b1 = grid_subset$b1[nearest_idx],
    c = fixed_c
  )
}


#' Snap Parameters to Nearest Grid Point
#'
#' Snaps arbitrary (a1, b1, c) values to the nearest valid HRF grid point.
#' Matches c value first, then finds the nearest (a1, b1) by normalized
#' Euclidean distance.
#'
#' @param a1 Numeric vector of a1 values to snap.
#' @param b1 Numeric vector of b1 values to snap.
#' @param c Numeric vector of c values to snap.
#' @param hrf_grid HRF grid data.frame with a1, b1, c columns.
#'
#' @return A data.frame with columns: a1, b1, c (snapped to grid).
#'
#' @keywords internal
snap_to_grid <- function(a1, b1, c, hrf_grid) {
  # Get unique grid combos
  grid_combos <- unique(hrf_grid[, c("a1", "b1", "c")])

  # Normalization ranges
  a1_range <- max(grid_combos$a1) - min(grid_combos$a1)
  b1_range <- max(grid_combos$b1) - min(grid_combos$b1)

  # Get unique c values in grid
  unique_c <- unique(grid_combos$c)

  n <- length(a1)
  a1_out <- numeric(n)
  b1_out <- numeric(n)
  c_out <- numeric(n)

  # Process by c value for efficiency
  for (c_val in unique_c) {
    # Which inputs match this c?
    c_match <- abs(c - c_val) < 1e-6
    if (!any(c_match)) next

    # Grid points for this c
    grid_c <- grid_combos[abs(grid_combos$c - c_val) < 1e-6, ]

    # Vectorized distance: inputs matching this c vs grid points for this c
    a1_sub <- a1[c_match]
    b1_sub <- b1[c_match]

    a1_dist <- outer(a1_sub, grid_c$a1, function(a, b) ((a - b) / a1_range)^2)
    b1_dist <- outer(b1_sub, grid_c$b1, function(a, b) ((a - b) / b1_range)^2)
    dist_matrix <- sqrt(a1_dist + b1_dist)

    nearest_idx <- apply(dist_matrix, 1, which.min)
    a1_out[c_match] <- grid_c$a1[nearest_idx]
    b1_out[c_match] <- grid_c$b1[nearest_idx]
    c_out[c_match] <- c_val
  }

  data.frame(a1 = a1_out, b1 = b1_out, c = c_out)
}


#' Extract a xifti geometry template from a workingHRF result
#'
#' Pulls the first successful subject's \code{pvalF_xii} so it can be reused as
#' a geometry template by \code{plot.regularizeHRFs} helpers without requiring
#' the caller to load a CIFTI separately. Called once at the end of
#' \code{regularize_allHRFs()} and stashed via \code{attr(result, "xii")}.
#'
#' @param workingHRF_results Output of \code{fit_workingHRF()}.
#' @return A \code{xifti} object, or \code{NULL} if no usable template is found.
#' @keywords internal
extract_xii_template <- function(workingHRF_results) {
  if (is.null(workingHRF_results$subject_results)) return(NULL)
  for (sr in workingHRF_results$subject_results) {
    if (!is.null(sr$status) && sr$status == "success") {
      xii <- sr$glm_results$pvalF_xii
      if (!is.null(xii)) return(sanitize_template_meta(xii))
    }
  }
  NULL
}


# TODO(remove-once-multiGLM-fixed): pvalF_xii (and Fstat_xii) in fit_workingHRF
# inherit `field_names` (length 8) and have empty intent because
# multiGLM.R:323 sets them unconditionally. Structurally they're 1-col
# dscalar-shaped, so we patch the column-level meta here to keep downstream
# consumers like smooth_cifti happy. When multiGLM is fixed, delete this
# helper and the call in extract_xii_template.
sanitize_template_meta <- function(xii) {
  xii$meta$cifti$names  <- "param"
  xii$meta$cifti$intent <- 3006L  # dscalar
  xii
}


#' Unmask pop_avg (fill + smooth t2p_mean / fwhm_mean) before grid snapping
#'
#' Builds full-cortex xifti maps from \code{pop_avg$t2p_mean} and
#' \code{pop_avg$fwhm_mean}, fills and smooths them via \code{unmask_xifti},
#' then returns an expanded \code{pop_avg} with one row per non-\code{NA}
#' cortex voxel.
#'
#' @param pop_avg The data.frame produced by Step 3 of
#'   \code{regularize_allHRFs} (cols: voxel, t2p_mean, fwhm_mean).
#' @param workingHRF_results Used to extract a xifti template.
#' @param mask_prop_NA Population mask vector from
#'   \code{workingHRF_results$activation_masks$mask_prop_NA}.
#' @param method,surf_FWHM,impute_FUN,impute_mask See \code{unmask_xifti}.
#' @param verbose Print step info.
#' @return The expanded \code{pop_avg} data.frame.
#' @keywords internal
unmask_pop_avg <- function(pop_avg, workingHRF_results, mask_prop_NA,
                            method = "median", surf_FWHM = 4,
                            impute_FUN = function(x) mean(x, na.rm = TRUE),
                            impute_mask = NULL,
                            verbose = 1) {
  if (verbose > 0) cat("\nUnmasking pop_avg (method=", method,
                       ", surf_FWHM=", surf_FWHM, "mm)...\n", sep = "")
  xii_template <- extract_xii_template(workingHRF_results)
  if (is.null(xii_template)) {
    stop("No xifti template available from workingHRF_results; cannot unmask.")
  }
  N <- length(mask_prop_NA)
  out <- data.frame(voxel = seq_len(N))
  for (col in c("t2p_mean", "fwhm_mean")) {
    v <- rep(NA, N); v[pop_avg$voxel] <- pop_avg[[col]]
    xii <- ciftiTools::newdata_xifti(xii_template, v)
    xii <- unmask_xifti(xii, method = method, surf_FWHM = surf_FWHM,
                        impute_FUN = impute_FUN, impute_mask = impute_mask)
    out[[col]] <- as.matrix(xii)[, 1]
  }
  out <- out[stats::complete.cases(out), ]
  if (verbose > 0) cat("Unmasked pop_avg covers", nrow(out), "voxels (was ",
                       nrow(pop_avg), ")\n", sep = "")
  out
}


#' Fill missing values in a xifti and smooth
#'
#' Method \code{"median"} replaces all \code{NA}s with the global median of
#' the in-mask values (one scalar sprayed everywhere). Method \code{"impute"}
#' uses \code{ciftiTools::impute_xifti} to fill from face-sharing surface
#' neighbors. After filling, the map is smoothed with
#' \code{ciftiTools::smooth_cifti}.
#'
#' @param xii A \code{xifti} object with \code{NA} values to fill.
#' @param method One of \code{"median"} (default), \code{"impute"}.
#' @param surf_FWHM Cortex smoothing FWHM in mm (required, no default — caller
#'   sets policy; \code{unmask_pop_avg} is the canonical entry point).
#'   \code{0} skips smoothing.
#' @param impute_FUN For \code{method="impute"}: function applied to neighbor
#'   values. Default is mean (ignoring \code{NA}).
#' @param impute_mask For \code{method="impute"}: optional logical mask of
#'   voxels to impute. Default \code{NULL} = impute every \code{NA} location.
#' @return The filled, smoothed \code{xifti}.
#' @keywords internal
unmask_xifti <- function(xii,
                          method = c("median", "impute"),
                          surf_FWHM,
                          impute_FUN = function(x) mean(x, na.rm = TRUE),
                          impute_mask = NULL) {
  method <- match.arg(method)
  if (method == "median") {
    m   <- as.matrix(xii)
    med <- stats::median(m, na.rm = TRUE)
    m[is.na(m)] <- med
    xii <- ciftiTools::newdata_xifti(xii, m)
  } else if (method == "impute") {
    if (is.null(xii$surf$cortex_left) || is.null(xii$surf$cortex_right)) {
      xii <- ciftiTools::add_surf(xii, surfL = "inflated", surfR = "inflated")
    }
    xii <- ciftiTools::impute_xifti(xii, mask = impute_mask, impute_FUN = impute_FUN)
  }
  if (surf_FWHM > 0) {
    xii <- ciftiTools::smooth_cifti(xii, surf_FWHM = surf_FWHM)
  }
  xii
}