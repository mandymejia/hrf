#' Build the Population HRF Template from fit_allHRFs Results
#'
#' Aggregates per-subject 25-candidate fits into a population template using
#' summed-RSS: for each voxel, sum per-candidate RSS across all activated
#' subjects, then pick the grid point with the lowest sum. Winner stays on
#' the grid (no snap, no smoothing). Non-activated voxels get filled with the
#' cortex-wide modal winning grid point.
#'
#' Requires \code{fit_allHRFs()} to have been run with \code{save_rss = TRUE}
#' and the per-subject RSS matrices to be retained on
#' \code{allHRF_results$subject_results$glm_result$mGLM0s$cortex[LR]$RSS}.
#'
#' @param workingHRF_results Output of \code{fit_workingHRF()}. Supplies
#'   the group activation mask (\code{$activation_masks}) and a xifti
#'   geometry template.
#' @param allHRF_results Output of \code{fit_allHRFs()}. Supplies the HRF
#'   grid and per-subject RSS matrices.
#' @param verbose Integer verbosity.
#'
#' @return A list with class \code{"regularizeHRFs"}:
#'   \describe{
#'     \item{pop_best}{Per-voxel winner: columns \code{voxel, a1,
#'       b1, c} (values come straight off the grid).}
#'     \item{winning_c}{Modal c across pop_best (scalar).}
#'     \item{hrf_grid}{HRF grid with \code{time_to_peak} + \code{FWHM} columns added.}
#'     \item{mask_prop_NA}{Population activation mask (passed through).}
#'     \item{winning_k}{Grid index per voxel picked by argmin; \code{mode_k}
#'       for non-activated voxels after modal fill.}
#'   }
#'
#' @export
regularize_allHRFs <- function(workingHRF_results, allHRF_results, verbose = 1) {

  # Step 1: attach t2p / FWHM columns to the HRF grid (downstream consumers use them).
  hrf_grid <- add_grid_metrics(allHRF_results$hrf_grid)

  # Step 2: stream per-candidate RSS sums across activated subjects.
  agg <- sum_rss_across_subjects(
    subject_results = allHRF_results$subject_results,
    subj_masks      = workingHRF_results$activation_masks$masks,
    n_candidates    = nrow(hrf_grid),
    verbose         = verbose
  )

  mask_prop_NA <- workingHRF_results$activation_masks$mask_prop_NA

  # Step 3: pick the winning grid point per pop-mask voxel (argmin summed RSS).
  winning_k <- pick_winning_grid_point(
    sum_rss        = agg$sum_rss,
    mask_prop_NA   = mask_prop_NA,
    n_contributing = agg$n_contributing,
    verbose        = verbose
  )

  # Step 4: fill non-activated voxels with the cortex-wide modal winner.
  winning_k <- modal_unmask(
    winning_k    = winning_k,
    mask_prop_NA = mask_prop_NA,
    hrf_grid     = hrf_grid,
    verbose      = verbose
  )

  # Step 5: assemble the regularize result (pop_best + xifti template).
  build_regularize_result(
    winning_k          = winning_k,
    hrf_grid           = hrf_grid,
    mask_prop_NA       = mask_prop_NA,
    workingHRF_results = workingHRF_results,
    verbose            = verbose
  )
}


# ---- Step helpers ------------------------------------------------------------


#' Step 1: attach t2p / FWHM columns to the HRF grid
#' @keywords internal
add_grid_metrics <- function(hrf_grid) {
  m <- get_hrf_metrics(hrf_grid$a1, hrf_grid$b1, hrf_grid$c)
  hrf_grid$time_to_peak <- m$time_to_peak
  hrf_grid$FWHM         <- m$FWHM
  hrf_grid
}


#' Step 2: sum per-candidate RSS across activated subjects
#'
#' Streams over subjects to build \code{sum_rss[voxel, candidate]}, restricted
#' to each subject's per-voxel activation mask. Also tracks
#' \code{n_contributing[voxel]} for diagnostics.
#'
#' @keywords internal
sum_rss_across_subjects <- function(subject_results, subj_masks, n_candidates, verbose = 1) {
  n_voxels   <- nrow(subj_masks)
  n_subjects <- ncol(subj_masks)
  stopifnot(length(subject_results) == n_subjects)

  if (verbose > 0) {
    cat(sprintf("Streaming RSS: %d subjects x %d voxels x %d candidates\n",
                n_subjects, n_voxels, n_candidates))
  }

  sum_rss        <- matrix(0, nrow = n_voxels, ncol = n_candidates)
  n_contributing <- integer(n_voxels)

  for (i in seq_len(n_subjects)) {
    if (verbose > 0 && (i %% 100 == 0 || i == n_subjects)) {
      cat(sprintf("  subject %d/%d\n", i, n_subjects))
    }
    subj <- subject_results[[i]]
    if (is.null(subj) || is.null(subj$glm_result)) next
    RSS <- rbind(subj$glm_result$mGLM0s$cortexL$RSS,
                 subj$glm_result$mGLM0s$cortexR$RSS)
    if (is.null(RSS) || nrow(RSS) == 0L) {
      warning("Subject ", i,
              " has no RSS. Was fit_allHRFs run with save_rss = TRUE? Skipping.")
      next
    }
    mask_v <- as.logical(subj_masks[, i])
    valid  <- which(mask_v & !is.na(RSS[, 1]))
    if (length(valid) == 0) next

    sum_rss[valid, ]      <- sum_rss[valid, ] + RSS[valid, ]
    n_contributing[valid] <- n_contributing[valid] + 1L
  }

  list(sum_rss = sum_rss, n_contributing = n_contributing)
}


#' Step 3: pick the winning grid point per voxel via argmin summed RSS
#'
#' Only voxels in the pop-mask that had at least one contributing subject get
#' a winner; the rest stay \code{NA} (they'll be filled by modal_unmask).
#'
#' @keywords internal
pick_winning_grid_point <- function(sum_rss, mask_prop_NA, n_contributing, verbose = 1) {
  n_voxels        <- length(mask_prop_NA)
  pop_mask_voxels <- which(!is.na(mask_prop_NA))
  eligible        <- pop_mask_voxels[n_contributing[pop_mask_voxels] > 0L]

  winning_k <- rep(NA_integer_, n_voxels)
  if (length(eligible) > 0) {
    # max.col of negated sum = argmin of sum, vectorized over rows.
    winning_k[eligible] <- max.col(-sum_rss[eligible, , drop = FALSE], ties.method = "first")
  }
  if (verbose > 0) {
    cat(sprintf("Argmin filled %d/%d pop-mask voxels\n",
                sum(!is.na(winning_k)), length(pop_mask_voxels)))
  }
  winning_k
}


#' Step 4: modal-fill non-activated voxels with the cortex-wide winner
#'
#' Takes the single most-common winning grid point across all activated
#' voxels and assigns it to every voxel outside the pop-mask.
#'
#' @keywords internal
modal_unmask <- function(winning_k, mask_prop_NA, hrf_grid, verbose = 1) {
  n_candidates <- nrow(hrf_grid)
  activated    <- winning_k[!is.na(winning_k)]
  if (length(activated) == 0) return(winning_k)

  tab      <- tabulate(activated, nbins = n_candidates)
  mode_k   <- which.max(tab)
  unmask_v <- which(is.na(mask_prop_NA))
  winning_k[unmask_v] <- mode_k

  if (verbose > 0) {
    cat(sprintf("Modal fill: %d voxels -> grid[%d] = (a1=%.2f, b1=%.2f, c=%.4f)\n",
                length(unmask_v), mode_k,
                hrf_grid$a1[mode_k], hrf_grid$b1[mode_k], hrf_grid$c[mode_k]))
  }
  winning_k
}


#' Step 5: assemble the pop_best + regularize result object
#'
#' @keywords internal
build_regularize_result <- function(winning_k, hrf_grid, mask_prop_NA, workingHRF_results, verbose = 1) {
  # t2p_mean / fwhm_mean columns are kept for downstream compatibility with
  # plot.regularizeHRFs (references them by name). Values are looked up from
  # the grid at each voxel's winning candidate -- semantically "winner's t2p"
  # rather than "mean t2p across subjects".
  pop_best <- data.frame(
    voxel      = seq_along(winning_k),
    a1 = hrf_grid$a1[winning_k],
    b1 = hrf_grid$b1[winning_k],
    c  = hrf_grid$c[winning_k],
    t2p_mean   = hrf_grid$time_to_peak[winning_k],
    fwhm_mean  = hrf_grid$FWHM[winning_k]
  )

  valid_c   <- pop_best$c[!is.na(pop_best$c)]
  winning_c <- if (length(valid_c) == 0) NA_real_
               else as.numeric(names(sort(table(valid_c), decreasing = TRUE))[1])

  result <- list(
    pop_best      = pop_best,
    winning_c    = winning_c,
    hrf_grid     = hrf_grid,
    mask_prop_NA = mask_prop_NA,
    winning_k    = winning_k
  )
  attr(result, "xii") <- extract_xii_template(workingHRF_results)
  class(result) <- "regularizeHRFs"

  if (verbose > 0) cat("regularize_allHRFs complete.\n")
  result
}


# ---- Supporting helpers ------------------------------------------------------


#' Compute HRF Metrics (Time-to-Peak and FWHM)
#'
#' Converts HRF parameters (a1, b1, c) to interpretable metrics:
#' time-to-peak and full-width-at-half-maximum.
#'
#' @param a1 First peak parameter
#' @param b1 First scale parameter
#' @param c Undershoot amplitude
#' @param tapered Whether to apply tapering for c > 0 HRFs
#'
#' @return Data.frame with columns a1, b1, c, time_to_peak, FWHM
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


#' Snap continuous (a1, b1, c) to nearest grid point
#'
#' Used by \code{fit_bestHRF} to round personalized-mode HRFs back to a grid
#' point. Not used by \code{regularize_allHRFs} (which produces on-grid winners
#' directly).
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


#' Fill + smooth a xifti map
#'
#' Used by \code{plot.regularizeHRFs} to render continuous param maps: fills
#' masked-out voxels via a chosen method, then optionally smooths. Not used
#' by \code{regularize_allHRFs} aggregation any more (pop_rss produces winners
#' directly on the grid).
#'
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
