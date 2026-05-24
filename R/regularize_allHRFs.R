#' Regularize All HRFs Using POCv2 Algorithm
#'
#' Creates a population-averaged HRF map, generates candidate offset maps,
#' and selects the best candidate per subject based on weighted RSS.
#'
#' @param workingHRF_results Results from \code{fit_workingHRF()}.
#' @param allHRF_results Results from \code{fit_allHRFs()}.
#' @param BOLD Character vector of CIFTI file paths. Required for refit mode.
#' @param EVs List of event data per subject. Required for refit mode.
#' @param nuisance Character vector of nuisance file paths (or NULL). For refit mode.
#' @param TR Numeric. Repetition time in seconds. Required for refit mode.
#' @param a1_offsets Numeric vector of a1 offsets for candidate maps.
#' @param b1_offsets Numeric vector of b1 offsets for candidate maps.
#' @param seffects Logical. If \code{TRUE} (default), fit all candidate maps
#'   per subject and pick a winning candidate per subject (the adapted-HRF
#'   workflow). If \code{FALSE}, return only the population average and skip
#'   per-subject candidate selection.
#' @param onsets Logical. Include onset regressors (refit mode only).
#' @param offsets Logical. Include offset regressors (refit mode only).
#' @param n_cores Integer. Number of parallel worker processes for refit mode
#'   (1 = sequential). Ignored in lookup mode (\code{save_rss = TRUE}).
#' @param log_dir Directory for cluster log files (refit + parallel only).
#' @param verbose Integer verbosity level.
#'
#' @return A list with class \code{"regularizeHRFs"} containing:
#'   \describe{
#'     \item{pop_avg}{Population average HRF map per voxel}
#'     \item{best_params_df}{Best HRF per subject per voxel (steps 1-2)}
#'     \item{candidate_maps}{List of 25 candidate maps}
#'     \item{subject_results}{Per-subject candidate selection results}
#'     \item{winning_c}{Winning c value}
#'     \item{c_votes}{Vote counts per c value}
#'     \item{hrf_grid}{HRF grid with t2p/fwhm columns}
#'     \item{mask_prop_NA}{Population activation mask}
#'   }
#'
#' @export
regularize_allHRFs <- function(workingHRF_results,
                                allHRF_results,
                                BOLD = NULL,
                                EVs = NULL,
                                nuisance = NULL,
                                TR = NULL,
                                a1_offsets = c(-2, -1, 0, 1, 2),
                                b1_offsets = c(-0.5, -0.25, 0, 0.25, 0.5),
                                seffects = TRUE,
                                onsets = FALSE,
                                offsets = FALSE,
                                n_cores = 1,
                                log_dir = "logs",
                                verbose = 1) {

  # Detect RSS mode
  save_rss <- allHRF_results$call_info$save_rss
  if (is.null(save_rss)) {
    if (verbose > 0) cat("No save_rss flag found in allHRF_results (old data). Defaulting to refit mode.\n")
    save_rss <- FALSE
  }

  if (save_rss) {
    if (verbose > 0) cat("RSS lookup mode: loading pre-computed RSS from .qs files (fast)\n")
  } else {
    if (verbose > 0) cat("Refit mode: will run multiGLM for each subject (slow, requires BOLD/EVs)\n")
    if (is.null(BOLD) || is.null(EVs)) stop("BOLD and EVs are required when save_rss = FALSE (refit mode)")
    if (is.null(TR)) stop("TR is required when save_rss = FALSE (refit mode)")
  }

  # Pack into session_data for internal use
  session_data <- if (!is.null(BOLD)) {
    list(BOLD_files = BOLD, EVs_list = EVs, nuisance_files = nuisance)
  } else {
    NULL
  }

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

  # Steps 5-6: Create candidate maps
  if (verbose > 0) cat("\nCreating candidate maps...\n")
  cm_result <- create_candidate_maps(pop_avg, hrf_grid, a1_offsets, b1_offsets, verbose = verbose)

  # Early return if seffects=FALSE (skip candidate fitting)
  if (!seffects) {
    if (verbose > 0) cat("\nseffects=FALSE: returning population average only.\n")
    result <- list(
      pop_avg = pop_avg,
      best_params_df = best_params_df,
      candidate_maps = cm_result$candidate_maps,
      offset_combos = cm_result$offset_combos,
      subject_results = NULL,
      winning_c = winning_result$winning_c,
      c_votes = winning_result$c_votes,
      hrf_grid = hrf_grid,
      mask_prop_NA = mask_prop_NA
    )
    class(result) <- "regularizeHRFs"
    if (verbose > 0) cat("regularize_allHRFs complete.\n")
    return(result)
  }

  # Steps 7-8: Fit candidates and pick winner per subject
  if (verbose > 0) cat("\nFitting candidates to subjects...\n")

  if (save_rss) {
    subject_results <- fit_candidate_maps_lookup(
      cm_result$candidate_maps, allHRF_results, workingHRF_results, verbose = verbose
    )
  } else {
    subject_results <- fit_candidate_maps_refit(
      cm_result$candidate_maps, hrf_grid, session_data,
      TR = TR, onsets = onsets, offsets = offsets,
      n_cores = n_cores, log_dir = log_dir, verbose = verbose
    )
  }

  # Build result
  result <- list(
    pop_avg = pop_avg,
    best_params_df = best_params_df,
    candidate_maps = cm_result$candidate_maps,
    offset_combos = cm_result$offset_combos,
    subject_results = subject_results,
    winning_c = winning_result$winning_c,
    c_votes = winning_result$c_votes,
    hrf_grid = hrf_grid,
    mask_prop_NA = mask_prop_NA
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


#' Create Candidate Maps from Population Average
#'
#' Generates candidate HRF maps by applying offset combinations to the
#' population average map. Each candidate is a shifted version of the
#' population map, snapped back to valid grid points.
#'
#' @param pop_avg Data.frame from population averaging with columns:
#'   voxel, a1_snapped, b1_snapped, c_snapped.
#' @param hrf_grid HRF grid with precomputed time_to_peak and FWHM columns.
#' @param a1_offsets Numeric vector of a1 offsets (default: c(-2,-1,0,1,2)).
#' @param b1_offsets Numeric vector of b1 offsets (default: c(-0.5,-0.25,0,0.25,0.5)).
#' @param verbose Integer verbosity level.
#'
#' @return A list with:
#'   \describe{
#'     \item{candidate_maps}{List of data.frames, one per offset combo. Each has
#'       columns: voxel, a1, b1, c, model_idx, time_to_peak, FWHM, offset_id,
#'       a1_offset, b1_offset.}
#'     \item{offset_combos}{Data.frame of offset combinations used.}
#'   }
#'
#' @keywords internal
create_candidate_maps <- function(pop_avg, hrf_grid,
                                   a1_offsets = c(-2, -1, 0, 1, 2),
                                   b1_offsets = c(-0.5, -0.25, 0, 0.25, 0.5),
                                   verbose = 1) {
  offset_combos <- expand.grid(a1_offset = a1_offsets, b1_offset = b1_offsets)
  n_candidates <- nrow(offset_combos)

  if (verbose > 0) cat("Creating", n_candidates, "candidate maps\n")

  # Valid grid ranges
  a1_min <- min(hrf_grid$a1)
  a1_max <- max(hrf_grid$a1)
  b1_min <- min(hrf_grid$b1)
  b1_max <- max(hrf_grid$b1)

  # Build model index lookup (unique HRF combos across all candidates)
  # We'll assign model indices after collecting all unique combos
  hrf_grid$key <- paste(hrf_grid$a1, hrf_grid$b1, hrf_grid$c, sep = "_")

  candidate_maps <- vector("list", n_candidates)

  for (i in seq_len(n_candidates)) {
    a1_off <- offset_combos$a1_offset[i]
    b1_off <- offset_combos$b1_offset[i]

    # Apply offset and clamp to valid range
    a1_shifted <- pmin(pmax(pop_avg$a1_snapped + a1_off, a1_min), a1_max)
    b1_shifted <- pmin(pmax(pop_avg$b1_snapped + b1_off, b1_min), b1_max)
    c_vals <- pop_avg$c_snapped

    # Snap to nearest valid grid point
    snapped <- snap_to_grid(a1_shifted, b1_shifted, c_vals, hrf_grid)

    # Look up model index and metrics from hrf_grid
    keys <- paste(snapped$a1, snapped$b1, snapped$c, sep = "_")
    grid_idx <- match(keys, hrf_grid$key)

    candidate_maps[[i]] <- data.frame(
      voxel = pop_avg$voxel,
      a1 = snapped$a1,
      b1 = snapped$b1,
      c = snapped$c,
      model_idx = grid_idx,
      time_to_peak = hrf_grid$time_to_peak[grid_idx],
      FWHM = hrf_grid$FWHM[grid_idx],
      offset_id = i,
      a1_offset = a1_off,
      b1_offset = b1_off
    )
  }

  if (verbose > 0) {
    # Count unique HRF combos across all candidates
    all_keys <- unique(unlist(lapply(candidate_maps, function(cm) {
      paste(cm$a1, cm$b1, cm$c, sep = "_")
    })))
    cat("Total unique HRF combos across candidates:", length(all_keys), "\n")
  }

  return(list(
    candidate_maps = candidate_maps,
    offset_combos = offset_combos
  ))
}


#' Fit Candidate Maps via RSS Lookup
#'
#' Uses pre-computed RSS matrices stored in .qs files to evaluate
#' candidate maps for each subject. Much faster than refitting.
#' Requires \code{fit_allHRFs} to have been run with \code{save_rss = TRUE}.
#'
#' @param candidate_maps List of candidate map data.frames from
#'   \code{create_candidate_maps()}.
#' @param allHRF_results Results from \code{fit_allHRFs(save_rss = TRUE)}.
#' @param workingHRF_results Results from \code{fit_workingHRF()}.
#' @param verbose Integer verbosity level.
#'
#' @return A data.frame with one row per subject, columns:
#'   subject, winning_candidate_id, winning_a1_offset, winning_b1_offset,
#'   winning_weighted_RSS, plus RSS for all candidates.
#'
#' @keywords internal
fit_candidate_maps_lookup <- function(candidate_maps, allHRF_results,
                                      workingHRF_results, verbose = 1) {
  result_paths <- attr(allHRF_results, "result_paths")
  n_subjects <- length(result_paths)
  n_candidates <- length(candidate_maps)

  if (verbose > 0) cat("RSS lookup mode:", n_subjects, "subjects,", n_candidates, "candidates\n")

  subject_results_list <- vector("list", n_subjects)

  for (i in seq_len(n_subjects)) {
    if (verbose > 1 && i %% 100 == 0) cat("  Subject", i, "/", n_subjects, "\n")

    # Load RSS from .qs file
    qs_obj <- qs2::qs_read(result_paths[i])
    RSS <- rbind(
      qs_obj$glm_result$mGLM0s$cortexL$RSS,
      qs_obj$glm_result$mGLM0s$cortexR$RSS
    )
    rm(qs_obj)

    # Fstat from fit_workingHRF (canonical HRF activation strength)
    Fstat <- c(
      workingHRF_results$subject_results[[i]]$glm_results$mGLM0s$cortexL$Fstat,
      workingHRF_results$subject_results[[i]]$glm_results$mGLM0s$cortexR$Fstat
    )

    # Subset Fstat to candidate map voxels only
    pop_voxels <- candidate_maps[[1]]$voxel
    fstat_sq <- Fstat[pop_voxels]^2
    fstat_sq_sum <- sum(fstat_sq, na.rm = TRUE)

    # Compute RSS per candidate
    rss_per_candidate <- numeric(n_candidates)
    weighted_rss_per_candidate <- numeric(n_candidates)

    for (j in seq_len(n_candidates)) {
      voxel_idx <- candidate_maps[[j]]$voxel
      model_idx <- candidate_maps[[j]]$model_idx
      voxel_rss <- RSS[cbind(voxel_idx, model_idx)]
      rss_per_candidate[j] <- sum(voxel_rss, na.rm = TRUE)
      weighted_rss_per_candidate[j] <- sum(voxel_rss * fstat_sq, na.rm = TRUE) / fstat_sq_sum
    }

    # Find winner
    winner <- which.min(weighted_rss_per_candidate)

    # Build row for this subject
    row <- data.frame(
      subject = i,
      winning_candidate_id = winner,
      winning_a1_offset = candidate_maps[[winner]]$a1_offset[1],
      winning_b1_offset = candidate_maps[[winner]]$b1_offset[1],
      winning_weighted_RSS = weighted_rss_per_candidate[winner]
    )

    # Add all candidate RSS values
    for (j in seq_len(n_candidates)) {
      row[[paste0("weighted_RSS_", j)]] <- weighted_rss_per_candidate[j]
    }

    subject_results_list[[i]] <- row
  }

  result <- do.call(rbind, subject_results_list)

  if (verbose > 0) {
    winner_counts <- table(result$winning_candidate_id)
    cat("Winner distribution:\n")
    print(winner_counts)
  }

  return(result)
}


#' Fit Candidate Maps via Refitting (multiGLM)
#'
#' Fits candidate HRF models to each subject's BOLD data via fresh multiGLM.
#' Used when \code{fit_allHRFs} was run with \code{save_rss = FALSE} and RSS
#' is not available in .qs files. Slower than lookup but doesn't require
#' stored RSS.
#'
#' @param candidate_maps List of candidate map data.frames from
#'   \code{create_candidate_maps()}.
#' @param hrf_grid HRF grid with a1, b1, c, a2, b2 columns.
#' @param session_data List with BOLD_files, EVs_list, nuisance_files.
#' @param TR Numeric. Repetition time in seconds.
#' @param brainstructures Character vector. Brain structures to model.
#' @param hpf Numeric. High-pass filter cutoff.
#' @param onsets Logical. Include onset regressors.
#' @param offsets Logical. Include offset regressors.
#' @param n_cores Integer. Number of parallel worker processes (1 = sequential).
#' @param log_dir Directory for cluster log files (parallel mode only).
#' @param verbose Integer verbosity level.
#'
#' @return A data.frame with same structure as \code{fit_candidate_maps_lookup}.
#'
#' @keywords internal
fit_candidate_maps_refit <- function(candidate_maps, hrf_grid,
                                     session_data, TR,
                                     brainstructures = c("left", "right"),
                                     hpf = 0.01, onsets = FALSE, offsets = FALSE,
                                     n_cores = 1, log_dir = "logs",
                                     verbose = 1) {
  n_subjects <- length(session_data$BOLD_files)
  n_candidates <- length(candidate_maps)

  # Collect unique HRF model indices across all candidate maps
  all_model_idx <- sort(unique(unlist(lapply(candidate_maps, function(cm) unique(cm$model_idx)))))
  unique_grid <- hrf_grid[all_model_idx, ]
  n_models <- nrow(unique_grid)
  idx_map <- setNames(seq_len(n_models), all_model_idx)
  pop_voxels <- candidate_maps[[1]]$voxel

  if (verbose > 0) cat("Refit mode:", n_subjects, "subjects,", n_candidates, "candidates,", n_models, "unique models using", n_cores, "cores\n")

  if (n_cores > 1) {
    subject_rows <- run_parallel_subjects_refit(
      candidate_maps, hrf_grid, session_data, TR,
      brainstructures, hpf, onsets, offsets,
      unique_grid, idx_map, pop_voxels, n_candidates,
      verbose, n_cores, log_dir
    )
  } else {
    if (verbose > 0) cat("Using sequential processing\n")
    subject_rows <- run_sequential_subjects_refit(
      candidate_maps, hrf_grid, session_data, TR,
      brainstructures, hpf, onsets, offsets,
      unique_grid, idx_map, pop_voxels, n_candidates,
      verbose
    )
  }

  result <- do.call(rbind, subject_rows)

  if (verbose > 0) {
    cat("Winner distribution:\n")
    print(table(result$winning_candidate_id))
  }

  return(result)
}


#' Process candidate subjects in parallel using cluster workers
#'
#' Sets up a parallel cluster and processes all subjects simultaneously using
#' load-balanced \code{parLapplyLB()}. Each worker calls
#' \code{process_candidate_subject_refit()} for one subject.
#'
#' @inheritParams fit_candidate_maps_refit
#' @param unique_grid Compact HRF grid restricted to model indices referenced by
#'   candidate maps.
#' @param idx_map Named integer vector mapping global model_idx (character key)
#'   to compact unique_grid row index.
#' @param pop_voxels Integer vector of population-mask voxels (rows of RSS that
#'   contribute to the weighted RSS sum).
#' @param n_candidates Integer. Number of candidate maps.
#'
#' @return List of single-row data.frames, one per subject.
#'
#' @keywords internal
run_parallel_subjects_refit <- function(candidate_maps, hrf_grid, session_data, TR,
                                        brainstructures, hpf, onsets, offsets,
                                        unique_grid, idx_map, pop_voxels, n_candidates,
                                        verbose, n_cores, log_dir) {

  if (verbose > 0) cat("Setting up parallel cluster with", n_cores, "cores.\n")

  n_workers <- n_cores

  if (!dir.exists(log_dir)) { dir.create(log_dir, recursive = TRUE) }
  log_file <- file.path(log_dir, sprintf("regularize_refit_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
  message("Cluster logs will be saved to: ", normalizePath(log_file, mustWork = FALSE))

  cl <- parallel::makeCluster(n_workers, outfile = log_file)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  vars_to_export <- c(
    "candidate_maps", "hrf_grid", "session_data", "TR",
    "brainstructures", "hpf", "onsets", "offsets",
    "unique_grid", "idx_map", "pop_voxels", "n_candidates", "verbose"
  )
  setup_parallel_cluster(cl, verbose, vars_to_export)

  if (verbose > 0) cat("Processing subjects in parallel\n")

  subject_rows <- parallel::parLapplyLB(cl, seq_len(length(session_data$BOLD_files)), function(i) {
    process_candidate_subject_refit(
      subject_idx = i,
      candidate_maps = candidate_maps, hrf_grid = hrf_grid, session_data = session_data,
      TR = TR, brainstructures = brainstructures, hpf = hpf,
      onsets = onsets, offsets = offsets,
      unique_grid = unique_grid, idx_map = idx_map,
      pop_voxels = pop_voxels, n_candidates = n_candidates
    )
  })

  if (verbose > 0) cat("Parallel processing completed\n")
  return(subject_rows)
}


#' Process candidate subjects sequentially (single-threaded)
#'
#' Iterates over subjects with \code{lapply()} and calls
#' \code{process_candidate_subject_refit()} for each. Identical math to the
#' parallel path.
#'
#' @inheritParams run_parallel_subjects_refit
#'
#' @return List of single-row data.frames, one per subject.
#'
#' @keywords internal
run_sequential_subjects_refit <- function(candidate_maps, hrf_grid, session_data, TR,
                                          brainstructures, hpf, onsets, offsets,
                                          unique_grid, idx_map, pop_voxels, n_candidates,
                                          verbose) {
  n_subjects <- length(session_data$BOLD_files)
  lapply(seq_len(n_subjects), function(i) {
    if (verbose > 0) cat("  Subject", i, "/", n_subjects, "\n")
    process_candidate_subject_refit(
      subject_idx = i,
      candidate_maps = candidate_maps, hrf_grid = hrf_grid, session_data = session_data,
      TR = TR, brainstructures = brainstructures, hpf = hpf,
      onsets = onsets, offsets = offsets,
      unique_grid = unique_grid, idx_map = idx_map,
      pop_voxels = pop_voxels, n_candidates = n_candidates
    )
  })
}


#' Process one subject through refit-mode candidate selection
#'
#' Loads BOLD + nuisance for one subject, builds the design tensor across the
#' compact set of unique HRF models, fits multiGLM, computes weighted RSS per
#' candidate map, and returns the winner row. Called by both parallel and
#' sequential refit drivers.
#'
#' @param subject_idx Integer subject index.
#' @inheritParams run_parallel_subjects_refit
#'
#' @return A single-row data.frame with columns subject, winning_candidate_id,
#'   winning_a1_offset, winning_b1_offset, winning_weighted_RSS, and
#'   weighted_RSS_1..n_candidates.
#'
#' @keywords internal
process_candidate_subject_refit <- function(subject_idx,
                                            candidate_maps, hrf_grid, session_data, TR,
                                            brainstructures, hpf, onsets, offsets,
                                            unique_grid, idx_map, pop_voxels, n_candidates) {
  i <- subject_idx
  n_models <- nrow(unique_grid)

  # Load BOLD + nuisance
  BOLD_xii <- ciftiTools::read_cifti(
    session_data$BOLD_files[[i]],
    brainstructures = brainstructures,
    resamp_res = 10000
  )
  nuisance <- as.matrix(utils::read.table(session_data$nuisance_files[[i]], header = FALSE))
  nT <- ncol(BOLD_xii)

  # Build canonical + per-model design matrices
  design_canonical <- make_design(
    EVs = session_data$EVs_list[[i]], nTime = nT, TR = TR, dHRF = 0,
    onset = onsets, offset = offsets
  )$design
  nK <- ncol(design_canonical)

  design_3D <- array(NA, dim = c(nT, nK, n_models))
  for (j in seq_len(n_models)) {
    taper_start_j <- get_taper_start(
      a1 = unique_grid$a1[j], b1 = unique_grid$b1[j],
      a2 = unique_grid$a2[j], b2 = unique_grid$b2[j],
      c = unique_grid$c[j], TR = TR, deriv = 0
    )
    design_3D[,,j] <- make_design(
      EVs = session_data$EVs_list[[i]], nTime = nT, TR = TR, dHRF = 0,
      onset = onsets, offset = offsets,
      taper_start = taper_start_j,
      a1 = unique_grid$a1[j], b1 = unique_grid$b1[j],
      c = unique_grid$c[j],
      a2 = unique_grid$a2[j], b2 = unique_grid$b2[j]
    )$design
  }

  # Fit multiGLM
  glm_result <- multiGLM(
    BOLD = BOLD_xii, brainstructures = brainstructures,
    resamp_res = NULL, design = design_3D,
    design_canonical = design_canonical,
    nuisance = nuisance, scrub = NULL, TR = TR, hpf = hpf
  )

  RSS <- rbind(glm_result$mGLM0s$cortexL$RSS, glm_result$mGLM0s$cortexR$RSS)
  Fstat <- c(glm_result$mGLM0s$cortexL$Fstat, glm_result$mGLM0s$cortexR$Fstat)

  fstat_sq <- Fstat[pop_voxels]^2
  fstat_sq_sum <- sum(fstat_sq, na.rm = TRUE)

  # Weighted RSS per candidate
  weighted_rss_per_candidate <- numeric(n_candidates)
  for (j in seq_len(n_candidates)) {
    compact_idx <- idx_map[as.character(candidate_maps[[j]]$model_idx)]
    voxel_rss <- RSS[cbind(pop_voxels, compact_idx)]
    weighted_rss_per_candidate[j] <- sum(voxel_rss * fstat_sq, na.rm = TRUE) / fstat_sq_sum
  }

  winner <- which.min(weighted_rss_per_candidate)
  row <- data.frame(
    subject = i,
    winning_candidate_id = winner,
    winning_a1_offset = candidate_maps[[winner]]$a1_offset[1],
    winning_b1_offset = candidate_maps[[winner]]$b1_offset[1],
    winning_weighted_RSS = weighted_rss_per_candidate[winner]
  )
  for (j in seq_len(n_candidates)) {
    row[[paste0("weighted_RSS_", j)]] <- weighted_rss_per_candidate[j]
  }
  row
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
                            method = "median", surf_FWHM = 5,
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