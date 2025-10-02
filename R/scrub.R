#' Convert logical vector to spike regressor matrix
#' @param scrub_logical Logical vector where TRUE = scrub this timepoint
#'
#' @inheritParams verbose_Param
#'
#' @return Matrix with one column per scrubbed timepoint
#' @keywords internal
logical_to_spike_regressors <- function(scrub_logical, verbose = 1) {
  n_scrub <- sum(scrub_logical)
  if (n_scrub == 0) {
    if (verbose > 1) cat("No timepoints marked for scrubbing\n")
    return(NULL)
  }

  if (verbose > 1) {
    cat("Converting logical scrub vector to spike regressors: ",
        n_scrub, " timepoints to scrub\n")
  }

  nT <- length(scrub_logical)
  spike_regressors <- matrix(0, nrow = nT, ncol = n_scrub)
  scrub_indices <- which(scrub_logical)

  for (i in seq_along(scrub_indices)) {
    spike_regressors[scrub_indices[i], i] <- 1
  }

  return(spike_regressors)
}

#' Validate spike regressor matrix format
#'
#' @inheritParams spike_matrix_Param
#' @inheritParams nT_Param
#' @keywords internal
validate_spike_regressors <- function(spike_matrix, nT) {
  if (nrow(spike_matrix) != nT) {
    stop("Spike regressor matrix has ", nrow(spike_matrix),
         " rows but should have ", nT, " rows to match timepoints.")
  }

  # Each column should sum to 1 and have exactly one 1 per column
  if (!all(colSums(spike_matrix) == 1)) {
    stop("Each column in spike regressor matrix should sum to 1")
  }

  if (!all(colSums(spike_matrix == 0) == nrow(spike_matrix) - 1)) {
    stop("Each column in spike regressor matrix should have exactly one 1 and the rest 0s")
  }

  return(TRUE)
}

#' Parse scrub input into spike regressor matrix
#'
#' @inheritParams scrub_Param
#' @inheritParams nT_Param
#' @inheritParams verbose_Param
#'
#' @return Spike regressor matrix or NULL.
#' @keywords internal
parse_scrub_input <- function(scrub, nT, verbose = 1) {
  if (is.null(scrub)) {
    return(NULL)
  }

  if (is.logical(scrub)) {
    if (length(scrub) != nT) {
      stop("Logical vector in `scrub` should indicate which volumes to remove. ",
           "But the length of this vector (", length(scrub), ") does not match ",
           "the number of timepoints (", nT, ").")
    }
    return(logical_to_spike_regressors(scrub, verbose))

  } else if (is.matrix(scrub) || is.data.frame(scrub)) {
    if (is.data.frame(scrub)) {
      scrub <- as.matrix(scrub)
    }

    validate_spike_regressors(scrub, nT)

    if (verbose > 1) {
      cat("Using provided spike regressor matrix: ", ncol(scrub),
          " spike regressors\n")
    }
    return(scrub)

  } else {
    stop("`scrub` must be a logical vector, matrix, data.frame, or NULL")
  }
}

#' Apply scrubbing by adding spike regressors to nuisance matrix
#'
#' This function handles scrubbing (removal of problematic timepoints) by adding
#' spike regressors to the nuisance matrix. The actual "removal" happens during
#' nuisance regression - we don't physically remove timepoints from data.
#'
#' Supports two input formats:
#' 1. Logical vector: TRUE/FALSE indicating which timepoints to scrub
#' 2. Spike regressor matrix: Pre-formatted matrix with one column per scrubbed timepoint
#'
#' @param nuisance Nuisance matrix (timepoints x nuisance_regressors) or NULL.
#' @inheritParams scrub_Param
#' @inheritParams nT_Param
#' @inheritParams verbose_Param
#'
#' @return Updated nuisance matrix with spike regressors appended.
#' @keywords internal
#'
#' @examples
#' # Example 1: Logical vector (scrub timepoints 5, 10, 15)
#' nT <- 100
#' scrub_logical <- rep(FALSE, nT)
#' scrub_logical[c(5, 10, 15)] <- TRUE
#'
#' # Example 2: Spike regressor matrix (3 columns for 3 scrubbed timepoints)
#' scrub_spikes <- matrix(0, nrow = nT, ncol = 3)
#' scrub_spikes[5, 1] <- 1   # First spike at timepoint 5
#' scrub_spikes[10, 2] <- 1  # Second spike at timepoint 10
#' scrub_spikes[15, 3] <- 1  # Third spike at timepoint 15
apply_scrubbing_single <- function(nuisance = NULL, scrub = NULL, nT, verbose = 1) {
  # Parse scrub input into spike regressor matrix
  spike_regressors <- parse_scrub_input(scrub, nT, verbose)

  # If no spike regressors, return nuisance unchanged
  if (is.null(spike_regressors)) return(nuisance)

  # Add spike regressors to nuisance matrix
  if (is.null(nuisance)) {
    updated_nuisance <- spike_regressors
  } else {
    updated_nuisance <- cbind(nuisance, spike_regressors)
  }

  if (verbose > 1) {
    cat("Added ", ncol(spike_regressors), " spike regressors to nuisance matrix\n")
    cat("Final nuisance matrix dimensions: ", paste(dim(updated_nuisance), collapse=" x "), "\n")
  }

  return(updated_nuisance)
}
