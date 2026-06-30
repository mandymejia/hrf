#' Mask out invalid data
#'
#' Mask out data locations that are invalid (missing data, low mean, or low
#'  variance) for any session.
#'
#' @param BOLD A session-length list of \eqn{T \times V} numeric BOLD data.
#' @param meanTol,varTol Tolerance for mean and variance of each data location.
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Defaults: \code{1e-6}.
#' @param verbose Print messages counting how many locations are removed?
#'  Default: \code{TRUE}.
#'
#' @importFrom matrixStats colVars
#' @return A logical vector indicating locations that are valid across all sessions.
#'
#' @examples
#' nT <- 30
#' nV <- 400
#' BOLD1 <- matrix(rnorm(nT*nV), nrow=nT)
#' BOLD1[,seq(30,50)] <- NA
#' BOLD2 <- matrix(rnorm(nT*nV), nrow=nT)
#' BOLD2[,65] <- BOLD2[,65] / 1e10
#' BOLD <- list(sess1=BOLD1, sess2=BOLD2)
#' do_QC(BOLD)
#'
#' @export
do_QC <- function(BOLD, meanTol=1e-6, varTol=1e-6, verbose=TRUE){

  nS <- length(BOLD)
  nV <- ncol(BOLD[[1]])

  # For each BOLD data matrix,
  mask_na <- mask_mean <- mask_var <- mask_snr <- rep(TRUE, nV)
  for (ss in seq(nS)) {
    # Mark columns with any NA or NaN values for removal.
    na_ss <- is.na(BOLD[[ss]]) | is.nan(BOLD[[ss]])
    mask_na[apply(na_ss, 2, any)] <- FALSE
    # Calculate means and variances of columns, except those with any NA or NaN.
    # Mark columns with mean/var falling under the thresholds for removal.
    means_ss <- colMeans(BOLD[[ss]][,mask_na,drop=FALSE])
    vars_ss <- matrixStats::colVars(BOLD[[ss]][,mask_na,drop=FALSE])
    snr_ss <- rep(NA, nV) # use `NA` values for `NA` columns.
    snr_ss[mask_na] <- means_ss/sqrt(vars_ss)
    mask_mean[mask_na][means_ss < meanTol] <- FALSE
    mask_var[mask_na][vars_ss < varTol] <- FALSE
    #mask_snr[mask_na][snr_ss < snrTol] <- FALSE
  }

  # Print counts of locations removed, for each reason.
  if (verbose) {
    warn_part1 <- " locations"
    warn_part2 <- if (nS > 1) { " in at least one session.\n" } else { ".\n" }
    if (any(!mask_na)) {
      cat(paste0(
        "\t", sum(!mask_na), warn_part1,
        " removed due to NA/NaN values", warn_part2
      ))
      warn_part1 <- " additional locations"
    }
    # Do not include NA locations in count.
    mask_mean2 <- mask_mean | (!mask_na)
    if (any(!mask_mean2)) {
      cat(paste0(
        "\t", sum(!mask_mean2), warn_part1,
        " removed due to low mean", warn_part2
      ))
      warn_part1 <- " additional locations"
    }
    # Do not include NA or low-mean locations in count.
    mask_var2 <- mask_var | (!mask_mean) | (!mask_na)
    if (any(!mask_var2)) {
      cat(paste0(
        "\t", sum(!mask_var2), warn_part1,
        " removed due to low variance", warn_part2
      ))
    }
    # # Do not include NA or low-mean or low-var locations in count.
    # mask_snr2 <- mask_snr | (!mask_mean) | (!mask_var) | (!mask_na)
    # if (any(!mask_snr2)) {
    #   cat(paste0(
    #     "\t", sum(!mask_snr2), warn_part1,
    #     " removed due to low SNR", warn_part2
    #   ))
    # }
  }

  # Return composite mask and other masks
  list(
    mask =  mask_na & mask_mean & mask_var, # & mask_snr,
    mask_na = mask_na,
    mask_mean = mask_mean,
    mask_var = mask_var,
    mask_snr = snr_ss #just return the actual SNR values for now
  )
}

#' Is a matrix or data.frame?
#'
#' Is this a matrix or data.frame?
#'
#' @param x The object
#' @return Length-one logical.
#' @keywords internal
is_matrix_or_df <- function(x){
  is.matrix(x) || is.data.frame(x)
}

#' `cbind` if first argument might be \code{NULL}
#'
#' `cbind`, but return the second argument if the first is \code{NULL}
#' @param mat_or_NULL \code{NULL} or a numeric matrix
#' @param to_add A numeric matrix with the same number of rows as \code{mat_or_NULL}
#' @return \code{cbind(mat_or_NULL, to_add)}, or just \code{to_add} if the first argument is NULL.
#' @keywords internal
cbind2 <- function(mat_or_NULL, to_add) {
  if (!is.null(mat_or_NULL)) {
    cbind(mat_or_NULL, to_add)
  } else {
    to_add
  }
}


#' Load BOLD data and extract dimensions
#'
#' Reads CIFTI file using ciftiTools, applies optional resampling, and extracts
#' key dimensional information needed for GLM setup.
#'
#' @param BOLD_file Character. File path to CIFTI data file.
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#' @inheritParams smoothing_Param
#' @inheritParams surf_FWHM_Param
#' @param scale Logical. If \code{TRUE}, convert BOLD to percent signal change
#'   per voxel: \code{(BOLD - mu) / mu * 100} where \code{mu} is the voxel's
#'   temporal mean. Errors out if any voxel has \code{|mu| < 1e-6} (PSC undefined).
#'   Default \code{FALSE}.
#'
#' @return List with elements:
#'   \item{BOLD_xii}{Loaded xifti object containing BOLD time-series}
#'   \item{nT}{Number of timepoints (volumes)}
#'   \item{n_locations}{Total number of brain locations}
#'
#' @keywords internal
load_bold_data <- function(BOLD_file, brainstructures, resamp_res, smoothing = TRUE, surf_FWHM = 5, scale = FALSE) {
  BOLD_xii <- ciftiTools::read_cifti(BOLD_file,
                                     brainstructures = brainstructures,
                                     resamp_res = resamp_res)
  if (isTRUE(smoothing)) {
    cat("Smoothing...\n")
    BOLD_xii <- ciftiTools::smooth_xifti(BOLD_xii, surf_FWHM = surf_FWHM, vol_FWHM = 3)
  }

  if (isTRUE(scale)) {
    mat <- as.matrix(BOLD_xii)        # n_voxels x nT
    mu  <- rowMeans(mat, na.rm = TRUE)
    if (any(!is.finite(mu) | abs(mu) < 1e-6)) {
      stop("Near-zero or non-finite voxel mean detected; PSC scaling invalid.")
    }
    mat <- ((mat - mu) / mu) * 100    # percent signal change
    BOLD_xii <- ciftiTools::newdata_xifti(BOLD_xii, mat)
  }

  nT <- ncol(BOLD_xii)
  n_locations <- nrow(as.matrix(BOLD_xii))

  return(list(
    BOLD_xii = BOLD_xii,
    nT = nT,
    n_locations = n_locations
  ))
}

#' Load nuisance regressors from file
#'
#' Reads nuisance regressor file and converts to matrix format.
#'
#' @inheritParams nuisance_file_Param
#'
#' @return Matrix of nuisance regressors or NULL if file not provided.
#'
#' @keywords internal
load_nuisance_regressors <- function(nuisance_file) {
  if(is.null(nuisance_file)) return(NULL)
  return(as.matrix(read.table(nuisance_file, header=FALSE)))
}

#' Generate time vector for HRF calculation
#'
#' Internal helper to create a time vector (in seconds) extending at least
#' to the specified duration, scaled by TR and upsampling factor.
#'
#' @param TR Numeric. Repetition time in seconds.
#' @param upsample Numeric. Temporal upsampling factor (default = 100).
#' @param duration Numeric. Desired total duration in seconds (default = 40).
#'
#' @return Numeric vector of time points (seconds).
#' @keywords internal
make_inds <- function(TR, upsample = 100, duration = 40) {
  max_index <- ceiling(duration / TR)
  seq(1 / upsample, max_index, 1 / upsample) * TR
}

#' Calculate taper start time for HRF
#'
#' Determines when to start tapering the hemodynamic response function (HRF)
#' based on the location of the undershoot peak. Returns NULL if no tapering
#' is needed (no undershoot or HRF resolves by 30 seconds).
#'
#' @param a1 Shape parameter for the positive gamma function
#' @param b1 Scale parameter for the positive gamma function
#' @param a2 Shape parameter for the negative gamma function
#' @param b2 Scale parameter for the negative gamma function
#' @param c Amplitude of the undershoot (0 = no undershoot)
#' @param TR Time repetition in seconds
#' @param deriv Derivative order (default NULL, passed to HRF_calc)
#'
#' @return Numeric taper start time in seconds, or NULL if no tapering needed
#'
#' @keywords internal
get_taper_start <- function(a1, b1, a2, b2, c, TR, deriv = NULL) {
  if (c == 0) return(NULL)     # Return NULL if no undershoot

  inds <- make_inds(TR)
  hrf_vals <- hrf::HRF_calc(   # Calculate HRF values
    t = inds,
    deriv = deriv,
    a1 = a1,
    b1 = b1,
    a2 = a2,
    b2 = b2,
    c = c
  )
  # Check if HRF at 30 seconds is below threshold
  if (abs(hrf_vals[which.min(abs(inds - 30))]) <= .HRF_THRESHOLD) { return(NULL)}

  peak2_time <- inds[which.min(hrf_vals)]   # Calculate taper_start
  taper_start <- min(peak2_time, 25)

  return(taper_start)
}

#' Convert design matrix to array format required by multiGLM
#'
#' Creates 3D array with noisy duplicate design to satisfy multiGLM requirements.
#'
#' @inheritParams design_matrix_Param
#'
#' @return 3D array (timepoints × regressors × 2)
#'
#' @keywords internal
convert_design_to_array <- function(design_matrix) {
  design_array <- array(design_matrix, dim=c(dim(design_matrix), 2))
  design_array[,,2] <- design_array[,,2] + rnorm(length(design_array[,,2]))
  return(design_array)
}

#' Fit GLM model for one subject
#'
#' Loads nuisance regressors and fits multiGLM to BOLD data using provided
#' design array. Handles filtering and scrubbing as specified.
#'
#' @inheritParams BOLD_xii_Param
#' @inheritParams design_array_Param
#' @inheritParams nuisance_file_Param
#' @inheritParams TR_Param
#' @inheritParams brainstructures_Param
#' @inheritParams hpf_Param
#' @inheritParams scrub_Param
#' @inheritParams verbose_Param
#' @inheritParams subject_idx_Param
#'
#' @return Output from multiGLM() containing GLM results and statistics.
#'
#' @keywords internal
fit_glm_model <- function(BOLD_xii, design_array, nuisance_file, TR, brainstructures,
                          hpf, scrub, verbose, subject_idx) {
  if(verbose > 1) cat("Subject", subject_idx, ": Fitting GLM model...\n")

  # Load nuisance regressors
  nuisance <- load_nuisance_regressors(nuisance_file)

  # Fit single GLM - multiGLM handles whatever design matrix we give it
  glm_result <- multiGLM(
    BOLD = BOLD_xii,
    brainstructures = brainstructures,
    resamp_res = NULL,
    design = design_array,  # Just pass the whole design array
    design_canonical = design_array[,,1],  # Use first design as canonical
    nuisance = nuisance,
    scrub = scrub,
    TR = TR,
    hpf = hpf
  )

  return(glm_result)  # Return single result,
}
