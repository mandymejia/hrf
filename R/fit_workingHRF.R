#' Fit a working‑HRF model (multi‑subject, parallel)
#'
#' Fits a general‑linear model (GLM) to BOLD CIFTI time‑series for one or
#' more subjects, using a double‑gamma parameterisation of the
#' haemodynamic response function (HRF). Subjects are processed either
#' sequentially or in parallel with load‑balanced \code{parLapplyLB()}.
#'
#' This function implements the "working HRF" approach, where a canonical
#' (SPM-style) HRF is used to identify brain regions showing significant
#' task-related activation. The results serve as a foundation for subsequent
#' HRF parameter optimization and spatial regularization steps.
#'
#' @section HRF parameterisation:
#' The HRF is modelled as a difference of two gamma functions:
#' \deqn{
#'   h(t) \;=\; c_1 \cdot \text{Gamma}(t; a_1, b_1) \;-\; c_2 \cdot \text{Gamma}(t; a_2, b_2)
#' }
#' 
#' This reproduces SPM's "canonical" HRF when the default parameters are
#' used: \code{a1 = 6}, \code{b1 = 1}, \code{c = 1/6}, \code{a2 = 16}, \code{b2 = 1}.
#'
#' @section Processing Pipeline:
#' For each subject, the function:
#' \enumerate{
#'   \item Loads BOLD data and applies optional resampling
#'   \item Creates design matrix from event definitions using specified HRF
#'   \item Optionally adds HRF derivatives and onset/offset regressors
#'   \item Fits GLM with nuisance regression and high-pass filtering
#'   \item Computes F-statistics for overall task activation
#' }
#' 
#' Across subjects, activation masks and proportion maps are generated
#' to identify consistently active brain regions.
#'
#' @inheritParams BOLD_Param
#' @inheritParams EVs_Param  
#' @inheritParams TR_Param
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#' @inheritParams hpf_Param
#' @inheritParams nuisance_Param
#' @inheritParams hrf_params_Param
#' @inheritParams derivatives_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams scrub_Param
#' @inheritParams alpha_Param
#' @inheritParams min_active_subjects_Param
#' @inheritParams verbose_Param
#' @inheritParams n_cores_Param
#' @inheritParams log_dir_Param
#' @param ... Reserved for future arguments (currently ignored).
#'
#' @return An object of class \code{"workingHRF"}; a named list with elements:
#' \describe{
#'   \item{\code{activation_masks}}{List containing:
#'     \itemize{
#'       \item \code{masks}: Logical matrix (locations × subjects) indicating
#'             significant activation at p < alpha
#'       \item \code{prop}: Numeric vector of activation proportions across subjects
#'       \item \code{alpha}: The threshold used
#'       \item \code{n_subjects}: Number of successfully processed subjects
#'     }
#'   }
#'   \item{\code{subject_results}}{List of per‑subject results, each containing:
#'     \itemize{
#'       \item \code{glm_results}: Output from \code{multiGLM()}
#'       \item \code{design_matrix}: The design matrix used
#'       \item \code{field_names}: Names of regressors
#'       \item \code{processing_time}: Time taken for this subject
#'       \item \code{status}: "success" or "error"
#'     }
#'   }
#'   \item{\code{call_info}}{Processing metadata:
#'     \itemize{
#'       \item \code{call}: Captured function call
#'       \item \code{n_subjects}: Total subjects provided
#'       \item \code{TR}: Repetition time used
#'       \item \code{alpha}: Statistical threshold
#'       \item \code{n_cores}: Cores used
#'       \item \code{parallel_method}: "parLapply" or "sequential"
#'       \item \code{completion_time}: When processing finished
#'     }
#'   }
#'   \item{\code{hrf_params}}{Echo of the supplied HRF parameter list}
#'   \item{\code{session_info}}{Processing settings:
#'     \itemize{
#'       \item \code{brainstructures}: Brain regions analyzed
#'       \item \code{resamp_res}: Resampling resolution
#'       \item \code{hpf}: High-pass filter frequency
#'       \item \code{derivatives}: Whether derivatives were included
#'       \item \code{onsets}, \code{offsets}: Onset/offset modeling flags
#'     }
#'   }
#' }
#'
#' @seealso 
#' \code{\link[ciftiTools]{read_cifti}} for CIFTI data loading,
#' \code{\link{multiGLM}} for the underlying GLM fitting,
#' \code{\link{make_design}} for design matrix creation,
#'
#' @examples
#' \dontrun{
#' # Basic usage with two subjects
#' BOLD_files <- c("subject1.dtseries.nii", "subject2.dtseries.nii")
#' 
#' # Event definitions - one per subject
#' events_subj1 <- data.frame(onset = c(10, 30, 50), duration = c(2, 2, 2))
#' events_subj2 <- data.frame(onset = c(12, 32, 52), duration = c(2, 2, 2))
#' EVs_list <- list(task1 = events_subj1, task1 = events_subj2)
#' 
#' # Fit working HRF model
#' result <- fit_workingHRF(
#'   BOLD = BOLD_files,
#'   EVs = EVs_list,
#'   TR = 0.72,
#'   n_cores = 4,
#'   verbose = 1
#' )
#' 
#' # Access subject-specific results
#' subj1_glm <- result$subject_results[[1]]$glm_results
#' }
#' 
#' @export
fit_workingHRF <- function(
    BOLD,
    EVs,
    TR,
    brainstructures = c("left", "right"),
    resamp_res = NULL,
    hpf = 0.01,
    nuisance = NULL,
    hrf_params = list(a1 = 6, b1 = 1, c = 1/6, a2 = 16, b2 = 1),
    derivatives = TRUE,
    onsets = TRUE,
    offsets = TRUE,
    scrub = NULL,
    alpha = 0.01,
    min_active_subjects = 20,
    verbose = 1,
    n_cores = 1, 
    log_dir = "logs", ...
) {
  call_match <- match.call()
  cat("========Version 3.4444=========\n")
  # Input validation
  validate_inputs(BOLD, EVs, nuisance, hrf_params, n_cores, onsets, offsets, verbose)

  if(verbose > 0) cat("fit_workingHRF called with", length(BOLD), "subjects using", n_cores, "cores\n")

  # Parallel processing setup for making the design matrices and fitting the GLMS
  if(n_cores > 1) {
    subject_results <- run_parallel_subjects_df(
      BOLD, EVs, nuisance, TR, brainstructures, resamp_res, hpf,
      hrf_params, derivatives, onsets, offsets, scrub, verbose, n_cores, log_dir
    )
  } else {
    if(verbose > 0) cat("Using sequential processing\n")
    subject_results <- run_sequential_subjects_df(
      BOLD, EVs, nuisance, TR, brainstructures, resamp_res, hpf,
      hrf_params, derivatives, onsets, offsets, scrub, verbose
    )
  }

  successful_subjects <- subject_results[sapply(subject_results, function(x) x$status == "success")]

  if (length(successful_subjects) == 0) {
    stop("No successful subjects available to create activation masks!")
  }

  report_design_fit_errors(successful_subjects, verbose)

  # Gather P-values for F-test and creat masks
  if(verbose > 0) cat("Creating activation masks and proportion maps...\n")
  mask_results <- create_activation_masks(subject_results, alpha, min_active_subjects, verbose)
  min_active_subjects <- mask_results$min_active_subjects
  mask_results$min_active_subjects <- NULL  # Remove temporary field

  result <- list(
    activation_masks = mask_results,
    subject_results = subject_results,
    call_info = list(
      call = call_match,
      n_subjects = length(BOLD),
      TR = TR,
      alpha = alpha,
      min_active_subjects = min_active_subjects,
      n_cores = n_cores,
      parallel_method = if(n_cores > 1) "parLapply" else "sequential",
      completion_time = Sys.time()
    ),
    hrf_params = hrf_params,
    session_info = list(
      brainstructures = brainstructures,
      resamp_res = resamp_res,
      hpf = hpf,
      derivatives = derivatives,
      onsets = onsets,
      offsets = offsets
    )
  )

  class(result) <- "workingHRF"

  return(result)
}

#' Input validation for fit_workingHRF
#'
#' Validates all input parameters for the main function, checking data types,
#' lengths, ranges, and consistency. Stops execution with informative error
#' messages if validation fails.
#'
#' @inheritParams BOLD_Param
#' @inheritParams EVs_Param
#' @inheritParams nuisance_Param
#' @inheritParams hrf_params_Param
#' @inheritParams n_cores_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams verbose_Param
#'
#' @return TRUE if all validations pass (called for side effects)
#' @keywords internal
validate_inputs <- function(BOLD, EVs, nuisance, hrf_params, n_cores, onsets, offsets, verbose) {
  # Basic input validation
  if(length(BOLD) == 0) stop("BOLD must contain at least one file path or xifti object")
  if(length(EVs) != length(BOLD)) stop("EVs must have the same length as BOLD")
  if(!is.null(nuisance) && length(nuisance) != length(BOLD)) {
    stop("nuisance must have the same length as BOLD or be NULL")
  }

  # Cores validation
  if(!is.numeric(n_cores) || n_cores < 1 || n_cores != round(n_cores)) {
    stop("n_cores must be a positive integer")
  }

  # Safely check available cores (may fail on some HPC systems)
  max_cores <- tryCatch({
    cores <- parallel::detectCores(logical = FALSE)
    ifelse(is.numeric(cores) && cores > 0, cores, NULL)
  }, error = function(e) NULL)

  if(!is.null(max_cores) && n_cores > max_cores) {
    warning("n_cores (", n_cores, ") exceeds available physical cores (", max_cores,
            "). Performance may degrade due to oversubscription.")
  }

  # HRF parameters validation
  required_params <- c("a1", "b1", "c")
  if(!all(required_params %in% names(hrf_params))) {
    stop("hrf_params must contain at least: ", paste(required_params, collapse=", "))
  }

  if(hrf_params$a1 <= 0) stop("hrf_params$a1 must be positive")
  if(hrf_params$b1 <= 0) stop("hrf_params$b1 must be positive")
  if(hrf_params$c < 0) stop("hrf_params$c must be non-negative")
  
  stopifnot(is.logical(onsets), length(onsets) == 1)
  stopifnot(is.logical(offsets), length(offsets) == 1)


  return(TRUE)
}

#' Process subjects in parallel using cluster workers
#'
#' Sets up a parallel cluster and processes all subjects simultaneously using
#' load-balanced \code{parLapplyLB()}.
#'
#' @inheritParams BOLD_Param
#' @inheritParams EVs_Param
#' @inheritParams nuisance_Param
#' @inheritParams TR_Param
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#' @inheritParams hpf_Param
#' @inheritParams hrf_params_Param
#' @inheritParams derivatives_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams scrub_Param
#' @inheritParams verbose_Param
#' @inheritParams n_cores_Param
#' @inheritParams log_dir_Param
#'
#' @return List of subject results, each containing GLM outputs, design matrices,
#'   processing times, and status information.
#'
#' @keywords internal
run_parallel_subjects_df <- function(BOLD, EVs, nuisance, TR, brainstructures, resamp_res,
                                     hpf, hrf_params, derivatives, onsets, offsets, scrub, verbose, n_cores, log_dir) {

  if(verbose > 0) cat("Setting up parallel cluster with", n_cores, "cores.\n")

  n_workers = n_cores

  if (!dir.exists(log_dir)) { dir.create(log_dir, recursive = TRUE)}
  log_file <- file.path(log_dir, sprintf("fit_workingHRF_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
  message("Cluster logs will be saved to: ", normalizePath(log_file, mustWork = FALSE))

  cl <- parallel::makeCluster(n_workers, outfile = log_file)

  on.exit(parallel::stopCluster(cl), add = TRUE)

  vars_to_export <- c("BOLD", "EVs", "nuisance", "scrub", "TR", "brainstructures", "resamp_res",
                      "hpf", "hrf_params", "derivatives", "onsets", "offsets", "verbose")
  setup_parallel_cluster(cl, verbose, vars_to_export)

  if(verbose > 0) cat("Processing subjects in parallel\n")

  subject_results <- parallel::parLapplyLB(cl, 1:length(BOLD), function(i) {
    process_single_subject_df(
      subject_idx = i,
      BOLD_file = BOLD[i],
      EVs = EVs[[i]],
      nuisance_file = if(!is.null(nuisance)) nuisance[i] else NULL,
      scrub = if(!is.null(scrub)) scrub[[i]] else NULL,
      TR = TR, brainstructures = brainstructures, resamp_res = resamp_res,
      hpf = hpf, hrf_params = hrf_params, derivatives = derivatives,
      onsets = onsets, offsets = offsets, verbose = verbose
    )
  })

  if(verbose > 0) cat("Parallel processing completed\n")
  return(subject_results)
}

#' Process subjects sequentially (single-threaded)
#'
#' Processes all subjects one at a time using \code{lapply()}. Used when
#' \code{n_cores = 1} or when parallel processing is not desired.
#' Simpler than parallel version but slower for multiple subjects.
#'
#' @inheritParams BOLD_Param
#' @inheritParams EVs_Param
#' @inheritParams nuisance_Param
#' @inheritParams TR_Param
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#' @inheritParams hpf_Param
#' @inheritParams hrf_params_Param
#' @inheritParams derivatives_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams scrub_Param
#' @inheritParams verbose_Param
#'
#' @return List of subject results, identical structure to parallel version.
#'
#' @keywords internal
run_sequential_subjects_df <- function(BOLD, EVs, nuisance, TR, brainstructures, resamp_res,
                                       hpf, hrf_params, derivatives, onsets, offsets, scrub, verbose) {

  lapply(1:length(BOLD), function(i) {
    process_single_subject_df(
      subject_idx = i,
      BOLD_file = BOLD[i],
      EVs = EVs[[i]],
      nuisance_file = if(!is.null(nuisance)) nuisance[i] else NULL,
      scrub = if(!is.null(scrub)) scrub[[i]] else NULL,
      TR = TR, brainstructures = brainstructures, resamp_res = resamp_res,
      hpf = hpf, hrf_params = hrf_params, derivatives = derivatives,
      onsets = onsets, offsets = offsets, verbose = verbose
    )
  })
}

#' Process one subject through complete pipeline
#'
#' Executes the full working-HRF pipeline for a single subject:
#' loads BOLD data, creates design matrix, fits GLM, and returns results.
#' Includes error handling and timing information. Called by both parallel
#' and sequential processing functions.
#'
#' @inheritParams subject_idx_Param
#' @param BOLD_file Character. File path to subject's CIFTI data.
#' @inheritParams EVs_Param
#' @param nuisance_file Character or NULL. File path to nuisance regressors.
#' @inheritParams TR_Param
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#' @inheritParams hpf_Param
#' @inheritParams hrf_params_Param
#' @inheritParams derivatives_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams scrub_Param
#' @inheritParams verbose_Param
#'
#' @return List with elements:
#'   \item{subject_idx}{Subject identifier}
#'   \item{n_locations}{Number of brain locations}
#'   \item{nT}{Number of timepoints}
#'   \item{design_dims}{Dimensions of design matrix}
#'   \item{design_matrix}{The design matrix used}
#'   \item{field_names}{Names of design matrix columns}
#'   \item{glm_results}{Output from multiGLM()}
#'   \item{processing_time}{Time taken in seconds}
#'   \item{status}{"success" or "error"}
#'   \item{error}{Error message (only if status = "error")}
#'
#' @keywords internal
process_single_subject_df <- function(subject_idx, BOLD_file, EVs, nuisance_file,
                                      TR, brainstructures, resamp_res, hpf,
                                      hrf_params, derivatives, onsets, offsets, scrub, verbose) {
  start_time <- Sys.time()
  if(verbose > 1) cat("Subject", subject_idx, ": Starting pipeline...\n")

  tryCatch({
    tictoc::tic()
    # Step 1a: Load data and create design matrix
    bold_and_design <- create_design_matrix(
      BOLD_file, EVs, TR, brainstructures, resamp_res, hrf_params, derivatives,
      onsets, offsets, verbose, subject_idx
    )
    cat("***Subject", subject_idx, ": Design matrix created in: ")
    tictoc::toc()

    # Step 1b: Fit GLM models
    tictoc::tic()
    glm_results <- fit_glm_model(
      bold_and_design$BOLD_xii, bold_and_design$design_array,
      nuisance_file, TR, brainstructures, hpf, scrub, verbose, subject_idx
    )
    cat("***Subject", subject_idx, ": GLM fitted in: ")
    tictoc::toc()

    bold_and_design$BOLD_xii <- NULL
    bold_and_design$design_array <- NULL
    gc()

    processing_time <- as.numeric(Sys.time() - start_time, units = "secs")
    if(verbose > 1) cat("Subject", subject_idx, ": Completed in", round(processing_time, 1), "s\n")

    return(list(
      subject_idx = subject_idx,
      n_locations = bold_and_design$n_locations,
      nT = bold_and_design$nT,
      design_dims = bold_and_design$design_dims,
      design_matrix = bold_and_design$design_matrix,
      field_names = bold_and_design$field_names,
      glm_results = glm_results,
      processing_time = processing_time,
      status = "success"
    ))

  }, error = function(e) {
    cat("Error processing subject", subject_idx, "\n")
    cat("Error message:", conditionMessage(e), "\n")
    processing_time <- as.numeric(Sys.time() - start_time, units = "secs")

    return(list(
      subject_idx = subject_idx,
      processing_time = processing_time,
      status = "error",
      error = e$message
    ))
  })
}

#' Create design matrix for one subject
#'
#' Loads BOLD data and creates GLM design matrix by convolving event definitions
#' with specified HRF. Handles derivatives, onset/offset modeling, and converts
#' to array format required by multiGLM().
#'
#' @param BOLD_file Character. File path to subject's CIFTI data.
#' @inheritParams EVs_Param
#' @inheritParams TR_Param
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#' @inheritParams hrf_params_Param
#' @inheritParams derivatives_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams verbose_Param
#' @inheritParams subject_idx_Param
#'
#' @return List with elements:
#'   \item{BOLD_xii}{Loaded xifti object}
#'   \item{n_locations}{Number of brain locations}
#'   \item{nT}{Number of timepoints}
#'   \item{design_matrix}{Design matrix (timepoints × regressors)}
#'   \item{design_array}{3D array version for multiGLM()}
#'   \item{field_names}{Names of design matrix columns}
#'   \item{design_dims}{Dimensions of design matrix}
#'
#' @keywords internal
create_design_matrix <- function(BOLD_file, EVs, TR, brainstructures, resamp_res,
                                 hrf_params, derivatives, onsets, offsets, verbose, subject_idx) {
  if(verbose > 1) cat("Subject", subject_idx, ": Loading BOLD and creating design matrix...\n")

  tictoc::tic()
  bold_data <- load_bold_data(BOLD_file, brainstructures, resamp_res)
  cat("***Subject", subject_idx, ": BOLD data loaded in: ")
  tictoc::toc()
  nT <- bold_data$nT
  n_locations <- bold_data$n_locations

  # Create design matrix
  dHRF <- if(derivatives) 2 else 0

  taper_start <- get_taper_start(
    a1 = hrf_params$a1,
    b1 = hrf_params$b1,
    a2 = hrf_params$a2,
    b2 = hrf_params$b2,
    c = hrf_params$c,
    TR = TR,
    deriv = 0 # We don't use derivative for when DECIDING to taper
  )
 
  cat(if(is.null(taper_start)) "No tapering needed (c=0 or HRF resolves by 30s)\n" else sprintf("Taper start: %.2f seconds\n", taper_start))

  tictoc::tic("***Made design matrix")
  design_result <- make_design(
    EVs = EVs, nTime = nT, TR = TR, dHRF = dHRF,
    onset = onsets,
    offset = offsets,
    taper_start = taper_start,
    a1 = hrf_params$a1, b1 = hrf_params$b1, c = hrf_params$c,
    a2 = hrf_params$a2, b2 = hrf_params$b2
  )
  tictoc::toc()

  # Convert to array format for multiGLM (research code pattern)
  design_array <- convert_design_to_array(design_result$design)

  if(verbose > 1) {
    cat("Subject", subject_idx, ": Design dims =", paste(dim(design_result$design), collapse="x"), "\n")
    cat("Subject", subject_idx, ": Field names =", paste(design_result$field_names, collapse=", "), "\n")
  }

  return(list(
    BOLD_xii = bold_data$BOLD_xii,
    n_locations = n_locations,
    nT = nT,
    design_matrix = design_result$design,
    design_array = design_array,
    field_names = design_result$field_names,
    design_dims = dim(design_result$design)
  ))
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

#' Load BOLD data and extract dimensions
#'
#' Reads CIFTI file using ciftiTools, applies optional resampling, and extracts
#' key dimensional information needed for GLM setup.
#'
#' @param BOLD_file Character. File path to CIFTI data file.
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#'
#' @return List with elements:
#'   \item{BOLD_xii}{Loaded xifti object containing BOLD time-series}
#'   \item{nT}{Number of timepoints (volumes)}
#'   \item{n_locations}{Total number of brain locations}
#'
#' @keywords internal
load_bold_data <- function(BOLD_file, brainstructures, resamp_res) {
  BOLD_xii <- ciftiTools::read_cifti(BOLD_file,
                                     brainstructures = brainstructures,
                                     resamp_res = resamp_res)
  nT <- ncol(BOLD_xii)
  n_locations <- nrow(as.matrix(BOLD_xii))

  return(list(
    BOLD_xii = BOLD_xii,
    nT = nT,
    n_locations = n_locations
  ))
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


#' Report subject processing errors
#'
#' Prints error messages for failed subjects and issues warnings. Only reports
#' errors without stopping execution - failed subjects are excluded from
#' subsequent cross-subject analyses.
#'
#' @inheritParams subject_results_Param  
#' @inheritParams verbose_Param
#'
#' @return NULL (called for side effects only)
#' @keywords internal
report_design_fit_errors <- function(subject_results, verbose = 1) {
  # Fatal errors
  errors <- sapply(subject_results, function(x) x$status == "error")
  if (any(errors)) {
    error_subjects <- which(errors)
    for (i in error_subjects) {
      cat("Subject", i, "ERROR:", subject_results[[i]]$error, "\n")
    }
    warning("Processing failed for ", sum(errors), " subjects. These subjects will not be
            used when making the activation masks")
  }

}

#' Create activation masks and proportion maps
#'
#' Generates cross-subject activation masks by thresholding F-test p-values
#' and computing proportion maps showing consistency of activation across subjects.
#' This implements the "1c functionality" from the research pipeline.
#'
#' @inheritParams subject_results_Param
#' @inheritParams alpha_Param
#' @param min_active_subjects Integer. Minimum number of subjects required for a voxel 
#'   to be considered active in the group mask
#' @inheritParams verbose_Param
#'
#' @return List with elements:
#'   \item{masks}{Logical matrix (locations × subjects) of significant activations}
#'   \item{mask_prop_NA}{Logical vector with NA for inactive voxels, TRUE for active}
#'   \item{prop}{Numeric vector of activation proportions across subjects}
#'   \item{alpha}{The threshold used}
#'   \item{n_subjects}{Number of successful subjects}
#'   \item{min_active_subjects}{Integer. Adjusted minimum number of active subjects 
#'     (capped at number of successful subjects if necessary)}
#'
#' @keywords internal
create_activation_masks <- function(subject_results, alpha = 0.01, min_active_subjects, verbose = 1) {
  successful_subjects <- which(sapply(subject_results, function(x) x$status == "success"))
  n_successful <- length(successful_subjects)

  if(n_successful == 0) {
    stop("No subjects processed successfully - cannot create activation masks")
  }

  # Get dimensions from first subject
  first_subject <- subject_results[[successful_subjects[1]]]
  first_pvals <- as.vector(as.matrix(first_subject$glm_results$pvalF_xii))
  n_locations <- length(first_pvals)

  # Collect p-values across all successful subjects
  pvals_matrix <- matrix(NA, nrow = n_locations, ncol = n_successful)

  for(i in seq_along(successful_subjects)) {
    subj_idx <- successful_subjects[i]
    pvals <- as.vector(as.matrix(subject_results[[subj_idx]]$glm_results$pvalF_xii))
    pvals_matrix[, i] <- pvals
  }

  # Create masks and proportions
  masks <- (pvals_matrix < alpha)
  colnames(masks) <- successful_subjects
  prop <- rowSums(masks, na.rm = TRUE) / rowSums(!is.na(masks))

  # Create mask_prop_NA
  masks_df <- masks_df <- as.data.frame(masks) # Reshape dataframe
  masks_df$voxel <- 1:nrow(masks_df)
  masks_df <- reshape2::melt(masks_df, id.vars = 'voxel', variable.name = 'subject', value.name = 'mask')

  masks_df_prop <- masks_df %>% group_by(voxel) %>% summarize(prop = mean(mask, na.rm=TRUE))
  
  subj_prop_result <- get_subj_prop(min_active_subjects, subject_results)
  subj_prop <- subj_prop_result$subj_prop
  min_active_subjects <- subj_prop_result$min_active_subjects


  mask_prop <- (masks_df_prop$prop > subj_prop)   # mask_prop <- (masks_df_prop$prop > 0.1)
  
  mask_prop_NA <- mask_prop; mask_prop_NA[!mask_prop_NA] <- NA

  # Check if mask has any active voxels
  n_active_voxels <- sum(mask_prop, na.rm = TRUE)
  total_voxels <- length(mask_prop)
  
  if (n_active_voxels == 0) {
    stop("No active voxels found in group mask! Try lowering alpha or min_active_subjects.")
  }
  
  if (n_active_voxels < total_voxels / 2) {
    warning("Less than half of voxels (", n_active_voxels, "/", total_voxels, 
            ") are active in the group mask. Consider adjusting alpha or min_active_subjects.")
  }


  return(list(
    masks = masks,
    mask_prop_NA = mask_prop_NA,
    prop = prop,
    alpha = alpha,
    n_subjects = n_successful,
    min_active_subjects = min_active_subjects
  ))
}

#' Calculate the Proportion of Active Subjects
#'
#' Computes the proportion threshold for determining active voxels based on
#' the minimum number of active subjects required.
#'
#' @param min_active_subjects Integer. Minimum number of subjects required 
#'   for a voxel to be considered active
#' @param subject_results List of subject results from HRF analysis
#' @return List with elements:
#'   \item{subj_prop}{Numeric proportion of active subjects (rounded to two decimals)}
#'   \item{min_active_subjects}{Integer. Adjusted minimum number of active subjects 
#'     (capped at number of successful subjects if necessary)}
#' @keywords internal
get_subj_prop <- function(min_active_subjects, subject_results) {
  # Count successful subjects
  n_successful <- sum(vapply(subject_results, \(x) identical(x[["status"]], "success"), logical(1)))
  stopifnot("No subjects processed successfully - cannot create activation masks" = n_successful > 0)

  calculated_min <- min(min_active_subjects, floor(n_successful * 0.80))
  
  if (min_active_subjects > calculated_min) {
    warning("min_active_subjects (", min_active_subjects, ") exceeds 80% of successful subjects (", 
            n_successful, "). Using ", calculated_min, " instead.")
    min_active_subjects <- calculated_min
  }
  subj_prop <- round(min_active_subjects / n_successful, 2)
  cat(sprintf("Using group mask threshold: %s (%d out of %d successful subjects)\n", 
            subj_prop, min_active_subjects, n_successful))
  return(list(subj_prop = subj_prop, min_active_subjects = min_active_subjects))
}
