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
#' @inheritParams verbose_Param
#' @inheritParams n_cores_Param
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
#' # Examine activation patterns
#' summary(result$activation_masks$prop)
#' plot(result$activation_masks$prop)  # Requires ciftiTools plotting
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
    verbose = 1,
    n_cores = 1, ...
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
      hrf_params, derivatives, onsets, offsets, scrub, verbose, n_cores
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
  mask_results <- create_activation_masks(subject_results, alpha, verbose)


  result <- list(
    activation_masks = mask_results,
    subject_results = subject_results,
    call_info = list(
      call = call_match,
      n_subjects = length(BOLD),
      TR = TR,
      alpha = alpha,
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
#'
#' @return List of subject results, each containing GLM outputs, design matrices,
#'   processing times, and status information.
#'
#' @keywords internal
run_parallel_subjects_df <- function(BOLD, EVs, nuisance, TR, brainstructures, resamp_res,
                                     hpf, hrf_params, derivatives, onsets, offsets, scrub, verbose, n_cores) {

  if(verbose > 0) cat("Setting up parallel cluster with", n_cores, "cores.\n")

  n_workers = n_cores
  cl <- parallel::makeCluster(n_workers, outfile = "hrf_log_500.txt") # Over subscription
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