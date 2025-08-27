#' Fit Grid-Based HRF Models Across Subjects
#'
#' Tests multiple HRF parameterizations across subjects to find optimal HRF
#' shapes at each brain location. Unlike \code{fit_workingHRF} which uses a
#' single canonical HRF, this function fits a grid of HRF parameter combinations
#' and identifies the best-fitting model for each voxel.
#'
#' @section Processing Pipeline:
#' For each subject:
#' \enumerate{
#'   \item Loads BOLD data
#'   \item Creates design matrices for all HRF parameter combinations (Step 2b)
#'   \item Fits multiGLM comparing all HRF models simultaneously (Step 2c)
#' }
#' Across subjects, extracts best-fitting parameters at each location (Step 2d).
#'
#' @section Computational Requirements:
#' This function is computationally intensive, fitting ~100+ HRF models per subject.
#' Parallel processing is strongly recommended for multiple subjects. Results may
#' be saved to disk in parallel mode to manage memory usage.
#'
#' @inheritParams BOLD_Param
#' @inheritParams EVs_Param
#' @inheritParams TR_Param
#' @inheritParams hrf_grid_Param
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#' @inheritParams hpf_Param
#' @inheritParams nuisance_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams scrub_Param
#' @inheritParams verbose_Param
#' @inheritParams n_cores_Param
#' @param ... Additional arguments passed to \code{hrf_grid} if it's a function.
#'
#' @return Object of class \code{"allHRFs"}; a named list with elements:
#' \describe{
#'   \item{\code{subject_results}}{List of per-subject GLM results from multiGLM()}
#'   \item{\code{best_params_results}}{Data frame with best-fitting HRF parameters
#'         per voxel across all subjects (columns: a1, b1, c, voxel, subject)}
#'   \item{\code{call_info}}{Processing metadata including number of HRF models tested}
#'   \item{\code{hrf_grid}}{The HRF parameter grid used}
#'   \item{\code{session_info}}{Processing settings}
#' }
#'
#' @seealso 
#' \code{\link{fit_workingHRF}} for single canonical HRF modeling,
#' \code{\link{generate_default_hrf_grid}} for creating HRF parameter grids,
#' \code{\link{multiGLM}} for the underlying GLM comparison
#'
#' @examples
#' \dontrun{
#' # Use default HRF grid
#' result <- fit_allHRFs(
#'   BOLD = c("subj1.dtseries.nii", "subj2.dtseries.nii"),
#'   EVs = list(subj1_events, subj2_events),
#'   TR = 0.72,
#'   n_cores = 4
#' )
#' 
#' # Custom HRF grid
#' custom_grid <- data.frame(
#'   a1 = c(4, 6, 8), b1 = c(1, 1, 1), c = c(0, 1/6, 0),
#'   a2 = c(12, 16, 20), b2 = c(1, 1, 1)
#' )
#' 
#' result <- fit_allHRFs(
#'   BOLD = BOLD_files,
#'   EVs = EVs_list,
#'   TR = 0.72,
#'   hrf_grid = custom_grid,
#'   n_cores = 8
#' )
#' }
#'
#' @export
fit_allHRFs <- function(
    BOLD,
    EVs,
    TR,
    hrf_grid = generate_default_hrf_grid,
    brainstructures = c("left", "right"),
    resamp_res = NULL,
    hpf = 0.01,
    nuisance = NULL,
    # Removed derivatives
    onsets = TRUE,
    offsets = TRUE,
    scrub = NULL,
    verbose = 1,
    n_cores = 1,
    ...
) {
  cat("***********Version 1.2222************\n")

  call_match <- match.call()

  hrf_grid <- set_hrf_grid(hrf_grid, ...)

  if (verbose > 0) cat("fit_allHRFs called with", length(BOLD), "subjects and", nrow(hrf_grid), "HRF models using", n_cores, "cores\n")

  # Parallel processing setup for subjects (2b + 2c per subject)
  if(n_cores > 1) {
    subject_results <- run_parallel_subjects_allHRFs(
      BOLD, EVs, nuisance, TR, brainstructures, resamp_res, hpf,
      hrf_grid, onsets, offsets, scrub, verbose, n_cores
    )
  } else {
    if(verbose > 0) cat("Using sequential processing\n")
    subject_results <- run_sequential_subjects_allHRFs(
      BOLD, EVs, nuisance, TR, brainstructures, resamp_res, hpf,
      hrf_grid, onsets, offsets, scrub, verbose
    )
  }

  # Only for parallel mode
  if(n_cores > 1) {
    # Get file paths from the results
    result_paths <- sapply(subject_results, function(x) if(!is.null(x$file_path)) x$file_path else NA)
    # Load the saved .qs files into a list (ignore errors)
    loaded_subject_results <- lapply(result_paths, function(fp)
      if(!is.na(fp) && file.exists(fp)) load_object(file_path = fp, delete_after_load = FALSE) else NULL
    )
  } else {
    # In sequential mode, subject_results is already a list of objects
    loaded_subject_results <- subject_results
  }
  # Needs to be modified for new save/load
  #   cat("Going into report_design_fit_errors")
  #   report_design_fit_errors(subject_results, verbose)

  # Step 2d: Extract best parameters across all subjects (DUMMY for now)
  if(verbose > 0) cat("Extracting best parameters across all subjects (2d step)...\n")
  tictoc::tic()
  best_params_results <- extract_best_params_all_subjects(loaded_subject_results, hrf_grid, verbose)
  tictoc::toc()



  result <- list(
    best_params_results = best_params_results,
    # subject_results = subject_results,
    subject_results = loaded_subject_results,
    call_info = list(
      call = call_match,
      n_subjects = length(BOLD),
      n_hrf_models = nrow(hrf_grid),
      TR = TR,
      n_cores = n_cores,
      parallel_method = if(n_cores > 1) "parLapply" else "sequential",
      completion_time = Sys.time()
    ),
    hrf_grid = hrf_grid,
    session_info = list(
      brainstructures = brainstructures,
      resamp_res = resamp_res,
      hpf = hpf,
      onsets = onsets,
      offsets = offsets
    )
  )

  class(result) <- "allHRFs"
  return(result)
}

#' Process subjects in parallel for allHRFs
#'
#' Sets up parallel cluster for grid-based HRF fitting. Uses save/load system
#' to manage large result objects in memory. Each worker processes one subject
#' through the complete allHRFs pipeline (design matrices + multiGLM).
#'
#' @inheritParams BOLD_Param
#' @inheritParams EVs_Param
#' @inheritParams nuisance_Param
#' @inheritParams TR_Param
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#' @inheritParams hpf_Param
#' @inheritParams hrf_grid_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams scrub_Param
#' @inheritParams verbose_Param
#' @inheritParams n_cores_Param
#'
#' @return List of processing summaries, each containing subject_idx, status 
#'   ("saved" or "error"), and file_path (for successful subjects).
#'
#' @keywords internal
run_parallel_subjects_allHRFs <- function(BOLD, EVs, nuisance, TR, brainstructures, resamp_res,
                                          hpf, hrf_grid, onsets, offsets, scrub, verbose, n_cores) {

  if(verbose > 0) cat("Setting up parallel cluster with", n_cores, "cores\n... Logfile is all_hrf_log.txt.\n")

  n_workers = n_cores
  cat("Processing with cores of:", n_cores)
  cl <- parallel::makeCluster(n_workers, outfile = "all_hrf_logv2.txt")
  on.exit(parallel::stopCluster(cl), add = TRUE)


  vars_to_export <- c(
    "BOLD", "EVs", "nuisance", "scrub", "TR", "brainstructures", "resamp_res",
    "hpf", "hrf_grid", "onsets", "offsets", "verbose"
  )

  setup_parallel_cluster(cl, verbose, vars_to_export)

  if(verbose > 0) cat("Processing subjects in parallel\n")


  subject_results <- parallel::parLapplyLB(cl, 1:length(BOLD), function(i) {
    tryCatch({
      cat(sprintf("###Worker starting subject %d\n", i))
      result <-  process_entire_subject(
        subject_idx = i,
        BOLD_file = BOLD[i],
        EVs = EVs[[i]],
        nuisance_file = if(!is.null(nuisance)) nuisance[i] else NULL,
        scrub = if(!is.null(scrub)) scrub[[i]] else NULL,
        TR = TR, brainstructures = brainstructures, resamp_res = resamp_res,
        hpf = hpf, hrf_grid = hrf_grid, onsets = onsets, offsets = offsets, verbose = verbose
      )


      file_path <- save_object(result, label = sprintf("subject_%03d", i), prefix = "allhrf_", tmp = FALSE)

      cat(sprintf("###Worker finished subject %d successfully\n", i))
      # result
      list(subject_idx = i, status = "saved", file_path = file_path)
    }, error = function(e) {
      cat(sprintf("###Worker FAILED subject %d: %s\n", i, e$message))
      list(subject_idx = i, status = "error", error = e$message)
    })
  })

  # subject_results <- parallel::parLapplyLB(cl, 1:length(BOLD), function(i) {
  #   process_entire_subject(
  #     subject_idx = i,
  #     BOLD_file = BOLD[i],
  #     EVs = EVs[[i]],
  #     nuisance_file = if(!is.null(nuisance)) nuisance[i] else NULL,
  #     scrub = if(!is.null(scrub)) scrub[[i]] else NULL,
  #     TR = TR, brainstructures = brainstructures, resamp_res = resamp_res,
  #     hpf = hpf, hrf_grid = hrf_grid, onsets = onsets, offsets = offsets, verbose = verbose
  #   )
  # })

  if(verbose > 0) cat("Parallel processing completed\n")
  return(subject_results)
}

#' Process subjects sequentially for allHRFs
#'
#' Processes subjects one at a time for grid-based HRF fitting. Used when
#' parallel processing is disabled. Returns full result objects (no save/load).
#'
#' @inheritParams BOLD_Param
#' @inheritParams EVs_Param
#' @inheritParams nuisance_Param
#' @inheritParams TR_Param
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#' @inheritParams hpf_Param
#' @inheritParams hrf_grid_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams scrub_Param
#' @inheritParams verbose_Param
#'
#' @return List of full subject results from \code{process_entire_subject}.
#'
#' @keywords internal
run_sequential_subjects_allHRFs <- function(BOLD, EVs, nuisance, TR, brainstructures, resamp_res,
                                            hpf, hrf_grid, onsets, offsets, scrub, verbose) {

  lapply(1:length(BOLD), function(i) {
    process_entire_subject(
      subject_idx = i,
      BOLD_file = BOLD[i],
      EVs = EVs[[i]],
      nuisance_file = if(!is.null(nuisance)) nuisance[i] else NULL,
      scrub = if(!is.null(scrub)) scrub[[i]] else NULL,
      TR = TR, brainstructures = brainstructures, resamp_res = resamp_res,
      hpf = hpf, hrf_grid = hrf_grid, onsets = onsets, offsets = offsets, verbose = verbose
    )
  })
}

#' Process entire subject through allHRFs pipeline
#'
#' Executes the complete grid-based HRF pipeline for one subject: loads BOLD data,
#' creates design matrices for all HRF parameter combinations (Step 2b), and fits
#' multiGLM comparing all models (Step 2c). Includes extensive memory monitoring
#' and timing diagnostics.
#'
#' @inheritParams subject_idx_Param
#' @param BOLD_file Character. File path to subject's CIFTI data.
#' @inheritParams EVs_Param
#' @inheritParams nuisance_file_Param
#' @inheritParams TR_Param
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#' @inheritParams hpf_Param
#' @inheritParams hrf_grid_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams scrub_Param
#' @inheritParams verbose_Param
#'
#' @return List with elements:
#'   \item{subject_idx}{Subject identifier}
#'   \item{nT}{Number of timepoints}
#'   \item{glm_result}{Complete multiGLM output with bestmodel_xii}
#'   \item{processing_time}{Total time in seconds}
#'   \item{status}{"success" or "error"}
#'   \item{error}{Error message (if status = "error")}
#'
#' @keywords internal
process_entire_subject <- function(subject_idx, BOLD_file, EVs, nuisance_file,
                                   TR, brainstructures, resamp_res, hpf,
                                   hrf_grid, onsets, offsets, scrub, verbose) {

  tictoc::tic()
  cat("***START OF PROCESS_ENTIRE_SUBJECT ZZZ!!!***\n")
  start_time <- Sys.time()
  if(verbose > 1) cat("Subject", subject_idx, ": Starting allHRFs pipeline...\n")

  tryCatch({
    # Debug ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------
    if (!requireNamespace("ps", quietly = TRUE)) stop("Please install the 'ps' package.")
    if (!requireNamespace("pryr", quietly = TRUE)) stop("Please install the 'pryr' package.")

    rss_report <- function(when, subject_idx) {
      rss <- ps::ps_memory_info(ps::ps_handle())[["rss"]]
      cat(sprintf("***Subject %d - %s - RSS: %.2f MB\n", subject_idx, when, rss / 1024^2))
    }
    rss_report("Start", subject_idx)
    mem_before <- pryr::mem_used()
    #------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------

    # Step 1: Load BOLD data (reuse from fit_workingHRF)
    tictoc::tic()
    bold_data <- load_bold_data(BOLD_file, brainstructures, resamp_res)
    cat("*****#%#time elapsed for load_bold_data: ")
    tictoc::toc()
    nT <- bold_data$nT
    if(verbose > 1) cat("Subject", subject_idx, ": BOLD loaded, nT =", bold_data$nT, "\n")

    # Step 2b: Create ALL design matrices for all HRF parameters
    rss_report("***Before design matrix", subject_idx)
    tictoc::tic()
    design_3D <- create_all_design_matrices(
      EVs, nT, TR, hrf_grid, onsets, offsets, verbose, subject_idx
    )
    rss_report("***After design matrix", subject_idx)
    cat("***design_3D RAM size:", format(object.size(design_3D), units = "MB"), "\n")
    cat("****all design matrices took: ")
    tictoc::toc()

    # Step 2c: Fit multiGLM with all designs
    rss_report("***Before GLM fit", subject_idx)
    tictoc::tic()
    glm_result <- fit_multiGLM_all_designs(
      bold_data$BOLD_xii, design_3D, nuisance_file, TR, brainstructures,
      hpf, scrub, resamp_res, verbose, subject_idx
    )
    cat("*****#%#time elapsed for multiGLM: ")
    tictoc::toc()
    rss_report("***After GLM fit", subject_idx)
    cat("***glm_result RAM size:", format(object.size(glm_result), units = "MB"), "\n")

    cat("\n--- MEMORY Z! REPORT for Subject", subject_idx, "---\n")
    cat("Size of BOLD data      :", format(object.size(bold_data), units = "MB"), "\n")
    cat("Size of design_3D      :", format(object.size(design_3D), units = "MB"), "\n")
    cat("Size of glm_result     :", format(object.size(glm_result), units = "MB"), "\n")
    cat("Total mem in function  :", format(object.size(environment()), units = "MB"), "\n")
    cat("--------------------------------------\n")
    cat("line261\n")
    # Clean up memory before returning
    # rm(bold_data, design_3D); gc()

    processing_time <- as.numeric(Sys.time() - start_time, units = "secs")
    if(verbose > 1) cat("Subject", subject_idx, ": Completed in", round(processing_time, 1), "s\n")


    mem_after <- pryr::mem_used()
    cat("line269\n")
    cat(sprintf("***Subject %d - Total mem change: %.2f MB\n", subject_idx,
                as.numeric(mem_after - mem_before) / 1024^2))
    cat("*******ls()***********\n")
    ls()
    cat("^^^***Subject", subject_idx, ": EVERYTHING DONE: ")
    tictoc::toc()
    return(list(
      subject_idx = subject_idx,
      nT = nT,
      glm_result = glm_result,  # Contains bestmodel_xii
      processing_time = processing_time,
      status = "success"
    ))

  }, error = function(e) {
    processing_time <- as.numeric(Sys.time() - start_time, units = "secs")

    return(list(
      subject_idx = subject_idx,
      processing_time = processing_time,
      status = "error",
      error = e$message
    ))
  })
}