#' Fit working HRF, grid of candidate HRFs, and regularize across subjects
#'
#' For each subject: fits a canonical working HRF (used to build the cross-
#' subject activation mask via F-test) and a grid of candidate HRF parameter
#' combinations (one best HRF per voxel). Across subjects: regularizes the
#' candidate results into a population-average HRF parameter map.
#'
#' @inheritParams BOLD_Param
#' @inheritParams EVs_Param
#' @inheritParams TR_Param
#' @inheritParams hrf_grid_Param
#' @inheritParams brainstructures_Param
#' @inheritParams resamp_res_Param
#' @inheritParams hpf_Param
#' @inheritParams nuisance_Param
#' @param working_hrf Named list of HRF params for the canonical working fit:
#'   \code{a1, b1, c, a2, b2}. Default is SPM canonical
#'   \code{list(a1=6, b1=1, c=1/6, a2=16, b2=1)}.
#' @param derivatives Logical. If \code{TRUE}, include temporal + dispersion
#'   derivatives in the working-HRF design (used for the F-test mask).
#' @param alpha Numeric. Per-subject p-value threshold (Bonferroni-corrected)
#'   for the working-HRF F-test that builds the activation mask.
#' @param min_active_subjects Integer. Minimum number of subjects that must
#'   pass the F-test at a voxel for it to enter the group activation mask.
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams scrub_Param
#' @inheritParams smoothing_Param
#' @inheritParams surf_FWHM_Param
#' @inheritParams verbose_Param
#' @inheritParams n_cores_Param
#' @inheritParams log_dir_Param
#' @inheritParams work_dir_Param
#' @param save_rss Logical. If \code{TRUE}, store per-voxel RSS for every HRF
#'   model in each subject's saved \code{.qs} so downstream callers can use
#'   lookup mode. If \code{FALSE} (default), the RSS is dropped on save.
#' @param ... Additional arguments passed to \code{hrf_grid} if it's a function.
#'
#' @return Object of class \code{"hrfs"}; a named list with three sub-objects
#'   (each carrying its own S3 class for plot dispatch) plus shared metadata:
#' \describe{
#'   \item{\code{fit_workingHRF}}{class \code{"workingHRF"}: activation_masks,
#'     per-subject working GLMs, hrf_params.}
#'   \item{\code{fit_allHRFs}}{class \code{"allHRFs"}: per-subject 25-candidate
#'     GLM cache (path attr "result_paths"), best_params_results, hrf_grid.}
#'   \item{\code{regularize_allHRFs}}{class \code{"regularizeHRFs"}: pop_avg,
#'     best_params_df, winning_c, c_votes, mask_prop_NA, hrf_grid.}
#'   \item{\code{session_info}}{Pipeline knobs used (brainstructures, resamp_res,
#'     hpf, onsets, offsets, smoothing, surf_FWHM, derivatives, alpha,
#'     min_active_subjects).}
#' }
#'   \code{attr(combo, "call_info")} holds the run metadata.
#'
#' @seealso
#' \code{\link{generate_hrf_grid}} for creating HRF parameter grids,
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
    hrf_grid = generate_hrf_grid,
    brainstructures = c("left", "right"),
    resamp_res = NULL,
    hpf = 0.01,
    nuisance = NULL,
    working_hrf = list(a1 = 6, b1 = 1, c = 1/6, a2 = 16, b2 = 1),
    derivatives = TRUE,
    alpha = 0.05,
    min_active_subjects = 20,
    onsets = TRUE,
    offsets = TRUE,
    scrub = NULL,
    smoothing = TRUE,
    surf_FWHM = 5,
    verbose = 1,
    n_cores = 1,
    save_rss = FALSE,
    log_dir = "logs", work_dir = NULL, ...
) {
  call_match <- match.call()

  hrf_grid <- set_hrf_grid(hrf_grid, ...)

  if (is.null(work_dir)) { work_dir <- "work" }
  work_dir <- file.path(work_dir, paste0("work_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample(1000:9999, 1)))

  if (verbose > 0) cat("fit_allHRFs called with", length(BOLD), "subjects and", nrow(hrf_grid), "HRF models using", n_cores, "cores\n")

  cfg <- list(
    BOLD = BOLD, EVs = EVs, nuisance = nuisance, scrub = scrub,
    TR = TR, brainstructures = brainstructures, resamp_res = resamp_res,
    hpf = hpf, smoothing = smoothing, surf_FWHM = surf_FWHM,
    onsets = onsets, offsets = offsets,
    hrf_grid = hrf_grid,
    working_hrf = working_hrf, derivatives = derivatives,
    alpha = alpha, min_active_subjects = min_active_subjects,
    verbose = verbose
  )

  # Parallel processing setup for subjects (2b + 2c per subject)
  if (n_cores > 1) {
    subject_results <- run_parallel_subjects_allHRFs(cfg, n_cores = n_cores,
                                                      log_dir = log_dir,
                                                      work_dir = work_dir,
                                                      save_rss = save_rss)
  } else {
    if (verbose > 0) cat("Using sequential processing\n")
    subject_results <- run_sequential_subjects_allHRFs(cfg, work_dir = work_dir,
                                                        save_rss = save_rss)
  }

 
  result_paths <- sapply(subject_results, function(x) if(!is.null(x$file_path)) x$file_path else NA)
  output_root <- getOption("hrf.output_root", default = NULL)
  result_paths_stored <- if (!is.null(output_root)) {
    sub(paste0("^", normalizePath(output_root, mustWork = FALSE), "/?"), "", result_paths)
  } else {
    result_paths
  }
  # Load the saved .qs files into a list and strip out design_3D to save memory
  loaded_subject_results <- lapply(result_paths, function(fp) {
    if(!is.na(fp) && file.exists(fp)) {
      obj <- load_object(file_path = fp, delete_after_load = FALSE)
      obj$design_3D <- NULL  # Remove the large 3D array from memory
      obj$glm_result$mGLM0s$cortexL$RSS <- NULL
      obj$glm_result$mGLM0s$cortexR$RSS <- NULL
      return(obj)
    } else {
      return(NULL)
    }
  })

  if(verbose > 0) cat("Extracting best parameters across all subjects (2d step)...\n")
  tictoc::tic()
  best_params_results <- extract_best_params_all_subjects(loaded_subject_results, hrf_grid, verbose)
  tictoc::toc()


  if (verbose > 0) cat("Aggregating working-HRF subject results into activation masks...\n")
  working_subject_results <- lapply(loaded_subject_results, `[[`, "working")
  mask_results <- create_activation_masks(working_subject_results, alpha, min_active_subjects, verbose)
  min_active_subjects <- mask_results$min_active_subjects
  mask_results$min_active_subjects <- NULL

  fit_workingHRF <- list(
    activation_masks = mask_results,
    subject_results  = working_subject_results,
    hrf_params       = working_hrf
  )
  class(fit_workingHRF) <- "workingHRF"

  fit_allHRFs <- list(
    best_params_results = best_params_results,
    subject_results     = loaded_subject_results,
    hrf_grid            = hrf_grid
  )
  attr(fit_allHRFs, "result_paths") <- result_paths_stored
  class(fit_allHRFs) <- "allHRFs"

  if (verbose > 0) cat("Running regularize_allHRFs aggregation...\n")
  reg_result <- regularize_allHRFs(fit_workingHRF, fit_allHRFs, verbose = verbose)

  combo <- list(
    fit_workingHRF     = fit_workingHRF,
    fit_allHRFs        = fit_allHRFs,
    regularize_allHRFs = reg_result,
    session_info = list(
      brainstructures = brainstructures, resamp_res = resamp_res, hpf = hpf,
      onsets = onsets, offsets = offsets,
      smoothing = smoothing, surf_FWHM = surf_FWHM,
      derivatives = derivatives,
      alpha = alpha, min_active_subjects = min_active_subjects
    )
  )
  attr(combo, "call_info") <- list(
    call = call_match,
    n_subjects = length(BOLD),
    n_hrf_models = nrow(hrf_grid),
    TR = TR,
    n_cores = n_cores,
    parallel_method = if (n_cores > 1) "parLapply" else "sequential",
    completion_time = Sys.time(),
    save_rss = save_rss
  )
  class(combo) <- "hrfs"
  combo
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
#' @inheritParams log_dir_Param
#' @inheritParams work_dir_Param
#'
#' @return List of processing summaries, each containing subject_idx, status 
#'   ("saved" or "error"), and file_path (for successful subjects).
#'
#' @keywords internal
run_parallel_subjects_allHRFs <- function(cfg, n_cores, log_dir, work_dir, save_rss) {

  if (cfg$verbose > 0) cat("Setting up parallel cluster with", n_cores, "cores\n")

  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
  log_file <- file.path(log_dir, sprintf("fit_allHRFs_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
  message("Cluster logs will be saved to: ", normalizePath(log_file, mustWork = FALSE))

  cl <- parallel::makeCluster(n_cores, outfile = log_file)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  setup_parallel_cluster(cl, cfg$verbose, c("cfg", "save_rss"))

  if (cfg$verbose > 0) cat("Processing subjects in parallel\n")

  parallel::parLapplyLB(cl, seq_along(cfg$BOLD), function(i) {
    tryCatch({
      cat(sprintf("[cluster] worker started subj %d\n", i))
      result <- process_entire_subject(subject_idx = i, cfg = cfg, save_rss = save_rss)
      file_path <- save_object(result, label = sprintf("subject_%03d", i),
                               prefix = "allhrf_", tmp = FALSE, work_dir = work_dir)
      cat(sprintf("[cluster] worker done    subj %d\n", i))
      list(subject_idx = i, status = "saved", file_path = file_path)
    }, error = function(e) {
      cat(sprintf("[cluster] worker FAILED  subj %d: %s\n", i, e$message))
      list(subject_idx = i, status = "error", error = e$message)
    })
  })
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
#' @inheritParams work_dir_Param
#'
#' @return List of full subject results from \code{process_entire_subject}.
#'
#' @keywords internal
run_sequential_subjects_allHRFs <- function(cfg, work_dir, save_rss) {
  lapply(seq_along(cfg$BOLD), function(i) {
    tryCatch({
      result <- process_entire_subject(subject_idx = i, cfg = cfg, save_rss = save_rss)
      file_path <- save_object(result, label = sprintf("subject_%03d", i),
                               prefix = "allhrf_", tmp = FALSE, work_dir = work_dir)
      list(subject_idx = i, status = "saved", file_path = file_path)
    }, error = function(e) {
      list(subject_idx = i, status = "error", error = e$message)
    })
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
#' @param cfg Named list of pipeline config (built by \code{fit_allHRFs}).
#' @param save_rss Logical. If \code{TRUE}, keep per-candidate RSS arrays in the
#'   saved \code{.qs} (needed for downstream lookup-mode scoring).
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
process_entire_subject <- function(subject_idx, cfg, save_rss = FALSE) {

  tictoc::tic()
  cat(sprintf("[subj %d] start\n", subject_idx))
  start_time <- Sys.time()

  BOLD_file     <- cfg$BOLD[subject_idx]
  EVs           <- cfg$EVs[[subject_idx]]
  nuisance_file <- if (!is.null(cfg$nuisance)) cfg$nuisance[subject_idx] else NULL
  scrub         <- if (!is.null(cfg$scrub))    cfg$scrub[[subject_idx]]  else NULL

  tryCatch({
    tictoc::tic()
    bold_data <- load_bold_data(BOLD_file, cfg$brainstructures, cfg$resamp_res, cfg$smoothing, cfg$surf_FWHM)
    cat(sprintf("[subj %d] load_bold_data: ", subject_idx)); tictoc::toc()
    nT <- bold_data$nT

    tictoc::tic()
    working_result <- fit_working_one_subject(
      BOLD_xii      = bold_data$BOLD_xii,
      nuisance_file = nuisance_file,
      scrub         = scrub,
      cfg           = cfg,
      subject_idx   = subject_idx
    )
    cat(sprintf("[subj %d] fit_working_one_subject: ", subject_idx)); tictoc::toc()

    tictoc::tic()
    design_3D <- create_all_design_matrices(EVs, nT, cfg$TR, cfg$hrf_grid,
                                            cfg$onsets, cfg$offsets, cfg$verbose, subject_idx)
    cat(sprintf("[subj %d] all design matrices: ", subject_idx)); tictoc::toc()

    tictoc::tic()
    glm_result <- fit_multiGLM_all_designs(
      bold_data$BOLD_xii, design_3D$array, nuisance_file, cfg$TR, cfg$brainstructures,
      cfg$hpf, scrub, cfg$resamp_res, cfg$verbose, subject_idx, EVs, cfg$onsets, cfg$offsets
    )
    cat(sprintf("[subj %d] multiGLM: ", subject_idx)); tictoc::toc()

    # Extract best (a1, b1, rss) per voxel for EACH unique c value
    # This reduces storage from ~17 MB/subject to ~0.9 MB/subject
    unique_c <- unique(cfg$hrf_grid$c)

    for (hemi in c("cortexL", "cortexR")) {
      RSS <- glm_result$mGLM0s[[hemi]]$RSS

      if (!is.null(RSS)) {
        best_per_c <- list()

        for (c_val in unique_c) {
          c_idx <- which(cfg$hrf_grid$c == c_val)
          RSS_c <- RSS[, c_idx]
          best_local <- apply(RSS_c, 1, which.min)
          best_global <- c_idx[best_local]

          best_per_c[[as.character(c_val)]] <- data.frame(
            a1 = cfg$hrf_grid$a1[best_global],
            b1 = cfg$hrf_grid$b1[best_global],
            c = c_val,
            rss = RSS_c[cbind(1:nrow(RSS), best_local)]
          )
        }

        glm_result$mGLM0s[[hemi]]$best_per_c <- best_per_c
        if (!save_rss) {
          glm_result$mGLM0s[[hemi]]$RSS <- NULL
        }
      }
    }

    cat(sprintf(
      "[subj %d] memory: bold=%s, design=%s, glm=%s, env=%s\n",
      subject_idx,
      format(object.size(bold_data),    units = "MB"),
      format(object.size(design_3D),    units = "MB"),
      format(object.size(glm_result),   units = "MB"),
      format(object.size(environment()), units = "MB")
    ))

    processing_time <- as.numeric(Sys.time() - start_time, units = "secs")
    cat(sprintf("[subj %d] done: ", subject_idx)); tictoc::toc()
    return(list(
      subject_idx = subject_idx,
      nT = nT,
      design_3D = design_3D,
      glm_result = glm_result,  # Contains bestmodel_xii
      working = working_result,
      processing_time = processing_time,
      status = "success"
    ))

  }, error = function(e) {
    processing_time <- as.numeric(Sys.time() - start_time, units = "secs")
    cat(sprintf("[subj %d] ERROR: %s\n", subject_idx, conditionMessage(e)))

    return(list(
      subject_idx = subject_idx,
      processing_time = processing_time,
      status = "error",
      error = e$message
    ))
  })
}

#' Create all design matrices for HRF parameter grid
#'
#' Generates design matrices for every HRF parameter combination in the grid.
#' Returns a 3D array where the third dimension indexes different HRF models.
#' This is the computationally intensive "Step 2b" of the allHRFs pipeline.
#'
#' @inheritParams EVs_Param
#' @inheritParams nT_Param
#' @inheritParams TR_Param
#' @inheritParams hrf_grid_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams verbose_Param
#' @inheritParams subject_idx_Param
#'
#' @return Named list with elements:
#'   \item{array}{3D array with dimensions (nT × nK × nModels), where nK is the
#'     number of regressors and nModels is the number of HRF parameter combinations}
#'   \item{list}{List of length nModels, where each element is a 2D design matrix
#'     (nT × nK) for one specific HRF parameter combination}
#'
#' @keywords internal
create_all_design_matrices <- function(EVs, nT, TR, hrf_grid, onsets, offsets, verbose, subject_idx) {
  if(verbose > 1) cat("Subject", subject_idx, ": Creating", nrow(hrf_grid), "design matrices...\n")

  # Get number of fields from first design (to determine array dimensions)
  hrf_params_1 <- extract_hrf_params(hrf_grid, 1)

  taper_start_1 <- get_taper_start(
    a1 = hrf_params_1$a1, b1 = hrf_params_1$b1,
    a2 = hrf_params_1$a2, b2 = hrf_params_1$b2,
    c = hrf_params_1$c, TR = TR, deriv = 0 )

  design_1 <- make_design(
    EVs = EVs, nTime = nT, TR = TR, dHRF = 0,  # No derivatives in allHRFs
    onset = onsets, offset = offsets,
    taper_start = taper_start_1,
    a1 = hrf_params_1$a1, b1 = hrf_params_1$b1, c = hrf_params_1$c,
    a2 = hrf_params_1$a2, b2 = hrf_params_1$b2
  )
  nK <- ncol(design_1$design)

  # Initialize 3D array: nT x nK x nModels
  design_3D <- array(NA, dim = c(nT, nK, nrow(hrf_grid)))
  design_list <- list()
  design_list[[1]] <- design_1$design 
  # design_3D <- array(NA, dim = c(nT, nK, 3))

  # Store first design
  design_3D[,,1] <- design_1$design

  # Create remaining designs
  for (pp in 2:nrow(hrf_grid)) {
    tictoc::tic()
    hrf_params <- extract_hrf_params(hrf_grid, pp)

    taper_start_pp <- get_taper_start(
      a1 = hrf_params$a1, b1 = hrf_params$b1,
      a2 = hrf_params$a2, b2 = hrf_params$b2,
      c = hrf_params$c, TR = TR, deriv = 0 )

    design_pp <- make_design(
      EVs = EVs, nTime = nT, TR = TR, dHRF = 0,
      onset = onsets, offset = offsets,
      taper_start = taper_start_pp,
      a1 = hrf_params$a1, b1 = hrf_params$b1, c = hrf_params$c,
      a2 = hrf_params$a2, b2 = hrf_params$b2
    )
    design_list[[pp]] <- design_pp$design

    design_3D[,,pp] <- design_pp$design
    cat(sprintf("[subj %d] make_design idx=%d: ", subject_idx, pp)); tictoc::toc()
  }

  if(verbose > 1) cat("Subject", subject_idx, ": Design array dims =", paste(dim(design_3D), collapse="x"), "\n")
  return(list(
    array = design_3D,      # The BIG 3D array
    list = design_list      # List of individual 2D matrices
  ))
}

#' Fit multiGLM with all HRF designs
#'
#' Fits multiGLM comparing all HRF parameter combinations simultaneously.
#' This is "Step 2c" of the allHRFs pipeline and typically the most 
#' computationally expensive step.
#'
#' @inheritParams BOLD_xii_Param
#' @param design_3D 3D array of design matrices (timepoints × regressors × models).
#' @inheritParams nuisance_file_Param
#' @inheritParams TR_Param
#' @inheritParams brainstructures_Param
#' @inheritParams hpf_Param
#' @inheritParams scrub_Param
#' @inheritParams resamp_res_Param
#' @inheritParams verbose_Param
#' @inheritParams subject_idx_Param
#'
#' @return Complete multiGLM result object with bestmodel_xii indicating
#'   which HRF model fits best at each brain location.
#'
#' @keywords internal
fit_multiGLM_all_designs <- function(BOLD_xii, design_3D, nuisance_file, TR, brainstructures,
                                     hpf, scrub, resamp_res, verbose, subject_idx, EVs, onsets, offsets) {
  if(verbose > 1) cat("Subject", subject_idx, ": Fitting multiGLM with", dim(design_3D)[3], "designs...\n")

  # Load nuisance regressors
  nuisance <- load_nuisance_regressors(nuisance_file)

  # Build canonical design (dHRF=0, default HRF) for Fstat computation
  design_canonical <- make_design(
    EVs = EVs, nTime = nrow(design_3D), TR = TR, dHRF = 0,
    onset = onsets, offset = offsets
  )$design

  # Fit multiGLM with all designs
  glm_result <- multiGLM(
    BOLD = BOLD_xii,
    brainstructures = brainstructures,
    resamp_res = NULL, # load_bold_data resamples, if asked.
    design = design_3D,  # 3D array with all designs
    design_canonical = design_canonical,
    nuisance = nuisance, # Can be NULL
    scrub = scrub,
    TR = TR,
    hpf = hpf
  )

  return(glm_result)
}

#' Extract HRF parameters from grid by index
#'
#' Extracts HRF parameter values from grid for a specific row index.
#' Computes a2, b2 if not present in the grid.
#'
#' @inheritParams hrf_grid_Param
#' @param idx Integer. Row index in the HRF grid.
#'
#' @return Named list with elements: a1, b1, c, a2, b2.
#'
#' @keywords internal
extract_hrf_params <- function(hrf_grid, idx) {
  # Support both row-indexing and field-indexing styles
  if (length(idx) != 1) stop("extract_hrf_params only supports a single index.")

  params <- as.list(hrf_grid[idx, ])

  # Calculate a2, b2 if not provided
  if (is.null(params$a2) || !("a2" %in% names(params))) {
    params$a2 <- (16 / sqrt(6)) * sqrt(params$a1) * sqrt(params$b1)
  }
  if (is.null(params$b2) || !("b2" %in% names(params))) {
    params$b2 <- params$b1
  }

  return(params)
}

#' Extract best parameters across all subjects
#'
#' Collects best-fitting HRF parameters from all subjects and combines
#' into a single data frame. This implements "Step 2d" of the allHRFs pipeline.
#'
#' @inheritParams subject_results_Param
#' @inheritParams hrf_grid_Param
#' @inheritParams verbose_Param
#'
#' @return Data frame with columns: a1, b1, c, voxel, subject.
#'
#' @keywords internal
extract_best_params_all_subjects <- function(subject_results, hrf_grid, verbose) {
  if(verbose > 0) cat("Extracting best parameters from", length(subject_results), "subjects\n")

  best_params_df <- NULL

  # For each subject:
  for(i in 1:length(subject_results)) {

    # Skip failed subjects
    if(subject_results[[i]]$status != "success") next
    #&# Store failed subjects for user to debug

    # 1. Get the best model indices from GLM results
    bestmodel_xii <- subject_results[[i]]$glm_result$bestmodel_xii
    bestmodel_mat <- as.matrix(bestmodel_xii)

    # 2. Convert model indices to parameter values using vectorized approach
    bestmodel_vec <- as.vector(bestmodel_mat)
    best_params_df_i <- data.frame(
      a1 = hrf_grid$a1[bestmodel_vec],
      b1 = hrf_grid$b1[bestmodel_vec],
      c = hrf_grid$c[bestmodel_vec],
      voxel = seq_along(bestmodel_vec),
      subject = i  # Subject index
    )

    # 3. Accumulate across subjects
    best_params_df <- rbind(best_params_df, best_params_df_i)
  }

  return(best_params_df)
}

#' Input validation for fit_allHRFs
#'
#' Validates inputs specific to the allHRFs pipeline, including HRF grid
#' structure and parameter combinations.
#'
#' @inheritParams BOLD_Param
#' @inheritParams EVs_Param
#' @inheritParams nuisance_Param
#' @inheritParams hrf_grid_Param
#' @inheritParams n_cores_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @inheritParams verbose_Param
#'
#' @return TRUE if validation passes (called for side effects).
#'
#' @keywords internal
validate_inputs_allHRFs <- function(BOLD, EVs, nuisance, hrf_grid, n_cores, onsets, offsets, verbose) {
  # Basic input validation (similar to fit_workingHRF)
  if(length(BOLD) == 0) stop("BOLD must contain at least one file path or xifti object")
  if(length(EVs) != length(BOLD)) stop("EVs must have the same length as BOLD")
  if(!is.null(nuisance) && length(nuisance) != length(BOLD)) {
    stop("nuisance must have the same length as BOLD or be NULL")
  }

  # HRF grid validation
  if(!is.data.frame(hrf_grid)) stop("hrf_grid must be a data frame")
  required_cols <- c("a1", "b1", "c", "a2", "b2")
  if(!all(required_cols %in% names(hrf_grid))) {
    stop("hrf_grid must contain columns: ", paste(required_cols, collapse=", "))
  }
  if(nrow(hrf_grid) == 0) stop("hrf_grid must contain at least one row")

  # Cores validation (reuse from fit_workingHRF logic)
  if(!is.numeric(n_cores) || n_cores < 1 || n_cores != round(n_cores)) {
    stop("n_cores must be a positive integer")
  }

  return(TRUE)
}