#' Fit Best HRF GLM for a Single Subject
#'
#' Fits a voxel-wise GLM under one or both HRF models per voxel: the
#' population-average HRF (\code{"adapted"}) and/or a subject-personalized
#' HRF chosen here by scoring 25 candidates around the pop_avg. Always
#' produces a canonical-HRF baseline (\code{"working"}) alongside.
#'
#' Candidate scoring for personalized mode:
#' \describe{
#'   \item{lookup}{If \code{allHRF_subject} is supplied, RSS is read from the
#'     subject's cached \code{fit_allHRFs(save_rss = TRUE)} qs file and scored
#'     in vectorized form (fast).}
#'   \item{refit}{Otherwise, the candidate designs are refit against the
#'     subject's BOLD via \code{multiGLM} (slow but works for new subjects
#'     that have no cached RSS).}
#' }
#'
#' @param fit_allHRFs_result An \code{hrfs} object returned by
#'   \code{\link{fit_allHRFs}}. Supplies \code{pop_avg} + \code{hrf_grid} via
#'   \code{fit_allHRFs_result$regularize_allHRFs}, and the per-subject qs
#'   cache via \code{fit_allHRFs_result$fit_allHRFs} (its \code{result_paths}
#'   attribute) for lookup-mode personalized scoring.
#' @param BOLD_file Character. Path to subject's CIFTI file.
#' @param EVs Event data for this subject.
#' @param nuisance_file Character. Path to nuisance regressor file (or NULL).
#' @param TR Numeric. Repetition time in seconds.
#' @param use Character vector. Subset of \code{c("personalized", "adapted")};
#'   pass either (or both) to request the corresponding HRF fit alongside the
#'   always-present working fit. Defaults to both.
#' @param subject_idx Integer or NULL. When non-NULL, personalized mode uses
#'   lookup scoring against \code{fit_allHRFs_result}'s cached RSS for that subject. When
#'   NULL, personalized mode falls back to refit scoring (multiGLM against
#'   this subject's BOLD) -- use NULL for new subjects not in the original
#'   \code{fit_allHRFs} run.
#' @param a1_offsets,b1_offsets Numeric vectors of pop_avg offsets to scan
#'   for personalized candidates. Defaults give 5x5 = 25 candidates.
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
#' @param log_file Character path or NULL. If set, per-mode progress lines are
#'   appended here (with timestamp + label + pid prefix). Useful for monitoring
#'   parallel mode workers.
#' @param n_cores Integer. Number of workers for mode-level parallelism
#'   (working / personalized / adapted run on independent workers via
#'   \code{parallel::parLapplyLB} on a PSOCK cluster, capped at the number
#'   of requested modes). Default \code{1} = pure sequential \code{lapply}.
#'   PSOCK works on macOS / Linux / Windows.
#' @param verbose Integer. Verbosity level.
#'
#' @return A list with class \code{"bestHRF"} containing:
#'   \describe{
#'     \item{working}{Canonical HRF fit: betas + contrasts as xifti (always present).}
#'     \item{adapted}{Population-average HRF fit (when \code{"adapted"} in use).}
#'     \item{personalized}{Personalized HRF fit (when \code{"personalized"} in
#'       use). Includes \code{$winning_candidate_id} and \code{$candidate_scores}.}
#'     \item{df}{Degrees of freedom (shared across modes).}
#'     \item{contrast_matrix}{The contrast matrix A used.}
#'   }
#'
#' @export
fit_bestHRF <- function(fit_allHRFs_result,
                        BOLD_file,
                        EVs,
                        nuisance_file = NULL,
                        TR,
                        use = c("personalized", "adapted"),
                        subject_idx = NULL,
                        a1_offsets = c(-2, -1, 0, 1, 2),
                        b1_offsets = c(-0.5, -0.25, 0, 0.25, 0.5),
                        contrasts = NULL,
                        working_hrf = list(a1 = 6, b1 = 1, c = 1/6),
                        brainstructures = c("left", "right"),
                        resamp_res = 10000,
                        hpf = 0.01,
                        onsets = FALSE,
                        offsets = FALSE,
                        p_adjust_method = "BH",
                        log_file = NULL,
                        n_cores = 1,
                        verbose = 1) {

  use <- match.arg(use, several.ok = TRUE)

  fit_allHRFs_result$fit_workingHRF <- NULL
  fit_allHRFs_result$session_info   <- NULL
  fit_allHRFs_result$fit_allHRFs$subject_results    <- NULL
  fit_allHRFs_result$fit_allHRFs$best_params_results <- NULL

  cfg <- list(
    fit_allHRFs_result = fit_allHRFs_result,
    BOLD_file       = BOLD_file,
    EVs             = EVs,
    nuisance_file   = nuisance_file,
    TR              = TR,
    use             = use,
    subject_idx     = subject_idx,
    a1_offsets      = a1_offsets,
    b1_offsets      = b1_offsets,
    contrasts       = contrasts,
    working_hrf     = working_hrf,
    brainstructures = brainstructures,
    resamp_res      = resamp_res,
    hpf             = hpf,
    onsets          = onsets,
    offsets         = offsets,
    p_adjust_method = p_adjust_method,
    verbose         = verbose
  )

  validate_bestHRF_inputs(cfg)

  pop_avg  <- cfg$fit_allHRFs_result$regularize_allHRFs$pop_avg
  hrf_grid <- cfg$fit_allHRFs_result$regularize_allHRFs$hrf_grid

  # Lookup mode is opt-in via subject_idx. The qs cache holds both per-model
  # RSS and canonical-HRF Fstat in mGLM0s (GLM_multi compares canonical vs null).
  allHRF_subject <- if ("personalized" %in% cfg$use && !is.null(cfg$subject_idx)) {
    list(qs_path = resolve_qs_path(attr(cfg$fit_allHRFs_result$fit_allHRFs, "result_paths")[cfg$subject_idx]))
  } else {
    NULL
  }

  # Load BOLD smoothed-unscaled (same as fit_allHRFs, so cached-RSS lookup is
  # consistent with refit scoring). PSC-scale inline for the GLM fits below
  # -- multiGLM rejects PSC-scaled input, so refit-mode scoring needs the
  # smoothed-unscaled xii.
  if (cfg$verbose > 0) cat("Loading BOLD data...\n")
  bold_data <- load_bold_data(cfg$BOLD_file, cfg$brainstructures, cfg$resamp_res,
                              smoothing = TRUE, scale = FALSE)
  BOLD_xii_unscaled <- bold_data$BOLD_xii

  mat <- as.matrix(BOLD_xii_unscaled)  # n_voxels x nT
  mu  <- rowMeans(mat, na.rm = TRUE)
  if (any(!is.finite(mu) | abs(mu) < 1e-6)) {
    stop("Near-zero or non-finite voxel mean detected; PSC scaling invalid.")
  }
  mat <- ((mat - mu) / mu) * 100
  y   <- t(mat)  # nT x nV
  xii <- BOLD_xii_unscaled  # template only; package_results overwrites data
  nT  <- nrow(y)
  nV  <- ncol(y)

  # Nuisance: raw matrix used by score_candidates_refit (multiGLM adds its own
  # DCT); nuisance_block (raw + DCT) is used by the per-mode GLM fits.
  nuisance_mat <- if (!is.null(cfg$nuisance_file)) {
    as.matrix(utils::read.table(cfg$nuisance_file, header = FALSE))
  } else {
    NULL
  }
  nuisance_block <- build_nuisance_block(nuisance_mat, nT, cfg$TR, cfg$hpf)

  # n_task / contrasts. Use working_hrf params for task derivation since task
  # structure depends on EVs, not HRF shape.
  first_td <- build_task_design(cfg$EVs, nT, cfg$TR,
                                cfg$working_hrf$a1, cfg$working_hrf$b1, cfg$working_hrf$c,
                                cfg$onsets, cfg$offsets)
  n_task <- ncol(first_td$design)
  A <- if (is.null(cfg$contrasts)) diag(n_task) else cfg$contrasts
  if (ncol(A) != n_task) stop("Contrast matrix has ", ncol(A), " columns but there are ", n_task, " task regressors")

  # Resolve HRF maps for each requested mode. Personalized may need BOLD_xii
  # (smoothed-UNSCALED) + raw nuisance (refit) or qs_path (lookup).
  hrf_resolutions <- setNames(
    lapply(cfg$use, resolve_hrf_map,
           cfg = cfg, pop_avg = pop_avg, hrf_grid = hrf_grid,
           allHRF_subject = allHRF_subject,
           BOLD_xii = BOLD_xii_unscaled, nuisance = nuisance_mat),
    cfg$use
  )

  if (cfg$verbose > 0) cat("fit_bestHRF: modes =", paste(cfg$use, collapse = ", "),
                       "(", nrow(hrf_resolutions[[1]]$map), "voxels )\n")

  # Job list: working (always) + each requested mode. The personalized job
  # carries winning_candidate_id + candidate_scores so packaging can stash
  # them on the result without re-touching resolve_hrf_map output.
  working_map <- data.frame(voxel = 1:nV, a1 = cfg$working_hrf$a1,
                            b1 = cfg$working_hrf$b1, c = cfg$working_hrf$c)
  jobs <- list(list(label = "working", map = working_map,
                    winning_candidate_id = NULL, candidate_scores = NULL))
  for (m in cfg$use) {
    res <- hrf_resolutions[[m]]
    jobs[[length(jobs) + 1]] <- list(
      label = m, map = res$map,
      winning_candidate_id = res$winning_candidate_id,
      candidate_scores     = res$candidate_scores
    )
  }

  # Concurrency layer: returns raw numeric matrices. xifti conversion happens
  # sequentially below so we don't push xifti objects through PSOCK exports.
  raw_results <- run_modes_parallel(
    jobs = jobs, n_cores = n_cores, cfg = cfg,
    y = y, nT = nT, nuisance_block = nuisance_block,
    A = A, n_task = n_task, nV = nV, log_file = log_file
  )

  sections <- mapply(package_mode_section, raw_results, jobs,
                     MoreArgs = list(xii = xii), SIMPLIFY = FALSE)
  names(sections) <- vapply(jobs, `[[`, character(1), "label")

  if (cfg$verbose > 0) cat("fit_bestHRF complete.\n")

  result <- c(sections, list(df = raw_results[[1]]$df, contrast_matrix = A))
  class(result) <- "bestHRF"
  result
}


#' Dispatch Mode Fits Sequentially or via a PSOCK Cluster
#'
#' Runs \code{\link{fit_one_mode}} once per job. With \code{n_cores = 1}
#' (default) uses \code{lapply} for behaviorally identical sequential fits.
#' With \code{n_cores > 1} builds a PSOCK cluster via
#' \code{\link{setup_parallel_cluster}} and dispatches with
#' \code{parallel::parLapplyLB}. PSOCK works on macOS / Linux / Windows.
#'
#' @param jobs List of jobs (label, map, winning_candidate_id, candidate_scores).
#' @param n_cores Integer. Workers, capped at \code{length(jobs)}.
#' @param cfg The fit_bestHRF cfg list; supplies EVs, TR, onsets, offsets,
#'   p_adjust_method, verbose.
#' @param y,nT,nuisance_block,A,n_task,nV See \code{\link{fit_one_mode}}.
#' @param log_file Character path or NULL. Also passed as the cluster
#'   \code{outfile} so worker stdout/stderr lands there alongside the
#'   per-mode timestamp lines.
#' @return List of raw fit results, one per job (same order as \code{jobs}).
#' @keywords internal
run_modes_parallel <- function(jobs, n_cores, cfg,
                               y, nT, nuisance_block,
                               A, n_task, nV,
                               log_file = NULL) {
  worker_fn <- function(j) {
    fit_one_mode(j, cfg, y, nT, nuisance_block, A, n_task, nV, log_file)
  }

  if (!isTRUE(n_cores > 1)) return(lapply(jobs, worker_fn))

  n_workers <- min(n_cores, length(jobs))
  if (cfg$verbose > 0) cat("run_modes_parallel: PSOCK cluster with", n_workers, "workers\n")

  cluster_outfile <- if (!is.null(log_file)) log_file else ""
  cl <- parallel::makeCluster(n_workers, outfile = cluster_outfile)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  setup_parallel_cluster(cl, cfg$verbose,
                         c("cfg", "y", "nT", "nuisance_block",
                           "A", "n_task", "nV", "log_file"))

  parallel::parLapplyLB(cl, jobs, worker_fn)
}


#' Fit One Mode (working / adapted / personalized) and Return Raw Matrices
#'
#' Wraps a single call to \code{\link{fit_all_voxels}} for one mode's HRF map.
#' Returns the raw numeric matrices only -- xifti conversion happens
#' sequentially in the caller so parallel workers don't have to push xifti
#' objects back through the fork pipe.
#'
#' @param j List with elements \code{label} (character) and \code{map}
#'   (HRF map data.frame).
#' @param cfg The fit_bestHRF cfg list; supplies EVs, TR, onsets, offsets,
#'   p_adjust_method, verbose.
#' @param y,nT,nuisance_block,A,n_task,nV See \code{\link{fit_all_voxels}}.
#' @param log_file Character path or NULL. If set, appends a timestamp +
#'   label + pid line at start and end of the fit.
#' @return List of raw numeric matrices from \code{fit_all_voxels()}.
#' @keywords internal
fit_one_mode <- function(j, cfg, y, nT, nuisance_block,
                         A, n_task, nV, log_file = NULL) {
  log_msg <- function(msg) {
    if (!is.null(log_file)) {
      cat(sprintf("[%s] [%s] [pid %d] %s\n",
                  format(Sys.time()), j$label, Sys.getpid(), msg),
          file = log_file, append = TRUE)
    }
  }
  log_msg("starting fit")
  raw <- fit_all_voxels(j$map, y, cfg$EVs, nT, cfg$TR, nuisance_block,
                        A, n_task, nV, j$label,
                        cfg$onsets, cfg$offsets, cfg$p_adjust_method, cfg$verbose)
  log_msg(sprintf("done (df = %d)", raw$df))
  raw
}


#' Resolve HRF Map for One Mode
#'
#' Returns the HRF map (and any per-mode diagnostics) for one requested mode.
#' \code{"adapted"} is a pure pop_avg lookup; \code{"personalized"} delegates
#' to \code{\link{resolve_personalized_hrf_map}} which scores candidates.
#'
#' @param mode One of \code{"personalized"} or \code{"adapted"}.
#' @param cfg The fit_bestHRF cfg list.
#' @param pop_avg,hrf_grid Pulled from \code{cfg$fit_allHRFs_result$regularize_allHRFs} by
#'   the caller (passed explicitly so this fn doesn't re-extract per call).
#' @param allHRF_subject,BOLD_xii,nuisance See
#'   \code{\link{resolve_personalized_hrf_map}}; ignored when
#'   \code{mode == "adapted"}.
#' @return List with elements \code{map} (data.frame with columns
#'   \code{voxel, a1, b1, c}), \code{winning_candidate_id} (integer or NULL),
#'   and \code{candidate_scores} (numeric vector or NULL).
#' @keywords internal
resolve_hrf_map <- function(mode, cfg, pop_avg, hrf_grid,
                            allHRF_subject = NULL,
                            BOLD_xii = NULL, nuisance = NULL) {
  if (mode == "adapted") {
    list(
      map = data.frame(voxel = pop_avg$voxel,
                       a1    = pop_avg$a1,
                       b1    = pop_avg$b1,
                       c     = pop_avg$c),
      winning_candidate_id = NULL,
      candidate_scores     = NULL
    )
  } else {
    resolve_personalized_hrf_map(cfg, pop_avg, hrf_grid,
                                 allHRF_subject, BOLD_xii, nuisance)
  }
}


#' Resolve Personalized HRF Map
#'
#' Builds 25 candidate maps around \code{pop_avg}, scores them via the
#' lookup path (when \code{allHRF_subject} supplies cached RSS + Fstat) or
#' the refit path (otherwise), and returns the winning candidate's map plus
#' diagnostics.
#'
#' @param cfg The fit_bestHRF cfg list; supplies EVs, TR, brainstructures,
#'   hpf, onsets, offsets, a1_offsets, b1_offsets, verbose.
#' @param pop_avg,hrf_grid Pulled from \code{cfg$fit_allHRFs_result$regularize_allHRFs}.
#' @param allHRF_subject Optional list with \code{qs_path} pointing to the
#'   subject's fit_allHRFs qs cache (which carries both per-model RSS and
#'   canonical-HRF Fstat). When supplied, lookup scoring is used; NULL
#'   forces the refit path.
#' @param BOLD_xii Smoothed/scaled BOLD as a xifti (from
#'   \code{load_bold_data}). Required for refit scoring.
#' @param nuisance Raw nuisance matrix (multiGLM adds its own DCT). May
#'   be NULL.
#' @return List with \code{map}, \code{winning_candidate_id},
#'   \code{candidate_scores}.
#' @keywords internal
resolve_personalized_hrf_map <- function(cfg, pop_avg, hrf_grid,
                                         allHRF_subject, BOLD_xii, nuisance) {
  cands <- build_candidate_maps(pop_avg, hrf_grid,
                                a1_offsets = cfg$a1_offsets,
                                b1_offsets = cfg$b1_offsets,
                                verbose    = max(0, cfg$verbose - 1))$candidate_maps

  scores <- if (!is.null(allHRF_subject)) {
    if (cfg$verbose > 0) cat("Scoring candidates via RSS lookup (qs: ",
                             basename(allHRF_subject$qs_path), ")\n", sep = "")
    cache <- load_qs_cache(allHRF_subject$qs_path)
    score_candidates_lookup(cands, cache$RSS, cache$Fstat)
  } else {
    if (cfg$verbose > 0) cat("Scoring candidates via multiGLM refit\n")
    score_candidates_refit(
      candidate_maps  = cands, hrf_grid = hrf_grid,
      BOLD_xii        = BOLD_xii, EVs = cfg$EVs, nuisance = nuisance, TR = cfg$TR,
      brainstructures = cfg$brainstructures, hpf = cfg$hpf,
      onsets = cfg$onsets, offsets = cfg$offsets
    )
  }

  winner <- pick_winning_candidate(scores)
  cm <- cands[[winner]]
  list(
    map = data.frame(voxel = cm$voxel, a1 = cm$a1, b1 = cm$b1, c = cm$c),
    winning_candidate_id = winner,
    candidate_scores     = scores
  )
}


#' Resolve a qs cache path from result_paths attr.
#' Absolute -> pass through. Relative -> prepend getOption("hrf.output_root").
#' @keywords internal
resolve_qs_path <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) return(path)
  if (file.exists(path)) return(path)
  root <- getOption("hrf.output_root", default = NULL)
  if (!is.null(root)) {
    if (startsWith(path, "/")) {
      # Stale absolute path -- extract the fit_allHRFs/... tail and re-anchor to root.
      tail <- regmatches(path, regexpr("fit_allHRFs/.*", path))
      if (length(tail) && nzchar(tail)) {
        candidate <- file.path(root, tail)
        if (file.exists(candidate)) return(candidate)
      }
    } else {
      candidate <- file.path(root, path)
      if (file.exists(candidate)) return(candidate)
    }
  }
  stop("qs cache not found at '", path,
       "' and could not resolve via options(hrf.output_root). ",
       "Either set the option to the current pipeline output dir, ",
       "or check the path.")
}


#' Load RSS + Canonical Fstat from a fit_allHRFs qs Cache
#'
#' Reads, for one subject, the stacked (cortexL + cortexR) per-model RSS
#' matrix \emph{and} the canonical-HRF Fstat vector written by
#' \code{fit_allHRFs(save_rss = TRUE)}.
#'
#' @param qs_path Character path to the subject's qs cache file.
#' @return List with elements \code{RSS} (numeric matrix
#'   \code{[n_voxels x n_models]}) and \code{Fstat} (numeric vector,
#'   length \code{n_voxels}).
#' @keywords internal
load_qs_cache <- function(qs_path) {
  qs_obj <- load_object(file_path = qs_path, delete_after_load = FALSE)
  RSS <- rbind(qs_obj$glm_result$mGLM0s$cortexL$RSS,
               qs_obj$glm_result$mGLM0s$cortexR$RSS)
  if (is.null(RSS) || nrow(RSS) == 0L) {
    stop("RSS not found in qs cache '", qs_path,
         "'. Was fit_allHRFs run with save_rss = TRUE?")
  }
  Fstat <- c(qs_obj$glm_result$mGLM0s$cortexL$Fstat,
             qs_obj$glm_result$mGLM0s$cortexR$Fstat)
  if (is.null(Fstat) || length(Fstat) == 0L) {
    stop("Fstat not found in qs cache '", qs_path, "' (mGLM0s$cortex{L,R}$Fstat).")
  }
  list(RSS = RSS, Fstat = Fstat)
}


#' Validate fit_bestHRF Inputs
#'
#' Checks: \code{cfg$fit_allHRFs_result} has \code{regularize_allHRFs} + \code{fit_allHRFs}
#' sub-objects with the expected columns, and (for personalized lookup) the
#' subject's qs cache path exists.
#'
#' @param cfg The fit_bestHRF cfg list.
#' @keywords internal
validate_bestHRF_inputs <- function(cfg) {
  reg <- cfg$fit_allHRFs_result$regularize_allHRFs
  if (is.null(reg$pop_avg)) {
    stop("fit_allHRFs_result$regularize_allHRFs$pop_avg is NULL. Pass an hrfs object from fit_allHRFs().")
  }
  if (is.null(reg$hrf_grid)) {
    stop("fit_allHRFs_result$regularize_allHRFs$hrf_grid is NULL. Pass an hrfs object from fit_allHRFs().")
  }
  required_pa <- c("voxel", "a1", "b1", "c")
  missing_pa <- setdiff(required_pa, names(reg$pop_avg))
  if (length(missing_pa) > 0L) {
    stop("fit_allHRFs_result$regularize_allHRFs$pop_avg is missing column(s): ",
         paste(missing_pa, collapse = ", "))
  }
  required_grid <- c("a1", "b1", "c", "time_to_peak", "FWHM")
  missing_grid <- setdiff(required_grid, names(reg$hrf_grid))
  if (length(missing_grid) > 0L) {
    stop("fit_allHRFs_result$regularize_allHRFs$hrf_grid is missing column(s): ",
         paste(missing_grid, collapse = ", "))
  }

  if ("personalized" %in% cfg$use && !is.null(cfg$subject_idx)) {
    paths <- attr(cfg$fit_allHRFs_result$fit_allHRFs, "result_paths")
    if (is.null(paths)) {
      stop("fit_allHRFs_result$fit_allHRFs has no 'result_paths' attribute -- ",
           "cannot run personalized lookup. Pass subject_idx = NULL to use refit mode.")
    }
    qs_path <- resolve_qs_path(paths[cfg$subject_idx])
    if (is.null(qs_path) || is.na(qs_path) || !nzchar(qs_path)) {
      stop("attr(fit_allHRFs_result$fit_allHRFs, 'result_paths')[", cfg$subject_idx,
           "] is empty -- did fit_allHRFs run with save_rss = TRUE and save subject ",
           cfg$subject_idx, "?")
    }
    if (!file.exists(qs_path)) {
      stop("qs cache for subject ", cfg$subject_idx, " not found at: ", qs_path)
    }
  }
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

  beta_out  <- matrix(NA_real_, nrow = n_vox, ncol = n_task)
  est_out   <- matrix(NA_real_, nrow = n_vox, ncol = n_contrasts)
  SE_out    <- matrix(NA_real_, nrow = n_vox, ncol = n_contrasts)
  tstat_out <- matrix(NA_real_, nrow = n_vox, ncol = n_contrasts)
  pval_out  <- matrix(NA_real_, nrow = n_vox, ncol = n_contrasts)

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

    beta_out[j, ]  <- beta_task
    est_out[j, ]   <- contr$est
    SE_out[j, ]    <- contr$SE
    tstat_out[j, ] <- contr$tstat
    pval_out[j, ]  <- contr$pval
  }

  list(beta = beta_out, est = est_out, SE = SE_out,
       tstat = tstat_out, pval = pval_out, df = df)
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


#' Package One Mode's Raw Fit Into a Section + Metadata
#'
#' Wraps \code{\link{package_results}} and attaches the per-mode metadata
#' (HRF assignment map; for personalized, also the winning candidate id and
#' the full vector of candidate scores).
#'
#' @param raw Raw fit list returned by \code{\link{fit_one_mode}}.
#' @param j Corresponding job entry built in \code{\link{fit_bestHRF}};
#'   carries \code{map}, and (for personalized) \code{winning_candidate_id}
#'   and \code{candidate_scores}.
#' @param xii Template xifti for \code{\link{package_results}}.
#' @return List in the shape of a single \code{fit_bestHRF} result section.
#' @keywords internal
package_mode_section <- function(raw, j, xii) {
  sec <- package_results(raw, xii)
  sec$hrf_assignments <- j$map
  if (!is.null(j$winning_candidate_id)) sec$winning_candidate_id <- j$winning_candidate_id
  if (!is.null(j$candidate_scores))      sec$candidate_scores      <- j$candidate_scores
  sec
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


#' Build Candidate HRF Maps for One Subject's Personalization
#'
#' Generates the candidate HRF maps that fit_bestHRF scores when picking a
#' subject's personalized HRF. For each combination of \code{a1_offset} and
#' \code{b1_offset}, shifts the population-average snapped a1/b1 by the
#' offset (clamped to the hrf_grid valid range) and snaps back to the nearest
#' grid point. \code{c} is taken verbatim from \code{pop_avg}.
#'
#' Port of the per-subject body of \code{regularize_allHRFs::create_candidate_maps}.
#' During the migration window this lives in both files; regularize keeps using
#' its copy until Phase 2 strips the candidate-fitting machinery from
#' \code{regularize_allHRFs.R}.
#'
#' @param pop_avg Data frame with columns \code{voxel}, \code{a1},
#'   \code{b1}, \code{c} (the regularize pop_avg).
#' @param hrf_grid HRF grid with at least \code{a1}, \code{b1}, \code{c},
#'   \code{time_to_peak}, \code{FWHM} columns.
#' @param a1_offsets Numeric vector of a1 offsets to scan (default
#'   \code{c(-2, -1, 0, 1, 2)}).
#' @param b1_offsets Numeric vector of b1 offsets to scan (default
#'   \code{c(-0.5, -0.25, 0, 0.25, 0.5)}).
#' @param verbose Integer verbosity level.
#'
#' @return List with two elements:
#'   \describe{
#'     \item{candidate_maps}{List of per-candidate data.frames with columns
#'       voxel, a1, b1, c, model_idx, time_to_peak, FWHM, offset_id,
#'       a1_offset, b1_offset.}
#'     \item{offset_combos}{Data.frame of the offset combinations scanned.}
#'   }
#' @keywords internal
build_candidate_maps <- function(pop_avg, hrf_grid,
                                 a1_offsets = c(-2, -1, 0, 1, 2),
                                 b1_offsets = c(-0.5, -0.25, 0, 0.25, 0.5),
                                 verbose = 1) {
  offset_combos <- expand.grid(a1_offset = a1_offsets, b1_offset = b1_offsets)
  n_candidates <- nrow(offset_combos)
  if (verbose > 0) cat("Creating", n_candidates, "candidate maps\n")

  a1_min <- min(hrf_grid$a1); a1_max <- max(hrf_grid$a1)
  b1_min <- min(hrf_grid$b1); b1_max <- max(hrf_grid$b1)
  hrf_grid$key <- paste(hrf_grid$a1, hrf_grid$b1, hrf_grid$c, sep = "_")

  candidate_maps <- vector("list", n_candidates)
  for (i in seq_len(n_candidates)) {
    a1_off <- offset_combos$a1_offset[i]
    b1_off <- offset_combos$b1_offset[i]

    a1_shifted <- pmin(pmax(pop_avg$a1 + a1_off, a1_min), a1_max)
    b1_shifted <- pmin(pmax(pop_avg$b1 + b1_off, b1_min), b1_max)
    c_vals     <- pop_avg$c

    snapped <- snap_to_grid(a1_shifted, b1_shifted, c_vals, hrf_grid)

    keys     <- paste(snapped$a1, snapped$b1, snapped$c, sep = "_")
    grid_idx <- match(keys, hrf_grid$key)

    candidate_maps[[i]] <- data.frame(
      voxel        = pop_avg$voxel,
      a1           = snapped$a1,
      b1           = snapped$b1,
      c            = snapped$c,
      model_idx    = grid_idx,
      time_to_peak = hrf_grid$time_to_peak[grid_idx],
      FWHM         = hrf_grid$FWHM[grid_idx],
      offset_id    = i,
      a1_offset    = a1_off,
      b1_offset    = b1_off
    )
  }

  if (verbose > 0) {
    all_keys <- unique(unlist(lapply(candidate_maps, function(cm) {
      paste(cm$a1, cm$b1, cm$c, sep = "_")
    })))
    cat("Total unique HRF combos across candidates:", length(all_keys), "\n")
  }

  list(candidate_maps = candidate_maps, offset_combos = offset_combos)
}


#' Pick Winning Candidate
#'
#' Selects the candidate with the lowest weighted RSS. Trivial \code{which.min}
#' wrapper kept as a named function for clarity at call sites.
#'
#' @param scores Numeric vector of weighted RSS scores (one per candidate).
#' @return Integer index of the candidate with the smallest score.
#' @keywords internal
pick_winning_candidate <- function(scores) {
  which.min(scores)
}


#' Score Candidate Maps for One Subject via RSS Lookup
#'
#' Computes Fstat-weighted RSS for each candidate map using pre-computed
#' per-voxel-per-model RSS. Fastest scoring path; usable when
#' \code{fit_allHRFs} was run with \code{save_rss = TRUE} so RSS matrices
#' are available (loaded either from the .qs cache or already in memory).
#'
#' Port of the per-subject body of
#' \code{regularize_allHRFs::fit_candidate_maps_lookup}, single-subject and
#' source-agnostic (caller decides how to source RSS + Fstat). During the
#' migration window the regularize copy is kept and still drives the
#' aggregate lookup path; the fit_bestHRF copy will drive Phase 3
#' personalization for new subjects.
#'
#' @param candidate_maps List of per-candidate data.frames from
#'   \code{\link{build_candidate_maps}}. Each must have a \code{voxel} and
#'   \code{model_idx} column.
#' @param RSS Numeric matrix \code{[n_voxels x n_models]}. Per-voxel
#'   residual sum of squares from \code{multiGLM} (cortexL + cortexR
#'   stacked).
#' @param fstat Numeric vector of canonical-HRF Fstat values, indexed by
#'   voxel (length = \code{nrow(RSS)}). Used to compute the Fstat-squared
#'   weights over the population-mask voxels.
#'
#' @return Numeric vector of length \code{length(candidate_maps)} giving
#'   Fstat-weighted RSS per candidate (lower is better).
#' @keywords internal
score_candidates_lookup <- function(candidate_maps, RSS, fstat) {
  pop_voxels <- candidate_maps[[1]]$voxel
  fstat_sq <- fstat[pop_voxels]^2
  fstat_sq_sum <- sum(fstat_sq, na.rm = TRUE)

  scores <- numeric(length(candidate_maps))
  for (j in seq_along(candidate_maps)) {
    voxel_idx <- candidate_maps[[j]]$voxel
    model_idx <- candidate_maps[[j]]$model_idx
    voxel_rss <- RSS[cbind(voxel_idx, model_idx)]
    scores[j] <- sum(voxel_rss * fstat_sq, na.rm = TRUE) / fstat_sq_sum
  }
  scores
}


#' Score Candidate Maps for One Subject by Refitting multiGLM
#'
#' Computes Fstat-weighted RSS for each candidate map by refitting all unique
#' HRF designs against the subject's BOLD via \code{multiGLM}. Slower than
#' \code{\link{score_candidates_lookup}} but doesn't require cached RSS, so
#' it's the path used for new subjects without a \code{save_rss = TRUE}
#' \code{fit_allHRFs} run.
#'
#' Port of the per-subject body of
#' \code{regularize_allHRFs::process_candidate_subject_refit}. Same arg
#' shapes as multiGLM expects (BOLD as xifti, raw nuisance matrix, TR/hpf),
#' so the math is bit-identical to regularize's refit path; future cleanup
#' could swap multiGLM for a y/nuisance_block GLM but parity comes first.
#'
#' Fstat for the weights comes from the SAME multiGLM call (canonical-design
#' mGLM0s output), not from a separately-run workingHRF result -- mirrors how
#' \code{process_candidate_subject_refit} sources it.
#'
#' @param candidate_maps List of per-candidate data.frames from
#'   \code{\link{build_candidate_maps}}.
#' @param hrf_grid HRF grid with at least \code{a1, b1, c, a2, b2} columns.
#' @param BOLD_xii Subject's BOLD as a xifti object (already loaded by caller).
#' @param EVs Event data for this subject.
#' @param nuisance Numeric matrix of nuisance regressors (DCT bases will be
#'   added internally by multiGLM via \code{hpf}).
#' @param TR Numeric. Repetition time in seconds.
#' @param brainstructures Character vector passed through to multiGLM.
#' @param hpf Numeric. High-pass filter cutoff (sets DCT basis count).
#' @param onsets,offsets Logical. Include onset/offset regressors.
#'
#' @return Numeric vector of length \code{length(candidate_maps)} giving
#'   Fstat-weighted RSS per candidate (lower is better).
#' @keywords internal
score_candidates_refit <- function(candidate_maps, hrf_grid,
                                   BOLD_xii, EVs, nuisance, TR,
                                   brainstructures = c("left", "right"),
                                   hpf = 0.01,
                                   onsets = FALSE, offsets = FALSE) {
  all_model_idx <- sort(unique(unlist(lapply(candidate_maps, function(cm) unique(cm$model_idx)))))
  unique_grid   <- hrf_grid[all_model_idx, ]
  n_models      <- nrow(unique_grid)
  idx_map       <- setNames(seq_len(n_models), all_model_idx)
  pop_voxels    <- candidate_maps[[1]]$voxel

  nT <- ncol(BOLD_xii)

  design_canonical <- make_design(
    EVs = EVs, nTime = nT, TR = TR, dHRF = 0,
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
      EVs = EVs, nTime = nT, TR = TR, dHRF = 0,
      onset = onsets, offset = offsets,
      taper_start = taper_start_j,
      a1 = unique_grid$a1[j], b1 = unique_grid$b1[j],
      c = unique_grid$c[j],
      a2 = unique_grid$a2[j], b2 = unique_grid$b2[j]
    )$design
  }

  glm_result <- multiGLM(
    BOLD = BOLD_xii, brainstructures = brainstructures,
    resamp_res = NULL, design = design_3D,
    design_canonical = design_canonical,
    nuisance = nuisance, scrub = NULL, TR = TR, hpf = hpf
  )

  RSS   <- rbind(glm_result$mGLM0s$cortexL$RSS,   glm_result$mGLM0s$cortexR$RSS)
  Fstat <- c(   glm_result$mGLM0s$cortexL$Fstat, glm_result$mGLM0s$cortexR$Fstat)

  fstat_sq     <- Fstat[pop_voxels]^2
  fstat_sq_sum <- sum(fstat_sq, na.rm = TRUE)

  scores <- numeric(length(candidate_maps))
  for (j in seq_along(candidate_maps)) {
    compact_idx <- idx_map[as.character(candidate_maps[[j]]$model_idx)]
    voxel_rss   <- RSS[cbind(pop_voxels, compact_idx)]
    scores[j]   <- sum(voxel_rss * fstat_sq, na.rm = TRUE) / fstat_sq_sum
  }
  scores
}
