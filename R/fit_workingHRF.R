
#' Working-HRF fit for a single subject given a preloaded BOLD xifti
#'
#' @param BOLD_xii Loaded \code{xifti} object (output of \code{load_bold_data}).
#' @param nuisance_file Character or NULL. Path to nuisance regressor TSV.
#' @inheritParams scrub_Param
#' @param cfg Named list of pipeline config (built by \code{fit_allHRFs}).
#' @inheritParams subject_idx_Param
#'
#' @return Named list with \code{glm_results}, \code{design_matrix},
#'   \code{field_names}, \code{design_dims}, \code{status = "success"}.
#'
#' @keywords internal
fit_working_one_subject <- function(BOLD_xii, nuisance_file, scrub, cfg, subject_idx) {

  if (cfg$verbose > 1) cat("Subject", subject_idx, ": working-HRF design + GLM\n")

  wh <- cfg$working_hrf
  taper_start <- get_taper_start(a1 = wh$a1, b1 = wh$b1, a2 = wh$a2, b2 = wh$b2,
                                 c = wh$c, TR = cfg$TR, deriv = 0)

  dHRF <- if (isTRUE(cfg$derivatives)) 2 else 0

  design <- make_design(
    EVs = cfg$EVs[[subject_idx]], nTime = ncol(BOLD_xii), TR = cfg$TR, dHRF = dHRF,
    onset = cfg$onsets, offset = cfg$offsets,
    taper_start = taper_start,
    a1 = wh$a1, b1 = wh$b1, c = wh$c, a2 = wh$a2, b2 = wh$b2
  )
  design_array <- convert_design_to_array(design$design)

  glm_results <- fit_glm_model(
    BOLD_xii, design_array, nuisance_file, cfg$TR, cfg$brainstructures,
    cfg$hpf, scrub, cfg$verbose, subject_idx
  )

  list(
    glm_results   = glm_results,
    design_matrix = design$design,
    field_names   = design$field_names,
    design_dims   = dim(design$design),
    status        = "success"
  )
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
    pvals_matrix[, i] <- stats::p.adjust(pvals, method = "bonferroni")
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
