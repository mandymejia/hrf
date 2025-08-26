#' Setup parallel cluster for subject processing
#'
#' Configures cluster workers with required packages, functions, and variables.
#' Synchronizes ciftiTools workbench path across workers and exports necessary
#' functions from calling environment.
#'
#' @param cl Cluster object from makeCluster().
#' @inheritParams verbose_Param
#' @param vars_to_export Character vector of variable names to export to workers.
#'
#'
#' @return NULL (called for side effects)
#'
#' @keywords internal
setup_parallel_cluster <- function(cl, verbose = 1, vars_to_export) {
  if (verbose > 0) cat("Exporting functions and packages to cluster workers\n")

  # Load necessary packages
  parallel::clusterEvalQ(cl, {
    library(hrf)
    library(ciftiTools)
    library(qs)
    library(pryr)# debugging remember to remove from DESCRIPTION
    library(ps) # debugging remember to remove from DESCRIPTION
    library(tictoc) #debugging remember to remove from DESCRIPTION
  })

  # Sync workbench path
  current_wb_path <- ciftiTools.getOption('wb_path')
  if (!is.null(current_wb_path)) {
    parallel::clusterCall(cl, function(path) {
      ciftiTools.setOption('wb_path', path)
    }, current_wb_path)

    if (verbose > 1) cat("Set wb_path on workers to:", current_wb_path, "\n")
  } else {
    if (verbose > 0) stop("No wb_path found in main session\n")
  }

  # Export necessary functions that arent apart of the NAMESPACE
  parallel::clusterExport(cl, c(
    "process_single_subject_df",
    "create_design_matrix",
    "fit_glm_model",
    "load_nuisance_regressors",
    "convert_design_to_array",
    "load_bold_data"
    # "get_onset_offset_tasks",
    # fit_allHRFs
    # "process_entire_subject",
    # "create_all_design_matrices",
    # "fit_multiGLM_all_designs",
    # "extract_hrf_params",
    # file handling
    # "save_object",
    # "load_object",
    # "try_save_debug"
  ), envir = parent.frame())  # ensures calling environment is used

  # Export vars if given
  if (!is.null(vars_to_export)) {
    existing_vars <- vars_to_export[sapply(vars_to_export, exists, envir = parent.frame())]
    if (length(existing_vars) > 0) {
      parallel::clusterExport(cl, varlist = existing_vars, envir = parent.frame())
      if (verbose > 1) cat("Exported variables to cluster:", paste(existing_vars, collapse=", "), "\n")
    }
  }

  parallel::clusterEvalQ(cl, Sys.sleep(runif(1, 0.1, 0.5))) # Stagger workers for stability
}
