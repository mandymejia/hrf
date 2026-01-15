
library(dplyr)
library(ciftiTools)
library(here)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

safe_path <- function(...) {
  path <- file.path(...)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  return(path)
}

root_dir <- "/Volumes/LaCie/root/hrf_adaptation"

sess_names <- c('tfMRI_MOTOR','tfMRI_WM','tfMRI_GAMBLING')
run_names <- c('LR','RL')


################################################################################
# WEIGHTS
################################################################################

plot_dir <- "~/Documents/Github/hrf-z/dev/exploration/plots/weights"

# Common zlims for each parameter (1st-99th percentile)
zlims <- list(
  a1 = c(0.29, 0.76),
  b1 = c(1.54, 2.25),
  c = c(11.48, 15.80)
)

# Define runs: (save_dir, label, suffix)
runs <- list(
  # list(save_dir = "validation_bonferroni_0.05/regularize_allHRFs", label = "Non-Smoothed", suffix = ""),
  # list(save_dir = "validation_bonferroni_0.05_smoothed/regularize_allHRFs", label = "Smoothed", suffix = "_smoothed"),
  list(save_dir = "validation/regularize_allHRFs", label = "Non-Smoothed (No MCC)", suffix = "_noMCC")
)

for (run in runs) {
  save_dir <- run$save_dir
  smooth_label <- run$label
  smooth_suffix <- run$suffix

  for (t in 1:3) {
    for (r in 1:2) {
      sess <- paste(sess_names[t], run_names[r], sep='_')
      cat("Processing", smooth_label, "-", sess, "...\n"); flush.console()
      session_data <- readRDS(file.path('~/Documents/Github/HRF-Adaptation-paper/Code/fixtures', paste0('session_data_', sess, '.rds')))
      regularize_allHRFs_result <- readRDS(file.path(root_dir, save_dir, paste0("regularize_allHRFs_result_", sess, "_", length(session_data$BOLD_files), "s.rds")))

      for (param in c("a1", "b1", "c")) {
        results_wls <- regularize_allHRFs_result[["regularized_params"]][[param]][["results_WLS"]]

        # Get unique voxel-weight pairs
        voxel_weights <- results_wls %>%
          distinct(voxel, weight_B)

        # Create full vector
        mask_prop_NA <- regularize_allHRFs_result[["mask_prop_NA"]]
        full_vector <- rep(NA, length(mask_prop_NA))
        full_vector[voxel_weights$voxel] <- voxel_weights$weight_B

        # Plot
        xii_template <- attr(regularize_allHRFs_result, "xii")
        weight_xii <- ciftiTools::newdata_xifti(xii_template, full_vector)

        plot_file <- file.path(plot_dir, paste0(sess, "_", param, smooth_suffix, ".png"))
        plot(
          weight_xii,
          title = paste0("Weights (Intercept-Only Model) - ", sess, " - ", param, " - ", smooth_label),
          color_mode = "sequential",
          NA_color = "#505560",
          zlim = zlims[[param]],
          fname = plot_file
        )
        cat("  Saved:", plot_file, "\n"); flush.console()
      }
    }
  }
}


################################################################################
# AFTER & BEFORE REGULARIZATION
################################################################################

plot_rounded_params <- function(rounded_params, xii_template, subject, sess, param, plot_dir, smooth_suffix) {
  # Parameter-specific zlims based on hrf_grid ranges
  zlims <- list(a1 = c(3, 12), b1 = c(0.5, 2.0), c = c(0, 0.17))

  # Filter to just this subject and parameter
  subject_params <- rounded_params %>%
    dplyr::filter(subject == !!subject) %>%
    dplyr::select(voxel, !!sym(param)) %>%
    dplyr::distinct()

  # Create xifti for plotting
  full_vector <- rep(NA, nrow(as.matrix(xii_template)))
  full_vector[subject_params$voxel] <- subject_params[[param]]
  param_xii <- ciftiTools::newdata_xifti(xii_template, full_vector)

  # Plot - organized by session subdirectory
  plot_dir_sess <- safe_path(plot_dir, "rounded_params", sess)
  plot_file <- file.path(plot_dir_sess, paste0(sess, "_", param, "_subject_", subject, smooth_suffix, ".png"))
  plot(param_xii,
       zlim = zlims[[param]],
       color_mode = "sequential",
       title = paste0("Parameter ", param, " (Rounded) - ", sess, " - Subject ", subject),
       fname = plot_file)
}

plot_raw_params <- function(raw_params, xii_template, subject, sess, param, plot_dir, smooth_suffix) {
  # Parameter-specific zlims based on hrf_grid ranges
  zlims <- list(a1 = c(3, 12), b1 = c(0.5, 2.0), c = c(0, 0.17))

  # Filter to just this subject and parameter
  subject_params_raw <- raw_params %>%
    dplyr::filter(subject == !!subject) %>%
    dplyr::select(voxel, !!sym(param)) %>%
    dplyr::distinct()

  # Create xifti for plotting
  full_vector_raw <- rep(NA, nrow(as.matrix(xii_template)))
  full_vector_raw[subject_params_raw$voxel] <- subject_params_raw[[param]]
  param_xii_raw <- ciftiTools::newdata_xifti(xii_template, full_vector_raw)

  # Plot - organized by session subdirectory
  plot_dir_sess <- safe_path(plot_dir, "raw_params", sess)
  plot_file <- file.path(plot_dir_sess, paste0(sess, "_", param, "_subject_", subject, smooth_suffix, ".png"))
  plot(param_xii_raw,
       zlim = zlims[[param]],
       color_mode = "sequential",
       title = paste0("Parameter ", param, " (Raw) - ", sess, " - Subject ", subject),
       fname = plot_file)
}

plot_raw_after_params <- function(regularize_result, xii_template, subject, sess, param, plot_dir, smooth_suffix) {
  # Parameter-specific zlims based on hrf_grid ranges
  zlims <- list(a1 = c(3, 12), b1 = c(0.5, 2.0), c = c(0, 0.17))

  # Extract continuous after-regularization values from WLS results
  wls_results <- regularize_result$regularized_params[[param]]$results_WLS

  subject_params_after <- wls_results %>%
    dplyr::filter(subject == !!subject) %>%
    dplyr::select(voxel, param_value = param_pred_B) %>%
    dplyr::distinct()

  # Create xifti for plotting
  full_vector_after <- rep(NA, nrow(as.matrix(xii_template)))
  full_vector_after[subject_params_after$voxel] <- subject_params_after$param_value
  param_xii_after <- ciftiTools::newdata_xifti(xii_template, full_vector_after)

  # Plot - organized by session subdirectory
  plot_dir_sess <- safe_path(plot_dir, "raw_after_params", sess)
  plot_file <- file.path(plot_dir_sess, paste0(sess, "_", param, "_subject_", subject, smooth_suffix, ".png"))

  plot(param_xii_after,
       zlim = zlims[[param]],
       color_mode = "sequential",
       title = paste0("Parameter ", param, " (After Reg., Continuous) - ", sess, " - Subject ", subject),
       fname = plot_file)
}

plot_dir <- "~/Documents/Github/hrf-z/dev/exploration/plots"

for (smoothed in c(FALSE, TRUE)) {
  save_dir <- if (smoothed) {
    "validation_bonferroni_0.05_smoothed/regularize_allHRFs"
  } else {
    "validation_bonferroni_0.05/regularize_allHRFs"
  }

  smooth_label <- if (smoothed) "Smoothed" else "Non-Smoothed"
  smooth_suffix <- if (smoothed) "_smoothed" else ""

  for(t in 1:3) {
    for(r in 1:2) {
      sess <- paste(sess_names[t], run_names[r], sep='_')
      session_data <- readRDS(file.path('~/Documents/Github/HRF-Adaptation-paper/Code/fixtures',
                                        paste0('session_data_', sess, '.rds')))
      subject_num <- length(session_data$BOLD_files)
      save_path <- file.path(root_dir, save_dir, paste0("regularize_allHRFs_result_", sess, "_", subject_num, "s.rds"))

      cat("Loading", smooth_label, "-", sess, "....\n"); flush.console()
      regularize_allHRFs_result <- readRDS(save_path)

      # Get both rounded and raw parameters
      rounded_params <- regularize_allHRFs_result$rounded_params$rounded_intercept_only
      raw_params <- regularize_allHRFs_result$best_params_df
      xii_template <- attr(regularize_allHRFs_result, "xii")

      for (subject in c("1", "2", "3", "4")) {
        for (param in c("a1", "b1", "c")) {
          cat("  Plotting subject", subject, "param", param, "...\n"); flush.console()
          # plot_rounded_params(rounded_params, xii_template, subject, sess, param, plot_dir, smooth_suffix)
          # plot_raw_params(raw_params, xii_template, subject, sess, param, plot_dir, smooth_suffix)
          plot_raw_after_params(regularize_allHRFs_result, xii_template, subject, sess, param, plot_dir, smooth_suffix)
        }
      }

      rm(regularize_allHRFs_result); gc()
    }
  }
}

cat("Done!\n")
