library(ciftiTools)
library(dplyr)
library(here)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

save_dir <- "/Volumes/LaCie/root/hrf_adaptation/validation/regularize_allHRFs"
plot_dir <- "~/Documents/Github/hrf-z/dev/exploration/plots"

sess_names <- c('tfMRI_MOTOR','tfMRI_WM','tfMRI_GAMBLING')
run_names <- c('LR','RL')


for(t in 1:3) {
  for(r in 1:2) {
    sess <- paste(sess_names[t], run_names[r], sep='_')
    session_data <- readRDS(file.path('~/Documents/Github/HRF-Adaptation-paper/Code/fixtures', paste0('session_data_', sess, '.rds')))
    subject_num <- length(session_data$BOLD_files)
    save_path <- file.path(save_dir, paste0("regularize_allHRFs_result_", sess, "_", subject_num, "s.rds"))

    cat("Loading....\n"); flush.console()
    regularize_allHRFs_result <- readRDS(save_path)
    cat("Loaded.\n"); flush.console()

    for (subject in c("1", "2", "3", "4")) {
      # Get rounded parameters
      rounded_params <- regularize_allHRFs_result$rounded_params$rounded_intercept_only

      # Filter to just this subject
      subject_params <- rounded_params %>%
        dplyr::filter(subject == !!subject) %>%
        dplyr::select(voxel, c) %>%
        dplyr::distinct()



      # Create xifti for plotting
      xii_template <- attr(regularize_allHRFs_result, "xii")
      full_vector <- rep(NA, nrow(as.matrix(xii_template)))
      full_vector[subject_params$voxel] <- subject_params$c

      c_xii <- ciftiTools::newdata_xifti(xii_template, full_vector)

      # Plot
      plot(c_xii,
           zlim = c(0, 0.17),
           color_mode = "sequential",
           title = paste0("Parameter c - ", sess, " - Subject ", subject, " (Intercept-Only Model)"),
           fname = file.path(plot_dir, "rounded_params", paste0(sess, "_subject_", subject, "_rounded_params.png")))
    }

  }
  rm(regularize_allHRFs_result); gc();
}




