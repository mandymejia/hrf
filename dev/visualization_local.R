library(ciftiTools)
library(tictoc)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)


safe_path <- function(...) {
  path <- file.path(...)
  dir <- dirname(path)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  return(path)
}

devtools::load_all("~/Documents/Github/hrf-z") # hrf-z is Zeshawn's branch with latest hrf package changes.
ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

sess_names <- c('tfMRI_MOTOR','tfMRI_WM','tfMRI_GAMBLING')
run_names <- c('LR','RL')


for(t in 1:3){
  for(r in 1){
    sess <- paste(sess_names[t], run_names[r], sep='_')
    session_data <- readRDS(file.path('~/Documents/Github/HRF-Adaptation-paper/Code/validation/fixtures', paste0('session_data_', sess, '.rds')))

    subj_counts <- c(length(session_data$BOLD_files), 500, 200, 100, 25)

    for (n_subj in subj_counts) {
      n_use <- min(n_subj, length(session_data$BOLD_files))
      # save_dir <- "/N/project/hrf_adaptation/validation/regularize_allHRFs"
      save_dir <- "/Volumes/LaCie/root/hrf_adaptation/validation/regularize_allHRFs"
      save_path <- safe_path(save_dir, paste0("regularize_allHRFs_result_", sess, "_", n_use, "s.rds"))
      regularize_allHRFs_result <- readRDS(save_path)


      ###############
      ## Mean Maps  #
      ###############

      for (param in c("a1", "b1", "c")) {
        # Mean across all subjects
        plot(regularize_allHRFs_result, type = "mean_all", param = param,
             fname = safe_path("~/Documents/Github/HRF-Adaptation-paper/plots/validation/regularize_allHRFs", sess, "mean_all", paste0(sess, "_", param, "_", n_use, "s.png")),
             title = paste0("Mean of ", param, " (over all subjects) – ", sess, " ", n_use, "s"))

        # Mean over activated voxels only
        plot(regularize_allHRFs_result, type = "mean", param = param,
             fname = safe_path("~/Documents/Github/HRF-Adaptation-paper/plots/validation/regularize_allHRFs", sess, "mean", paste0(sess, "_", param, "_", n_use, "s.png")),
             title = paste0("Mean of ", param, " (over subjects with activation) – ", sess, " ", n_use, "s"))
      }

      rm(regularize_allHRFs_result)
      gc()

    }
  }
}

