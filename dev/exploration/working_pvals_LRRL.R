library(here)
library(ciftiTools)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

save_dir <- "/Volumes/LaCie/root/hrf_adaptation/validation/fit_workingHRF"
plot_dir <- "~/Documents/Github/hrf-z/dev/exploration/plots/pfvals"

sess_names <- c("tfMRI_MOTOR")   # t in 1:1
run_names  <- c("LR", "RL")      # r in 1:2

subject_idx <- 789 # Deprecated, ignore.

num_subjects <- 200

for (t in 1:1) {
  for (r in 1:2) {

    sess <- paste(sess_names[t], run_names[r], sep = "_")

    # same pattern as your working script: get subject_num from fixtures
    session_data <- readRDS(file.path(
      "~/Documents/Github/HRF-Adaptation-paper/Code/fixtures",
      paste0("session_data_", sess, ".rds")
    ))
    subject_num <- length(session_data$BOLD_files)

    save_path <- file.path(save_dir, paste0("fit_workingHRF_result_", sess, "_", subject_num, "s.rds"))

    cat("Loading ", save_path, "...\n", sep=""); flush.console()
    workingHRF_result <- readRDS(save_path)
    cat("Loaded.\n"); flush.console()

    f0 <- workingHRF_result[["subject_results"]][[1]][["glm_results"]][["Fstat_xii"]]
    n_vox <- nrow(as.matrix(f0))

    F_all <- matrix(NA_real_, nrow = n_vox, ncol = num_subjects)
    for (i in seq_len(num_subjects)) {
      fstat_xii <- workingHRF_result[["subject_results"]][[i]][["glm_results"]][["Fstat_xii"]]
      # pval_xii  <- workingHRF_result[["subject_results"]][[subject_idx]][["glm_results"]][["pvalF_xii"]]

      F_all[, i] <- as.vector(as.matrix(fstat_xii))
      stopifnot(length(as.vector(as.matrix(fstat_xii))) == n_vox)
    }

    F_avg_mat <- matrix(rowMeans(F_all, na.rm = TRUE), ncol = 1)

    F_avg_xii <- ciftiTools::newdata_xifti(f0, F_avg_mat)

    # F-stat plot
    plot(
      F_avg_xii,
      title = paste0(sess, " - ", num_subjects , "Subjects - F-statistics"),
      color_mode = "sequential",
      zlim = c(0, 3.96),
      fname = file.path(plot_dir, paste0(sess, "_avg_", num_subjects, "_subjects_Fstat.png"))
    )

    # p-value plot
    # plot(
    #   pval_xii,
    #   title = paste0(sess, " - Subject ", subject_idx, " - P-values"),
    #   color_mode = "sequential",
    #   zlim = c(0, 0.755),
    #   fname = file.path(plot_dir, paste0(sess, "_subject_", subject_idx, "_pval.png"))
    # )

    rm(workingHRF_result); gc()
  }
}


