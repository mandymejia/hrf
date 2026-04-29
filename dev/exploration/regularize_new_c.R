source("~/Documents/Github/HRF-Adaptation-paper/Code/0_analysis.R")

devtools::load_all("~/Documents/Github/hrf-z")
ciftiTools.setOption('wb_path', if (Sys.info()[['sysname']] == 'Darwin') '/Applications/workbench/bin_macosxub/wb_command' else '~/workbench/bin_rh_linux64/wb_command')

# Base directories - change these to switch between validation and smoothed
# Option 1: validation_bonferroni_0.05
# Option 2: retest_bonferroni_0.05_smoothed
base_dir <- "/N/project/hrf_adaptation/validation_bonferroni_0.05"
# base_dir <- "/N/project/hrf_adaptation/retest_bonferroni_0.05_smoothed"

load_dir_workingHRF <- file.path(base_dir, "fit_workingHRF")
load_dir_allHRFs <- file.path(base_dir, "fit_allHRFs")
save_dir <- file.path(base_dir, "regularize_allHRFs_winnerc")

sess_names <- c('tfMRI_MOTOR','tfMRI_WM','tfMRI_GAMBLING')
run_names <- c('LR','RL')

for(t in 1:3){
  for(r in 1:2){
    sess <- paste(sess_names[t], run_names[r], sep='_')
    session_data <- readRDS(file.path('~/Documents/Github/HRF-Adaptation-paper/Code/fixtures', paste0('session_data_', sess, '.rds')))

    # Only use max subjects
    n_use <- length(session_data$BOLD_files)

    cat("\n========================================\n")
    cat("Processing:", sess, "with", n_use, "subjects\n")
    cat("========================================\n")

    # Load results
    workingHRF_file <- file.path(load_dir_workingHRF, paste0("fit_workingHRF_result_", sess, "_", n_use, "s.rds"))
    allHRF_file <- file.path(load_dir_allHRFs, paste0("fit_allHRFs_result_", sess, "_", n_use, "s.rds"))

    cat("Loading:", workingHRF_file, "\n")
    workingHRF_results <- readRDS(workingHRF_file)

    cat("Loading:", allHRF_file, "\n")
    allHRF_results <- readRDS(allHRF_file)

    # Run regularization with new c winner-take-all logic
    cat("Running regularize_allHRFs...\n")
    tic(sprintf("regularize_allHRFs | session=%s | n_subj=%d", sess, n_use))
    regularize_allHRFs_result <- regularize_allHRFs(workingHRF_results, allHRF_results)
    toc()

    # Save results
    save_path <- safe_path(save_dir, paste0("regularize_allHRFs_result_", sess, "_", n_use, "s.rds"))
    saveRDS(regularize_allHRFs_result, save_path)
    cat("Saved to:", save_path, "\n")

    rm(workingHRF_results, allHRF_results, regularize_allHRFs_result); gc()
  }
}

cat("\n========================================\n")
cat("All sessions completed!\n")
cat("Results saved to:", save_dir, "\n")
cat("========================================\n")
