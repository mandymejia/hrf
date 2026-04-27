# Quick HPC sample: run fit_bestHRF on MOTOR_LR (max subjects, population mode)
# for subject 1 and save the working + population beta plots as PNGs.
#
# See dev/hpc/README.md for setup (RStudio launch, dependency install).

library(ciftiTools)
devtools::load_all("~/Documents/Github/hrf-z")

ciftiTools.setOption('wb_path', '~/workbench/bin_rh_linux64/wb_command')

# ** Paths (mirror the pipeline's output layout) **
run_name        <- "validation_final"
output_root     <- file.path("/N/project/hrf_adaptation", run_name)
sess            <- "tfMRI_MOTOR_LR"
n_use           <- 1080  # max subjects for MOTOR_LR

regularize_path <- file.path(output_root, "regularize_allHRFs",
  paste0("regularize_allHRFs_result_", sess, "_", n_use, "s.rds"))
session_path    <- file.path("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures",
  paste0("session_data_", sess, ".rds"))

plot_dir <- file.path("~/Documents/Github/hrf-z/dev/hpc/plots", run_name, sess)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ** Load **
cat("Loading regularize result...\n")
regularize_result <- readRDS(regularize_path)

cat("Loading session_data...\n")
session_data <- readRDS(session_path)

# ** Run fit_bestHRF for subject 1 in population mode **
cat("Running fit_bestHRF (population mode)...\n")
result <- fit_bestHRF(
  regularize_result = regularize_result,
  BOLD_file = session_data$BOLD_files[1],
  EVs       = session_data$EVs_list[[1]],
  nuisance_file = session_data$nuisance_files[1],
  TR = 0.72,
  use = "population",
  onsets  = TRUE,
  offsets = TRUE,
  verbose = 1
)

# ** Plot betas — save to PNG so this runs headless on HPC **
cat("Saving plots to", plot_dir, "\n")
plot(result$working$betas,    idx = 1,
     fname = file.path(plot_dir, "subject1_working_betas.png"),
     title = "Subject 1 - Working HRF betas (regressor 1)")

plot(result$population$betas, idx = 1,
     fname = file.path(plot_dir, "subject1_population_betas.png"),
     title = "Subject 1 - Population HRF betas (regressor 1)")

cat("Done.\n")
