# Quick HPC sample: same as fit_bestHRF_sample.R but with a custom contrast
# (hands - feet on MOTOR). Plots the contrast t-stat instead of plain betas.
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

# ** Build the contrast matrix **
# MOTOR EV order (HCP convention): cue, lf, lh, rf, rh, t  (6 main regressors)
# With onsets=TRUE & offsets=TRUE, additional onset/offset regressors come
# after the main 6. Contrast columns past the 6 main ones stay 0.
ev_names <- names(session_data$EVs_list[[1]])
n_main   <- length(ev_names)
n_task   <- n_main + 2  # main EVs + 1 shared onset + 1 shared offset
A <- matrix(0, nrow = 1, ncol = n_task)
colnames(A) <- c(ev_names, "onset", "offset")
rownames(A) <- "hands_minus_feet"
A[1, "left_hand"]  <- 0.5
A[1, "right_hand"] <- 0.5
A[1, "left_foot"]  <- -0.5
A[1, "right_foot"] <- -0.5

cat("Contrast matrix:\n"); print(A)

# ** Run fit_bestHRF for subject 1 in population mode **
cat("Running fit_bestHRF (population mode)...\n")
result <- fit_bestHRF(
  regularize_result = regularize_result,
  BOLD_file = session_data$BOLD_files[1],
  EVs       = session_data$EVs_list[[1]],
  nuisance_file = session_data$nuisance_files[1],
  TR = 0.72,
  use = "population",
  contrasts = A,
  onsets  = TRUE,
  offsets = TRUE,
  verbose = 1
)

# ** Plot contrast t-stats for the hands-minus-feet contrast (idx = 1) **
cat("Saving plots to", plot_dir, "\n")
plot(result$working$contrasts$tstat,    idx = 1,
     fname = file.path(plot_dir, "subject1_working_hands_minus_feet_tstat.png"),
     title = "Subject 1 - Working HRF - hands - feet (t-stat)")

plot(result$population$contrasts$tstat, idx = 1,
     fname = file.path(plot_dir, "subject1_population_hands_minus_feet_tstat.png"),
     title = "Subject 1 - Population HRF - hands - feet (t-stat)")

cat("Done.\n")
