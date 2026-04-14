library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z")

session_data <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")

# Canonical HRF baseline: every voxel gets a1=6, b1=1, c=1/6
nV <- 18792
canonical_reg <- list(
  pop_avg = data.frame(
    voxel = 1:nV,
    a1_snapped = rep(6, nV),
    b1_snapped = rep(1, nV),
    c_snapped = rep(1/6, nV)
  ),
  subject_results = NULL
)

result_canonical <- fit_bestHRF(
  canonical_reg,
  BOLD_file = session_data$BOLD_files[1],
  EVs = session_data$EVs_list[[1]],
  nuisance_file = session_data$nuisance_files[1],
  TR = 0.72, onsets = TRUE, offsets = TRUE, verbose = 1
)

# Population-adapted HRF
reg_result <- readRDS("dev/fixtures/regularize_result_motorlr_1125s_slim.rds")
reg_result$subject_results <- NULL

result_adapted <- fit_bestHRF(
  reg_result,
  BOLD_file = session_data$BOLD_files[1],
  EVs = session_data$EVs_list[[1]],
  nuisance_file = session_data$nuisance_files[1],
  TR = 0.72, onsets = TRUE, offsets = TRUE, verbose = 1
)

# Compare cue betas side by side
plot(result_canonical$contrasts$tstat, idx = 1, zlim = c(-8, 8), title = "Canonical HRF - cue t-stat")
plot(result_adapted$contrasts$tstat, idx = 1, zlim = c(-8, 8), title = "Adapted HRF - cue t-stat")

plot(result_canonical$betas, idx = 1, title = "Canonical HRF - cue betas")
plot(result_adapted$betas, idx = 1, title = "Pop Adapted HRF - cue betas")

# Subject-adapted HRF (seffects): use candidate map for subject 1
reg_subj <- readRDS("dev/fixtures/regularize_result_motorlr_1125s_slim.rds")

result_subject <- fit_bestHRF(
  reg_subj,
  BOLD_file = session_data$BOLD_files[1],
  EVs = session_data$EVs_list[[1]],
  nuisance_file = session_data$nuisance_files[1],
  TR = 0.72, subject_idx = 1,
  onsets = TRUE, offsets = TRUE, verbose = 1
)

plot(result_subject$contrasts$tstat, idx = 1, zlim = c(-8, 8), title = "Subject Adapted HRF - cue t-stat")
plot(result_subject$betas, idx = 1, title = "Subject Adapted HRF - cue betas")
