library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z")

session_data <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")
reg_result   <- readRDS("dev/fixtures/regularize_allHRFs_result_motorlr_4s.rds")
allHRF_res   <- readRDS("dev/fixtures/fit_allHRFs_result_motorlr_4s.rds")

subj <- 1L

# Both modes in one call (lookup-mode personalized via allHRF cache).
result <- fit_bestHRF(
  regularize_result = reg_result,
  BOLD_file         = session_data$BOLD_files[subj],
  EVs               = session_data$EVs_list[[subj]],
  nuisance_file     = session_data$nuisance_files[subj],
  TR                = 0.72,
  use               = c("personalized", "adapted"),
  allHRF_result     = allHRF_res,
  subject_idx       = subj,
  onsets = TRUE, offsets = TRUE,
  verbose = 1
)

names(result)
result$df

# Adapted (population-average HRF applied to this subject)
plot(result$adapted$betas,             idx = 1, title = "adapted: cue betas")
plot(result$adapted$contrasts$tstat,   idx = 1, title = "adapted: cue t-stat")
plot(result$adapted$contrasts$pval_adj, idx = 1, title = "adapted: cue p-val (BH adjusted)")

# Personalized (per-voxel best HRF from the candidate scan)
plot(result$personalized$betas,             idx = 1, title = "personalized: cue betas")
plot(result$personalized$contrasts$tstat,   idx = 1, title = "personalized: cue t-stat")
plot(result$personalized$contrasts$pval_adj, idx = 1, title = "personalized: cue p-val (BH adjusted)")
