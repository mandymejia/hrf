# Probe: run fit_bestHRF on the 4-subject MOTOR_LR fixture.
# This is what real usage looks like — fit_allHRFs() returns the combo,
# fit_bestHRF() takes it directly. No object construction needed.

library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z")

session_data       <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")
fit_allHRFs_result <- readRDS("dev/fixtures/fit_allHRFs_combo_motorlr_4s.rds")

subj <- 1L

result <- fit_bestHRF(
  fit_allHRFs_result,
  BOLD_file     = session_data$BOLD_files[subj],
  EVs           = session_data$EVs_list[[subj]],
  nuisance_file = session_data$nuisance_files[subj],
  TR            = 0.72,
  use           = c("personalized", "adapted"),
  subject_idx   = subj,
  onsets = TRUE, offsets = TRUE,
  verbose = 1
)

names(result)
result$df

plot(result$adapted$betas,              idx = 1, title = "adapted: cue betas")
plot(result$adapted$contrasts$tstat,    idx = 1, title = "adapted: cue t-stat")
plot(result$adapted$contrasts$pval_adj, idx = 1, title = "adapted: cue p-val (BH adjusted)")

plot(result$personalized$betas,              idx = 1, title = "personalized: cue betas")
plot(result$personalized$contrasts$tstat,    idx = 1, title = "personalized: cue t-stat")
plot(result$personalized$contrasts$pval_adj, idx = 1, title = "personalized: cue p-val (BH adjusted)")
