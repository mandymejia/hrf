library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z")

session_data <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")
reg_result   <- readRDS("dev/fixtures/regularize_allHRFs_result_motorlr_4s.rds")
allHRF_res   <- readRDS("dev/fixtures/fit_allHRFs_result_motorlr_4s.rds")

# Synthesize a combo for the probe. In real use, this comes straight from
# fit_allHRFs() and you just pass `combo` directly.
combo <- list(
  fit_workingHRF     = NULL,
  fit_allHRFs        = allHRF_res,
  regularize_allHRFs = reg_result
)
class(combo) <- "hrfs"

subj <- 1L

# Both modes in one call (lookup-mode personalized via subject_idx).
result <- fit_bestHRF(
  combo         = combo,
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

plot(result$adapted$betas,             idx = 1, title = "adapted: cue betas")
plot(result$adapted$contrasts$tstat,   idx = 1, title = "adapted: cue t-stat")
plot(result$adapted$contrasts$pval_adj, idx = 1, title = "adapted: cue p-val (BH adjusted)")

plot(result$personalized$betas,             idx = 1, title = "personalized: cue betas")
plot(result$personalized$contrasts$tstat,   idx = 1, title = "personalized: cue t-stat")
plot(result$personalized$contrasts$pval_adj, idx = 1, title = "personalized: cue p-val (BH adjusted)")
