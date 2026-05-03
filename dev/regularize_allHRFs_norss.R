library(here)
library(ciftiTools)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z")

fit_workingHRF_result <- readRDS(here("dev", "fixtures", "fit_workingHRF_result_motorlr_4s.rds"))
fit_allHRFs_result    <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_4s_norss.rds"))
session_data          <- readRDS(here("dev", "fixtures", "session_data_4s", "session_data_motor_lr_4s.rds"))

regularize_result <- regularize_allHRFs(
  workingHRF_results = fit_workingHRF_result,
  allHRF_results     = fit_allHRFs_result,
  BOLD     = session_data$BOLD_files,
  EVs      = session_data$EVs_list,
  nuisance = session_data$nuisance_files,
  TR = 0.72,
  onsets = TRUE, offsets = TRUE,
  seffects = TRUE,
  n_cores  = 4,
  log_dir  = "dev/logs",
  verbose  = 1
)

saveRDS(
  regularize_result,
  file = here("dev", "fixtures", "regularize_allHRFs_result_motorlr_4s_norss.rds")
)
