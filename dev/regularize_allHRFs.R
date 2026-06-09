library(here)
library(ciftiTools)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z")

fit_workingHRF_result <- readRDS(here("dev", "fixtures", "fit_workingHRF_result_motorlr_4s.rds"))
fit_allHRFs_result    <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_4s.rds"))

regularize_result <- regularize_allHRFs(
  workingHRF_results = fit_workingHRF_result,
  allHRF_results     = fit_allHRFs_result,
  verbose            = 1
)

saveRDS(
  regularize_result,
  file = here("dev", "fixtures", "regularize_allHRFs_result_motorlr_4s.rds")
)
