library(here)
library(ciftiTools)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z")

fit_workingHRF_result <- readRDS(here("dev", "fixtures", "fit_workingHRF_result_motorlr_4s.rds"))
fit_allHRFs_result    <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_4s.rds"))
session_data          <- readRDS(here("dev", "fixtures", "session_data_4s", "session_data_motor_lr_4s.rds"))

regularize_result <- regularize_allHRFs(
  workingHRF_results = fit_workingHRF_result,
  allHRF_results     = fit_allHRFs_result,
  seffects           = TRUE,
  verbose            = 1
)

saveRDS(
  regularize_result,
  file = here("dev", "fixtures", "regularize_allHRFs_result_motorlr_4s.rds")
)


# ------------------------------------------------------------------------------
#                                   PLOTS
# ------------------------------------------------------------------------------
# Heatmap: frequency of best per-subject params on the HRF grid.
plot(regularize_result, type = "param_heatmap")
plot(regularize_result, type = "param_heatmap", mask = TRUE)

# Brain-surface population-average HRF parameter maps.
# (xii template auto-attached by regularize_allHRFs; no template arg needed.)
plot(regularize_result, type = "pop_avg", param = "a1")
plot(regularize_result, type = "pop_avg", param = "b1")
plot(regularize_result, type = "pop_avg", param = "c")

