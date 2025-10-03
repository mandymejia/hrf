library(here)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z") # hrf-z is Zeshawn's branch with latest hrf package changes.
# devtools::load_all("~/Documents/Github/hrf-HRFcalc-mods") # Temporary while cleaning up hrf-z branch

# workingHRF_results <- readRDS(here("dev", "fixtures", "fit_workingHRF_result_motorlr_4s.rds"))
# allHRF_results <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_4s.rds"))

# workingHRF_results <- readRDS("/Users/zeshawnzahid/Downloads/19s_fit_workingHRF_result.rds")
# allHRF_results <- readRDS("/Users/zeshawnzahid/Downloads/19s_fit_allHRFs_result.rds")

workingHRF_results <- readRDS(here("dev", "fixtures", "fit_workingHRF_result_motorlr_500s.rds"))
allHRF_results <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_500s.rds"))

regularize_allHRFs_result <- regularize_allHRFs(
  workingHRF_results,
  allHRF_results,
  rounding = TRUE#,
  # log = TRUE,
  # TRUNCATE = TRUE,
  # mask = TRUE,
)

saveRDS(
  regularize_allHRFs_result,
  file = here("dev", "fixtures", "regularize_allHRFs_result_motorlr_500s.rds")
)


plot(regularize_allHRFs_result, type = "mask", fname = here("dev", "test_plots", "regularize_allHRFs", "mask_prop_s"))
