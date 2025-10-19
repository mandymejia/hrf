library(here)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z") # hrf-z is Zeshawn's branch with latest hrf package changes.
# devtools::load_all("~/Documents/Github/hrf-HRFcalc-mods") # Temporary while cleaning up hrf-z branch

workingHRF_results <- readRDS(here("dev", "fixtures", "fit_workingHRF_result_motorlr_4s.rds"))
allHRF_results <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_4s.rds"))

# workingHRF_results <- readRDS("/Users/zeshawnzahid/Downloads/19s_fit_workingHRF_result.rds")
# allHRF_results <- readRDS("/Users/zeshawnzahid/Downloads/19s_fit_allHRFs_result.rds")

workingHRF_results <- readRDS(here("dev", "fixtures", "fit_workingHRF_result_motorlr_500s.rds"))
allHRF_results <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_500s.rds"))

allHRF_results <- readRDS("~/Downloads/fit_allHRFs_result_tfMRI_MOTOR_LR_1000s.rds")
workingHRF_results <- readRDS("~/Downloads/fit_workingHRF_result_tfMRI_MOTOR_LR_1000s.rds")

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

regularize_allHRFs_result <- readRDS(here("dev", "fixtures", "regularize_allHRFs_result_motorlr_4s.rds"))




plot(regularize_allHRFs_result, type = "variance", param = "a1", model = "A", method = "OLS", fname = here("dev", "test_plots", "regularize_allHRFs", "OLS_A_a1"))
plot(regularize_allHRFs_result, type = "variance", param = "a1", model = "A", method = "WLS", fname = here("dev", "test_plots", "regularize_allHRFs", "WLS_A_a1"))

plot(regularize_allHRFs_result, type = "mean", param = "b1", fname = here("dev", "test_plots", "regularize_allHRFs", "mean_b1.png"))

plot(regularize_allHRFs_result,  type = "mean_all",  param = "b1",  fname = here("dev", "test_plots", "regularize_allHRFs", "mean_b1_filtered.png"))

# DEBUG AREA

var_range_A <- range(regularize_allHRFs_result[["regularized_params"]][["a1"]][["results_OLS"]][["variance_A"]], na.rm = TRUE)
var_range_B <- range(regularize_allHRFs_result[["regularized_params"]][["a1"]][["results_OLS"]][["variance_B"]], na.rm = TRUE)
range(c(var_range_A, var_range_B))


plot(regularize_allHRFs_result[["regularized_params"]][["a1"]][["residual_variance_xii_A_OLS"]],
     title = 'Variance of Model A Predictions for a1 using OLS method',
     color_mode = 'sequential',
     zlim = range(c(regularize_allHRFs_result[["regularized_params"]][["a1"]][["residual_variance_xii_A_OLS"]]$data,
                    regularize_allHRFs_result[["regularized_params"]][["a1"]][["residual_variance_xii_B_OLS"]]$data),
                  na.rm = TRUE),
     fname = NULL,
     legend_fname = NULL)


package_best_params_df <- readRDS("~/Downloads/package_best_params_df.rds")
research_best_params_df <- readRDS("~/Downloads/research_best_params_df.rds")
