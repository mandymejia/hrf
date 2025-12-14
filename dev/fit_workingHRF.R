library(here)
library(tictoc)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z") # hrf-z is Zeshawn's branch with latest hrf package changes.
# devtools::load_all("~/Documents/Github/hrf-HRFcalc-mods") # Temporary while cleaning up hrf-z branch


session_data <- readRDS(here("dev", "fixtures", "session_data_4s", "session_data_motor_lr_4s.rds"))

t <- 1 # Temporary till we loop

tic("fit_workingHRF done in:")
fit_workingHRF_result <- fit_workingHRF(
  BOLD = session_data$BOLD_files,
  EVs = session_data$EVs_list,
  TR = 0.72,
  brainstructures = c('left', 'right'),
  resamp_res = 10000, # 10,000
  hpf = 0.01,
  nuisance = session_data$nuisance_files,
  onsets = if(t == 3) FALSE else TRUE,    # No onsets for gambling
  offsets = if(t == 3) FALSE else TRUE,   # No offsets for gambling
  hrf_params = list(a1 = 6, b1 = 1, c = 1/6, a2 = 16, b2 = 1),
  derivatives = FALSE,
  scrub = NULL,
  alpha = 0.01,
  verbose = 2,
  n_cores = 1,
  min_active_subjects = 20,
  log_dir = "dev/logs"
)
toc()

saveRDS(
  fit_workingHRF_result,
  file = here("dev", "fixtures", "fit_workingHRF_result_motorlr_4s.rds")
)

statuses <- sapply(fit_workingHRF_result[["subject_results"]], function(x) x[["status"]])
table(statuses)



fit_workingHRF_result <- readRDS(here("dev", "fixtures", "fit_workingHRF_result_motorlr_4s.rds"))
fit_workingHRF_result <- readRDS(here("dev", "fixtures", "fit_workingHRF_result_motorlr_500s.rds"))

plot(fit_workingHRF_result, subject = 2)

plot(fit_workingHRF_result, type = "mask", fname = here("dev", "test_plots", "fit_workingHRF", "mask_prop"))

plot(fit_workingHRF_result,
     type = "proportion",
     fname = here("dev", "test_plots", "proportion_0.01"))
plot(fit_workingHRF_result,
     type = "proportion",
     alpha = 0.001,
     fname = here("dev", "test_plots", "proportion_0.001"))



plot(fit_workingHRF_result,
     type = "binary",
     fname = here("dev", "test_plots", "fit_workingHRF", "binary_0.01_10pct_test"))
plot(fit_workingHRF_result,
     type = "binary", threshold = 0.5,
     fname = here("dev", "test_plots", "fit_workingHRF", "binary_0.01_50pct"))
plot(fit_workingHRF_result,
     type = "binary", alpha = 0.001,
     fname = here("dev", "test_plots", "fit_workingHRF", "binary_0.001_10pct"))
plot(fit_workingHRF_result,
     type = "binary", alpha = 0.001, threshold = 0.5,
     fname = here("dev", "test_plots", "fit_workingHRF", "binary_0.001_50pct"))
