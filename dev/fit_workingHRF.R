library(here)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z") # hrf-z is Zeshawn's branch with latest hrf package changes.
# devtools::load_all("~/Documents/Github/hrf-HRFcalc-mods") # Temporary while cleaning up hrf-z branch


session_data <- readRDS(here("dev", "fixtures", "session_data_motorlr_4s.rds"))


fit_workingHRF_result <- fit_workingHRF(
  BOLD = session_data$BOLD_files,
  EVs = session_data$EVs_list,
  TR = 0.72,
  brainstructures = c('left', 'right'),
  resamp_res = 10000, # 10,000
  hpf = 0.01,
  nuisance = session_data$nuisance_files,
  onsets = TRUE,
  offsets = TRUE,
  hrf_params = list(a1 = 6, b1 = 1, c = 1/6, a2 = 16, b2 = 1),
  derivatives = FALSE,
  scrub = NULL,
  alpha = 0.01,
  verbose = 2,
  n_cores = 2,
  log_dir = "dev/logs"
)

statuses <- sapply(fit_workingHRF_result[["subject_results"]], function(x) x[["status"]])
table(statuses)
