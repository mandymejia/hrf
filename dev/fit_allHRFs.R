library(here)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z") # hrf-z is Zeshawn's branch with latest hrf package changes.
# devtools::load_all("~/Documents/Github/hrf-HRFcalc-mods") # Temporary while cleaning up hrf-z branch


session_data <- readRDS(here("dev", "fixtures", "session_data_motorlr_4s.rds"))

data("default_hrf_grid")
head(default_hrf_grid)

t <- 1 # Temporary till we loop

fit_allHRFs_result <- fit_allHRFs(
  BOLD = session_data$BOLD_files,
  EVs = session_data$EVs_list,
  TR = 0.72,
  hrf_grid = generate_default_hrf_grid,
  brainstructures = c('left', 'right'),
  resamp_res = 10000,
  hpf = 0.01,
  nuisance = session_data$nuisance_files,
  onsets = if(t == 3) FALSE else TRUE,    # No onsets for gambling
  offsets = if(t == 3) FALSE else TRUE,   # No offsets for gambling
  scrub = NULL,
  verbose = 2,
  n_cores = 1,
  log_dir = "dev/logs"
)

saveRDS(
  fit_allHRFs_result,
  file = here("dev", "fixtures", "fit_allHRFs_result_motorlr_4s.rds")
)

statuses <- sapply(fit_workingHRF_result$subject_results, function(x) x$status)
table(statuses)
