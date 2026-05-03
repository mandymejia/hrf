# Same as dev/fit_allHRFs.R but with save_rss = FALSE so regularize_allHRFs
# will go into refit mode. Used to test refit-mode parallelization locally.

library(here)
library(qs2)
library(ciftiTools)

ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z")

my_hrf_grid <- generate_hrf_grid()

session_data <- readRDS(here("dev", "fixtures", "session_data_4s", "session_data_motor_lr_4s.rds"))

t <- 1

fit_allHRFs_result <- fit_allHRFs(
  BOLD = session_data$BOLD_files,
  EVs = session_data$EVs_list,
  TR = 0.72,
  hrf_grid = my_hrf_grid,
  brainstructures = c('left', 'right'),
  resamp_res = 10000,
  hpf = 0.01,
  nuisance = session_data$nuisance_files,
  onsets = if(t == 3) FALSE else TRUE,
  offsets = if(t == 3) FALSE else TRUE,
  scrub = NULL,
  save_rss = FALSE,
  verbose = 2,
  n_cores = 1,
  log_dir = "dev/logs",
  work_dir = "dev/work/"
)

saveRDS(
  fit_allHRFs_result,
  file = here("dev", "fixtures", "fit_allHRFs_result_motorlr_4s_norss.rds")
)

statuses <- sapply(fit_allHRFs_result$subject_results, function(x) x$status)
table(statuses)
