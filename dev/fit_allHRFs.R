library(here)
library(qs)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z") # hrf-z is Zeshawn's branch with latest hrf package changes.
# devtools::load_all("~/Documents/Github/hrf-HRFcalc-mods") # Temporary while cleaning up hrf-z branch

# ------------------------------------------------------------------------------
#                              HRF GRID GENERATION
# ------------------------------------------------------------------------------

my_hrf_grid <- generate_hrf_grid()


# ------------------------------------------------------------------------------
#                                 RUN fit_allHRFs
# ------------------------------------------------------------------------------

session_data <- readRDS(here("dev", "fixtures", "session_data_4s", "session_data_motor_lr_4s.rds"))

data("default_hrf_grid")
head(default_hrf_grid)

t <- 1 # Temporary till we loop

fit_allHRFs_result <- fit_allHRFs(
  BOLD = session_data$BOLD_files,
  EVs = session_data$EVs_list,
  TR = 0.72,
  hrf_grid = my_hrf_grid,
  # hrf_grid = default_hrf_grid,
  brainstructures = c('left', 'right'),
  resamp_res = 10000,
  hpf = 0.01,
  nuisance = session_data$nuisance_files,
  onsets = if(t == 3) FALSE else TRUE,    # No onsets for gambling
  offsets = if(t == 3) FALSE else TRUE,   # No offsets for gambling
  scrub = NULL,
  verbose = 2,
  n_cores = 1,
  log_dir = "dev/logs",
  work_dir = "dev/work/"
)

saveRDS(
  fit_allHRFs_result,
  file = here("dev", "fixtures", "fit_allHRFs_result_motorlr_4s.rds")
)

statuses <- sapply(fit_allHRFs_result$subject_results, function(x) x$status)
table(statuses)

# fit_allHRFs_result <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_500s.rds"))
fit_allHRFs_result <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_4s.rds"))


plot_hrf_preview(default_hrf_grid, EVs = session_data$EVs_list[[1]],
                 TR = 0.72, nT = 284, hrf_idx = 80)

plot(fit_allHRFs_result, subject = 1, hrf_idx = 12)
plot(fit_allHRFs_result, type = "hrfs")
plot(fit_allHRFs_result, type = "hrfs", tapered = FALSE)
plot(fit_allHRFs_result, type = "param_grid")
plot(fit_allHRFs_result, type = "single_hrf", hrf_idx = 40)
plot(fit_allHRFs_result, type = "single_hrf", hrf_idx = 40, tapered = FALSE)

plot(fit_allHRFs_result, type = "multiple_hrf", hrf_idx = c(1, 5, 40),
     colors = c("#2c7fb8", "#d95f02", "#7570b3"), tapered = FALSE)


# ------------------------------------------------------------------------------
#                                   DEBUG AREA
# ------------------------------------------------------------------------------
file_paths <- attr(fit_allHRFs_result, "result_paths")
work1 <- qs::qread(file_paths[[1]])

