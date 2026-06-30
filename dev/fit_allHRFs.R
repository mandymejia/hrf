library(here)
library(qs2)
library(ciftiTools) # make hrf package load ciftitools

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

combo <- fit_allHRFs(
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
  save_rss = TRUE,
  verbose = 2,
  n_cores = 4,
  working_hrf = list(a1 = 6, b1 = 1, c = 1/6, a2 = 16, b2 = 1),
  derivatives = FALSE,
  alpha = 0.01,
  min_active_subjects = 20,
  log_dir = "dev/logs",
  work_dir = "dev/work/"
)

class(combo)              # "hrfs"
names(combo)              # fit_workingHRF / fit_allHRFs / regularize_allHRFs / session_info
str(combo$session_info)

statuses <- sapply(combo$fit_allHRFs$subject_results, function(x) x$status)
table(statuses)


# ------------------------------------------------------------------------------
#                              PLOT.HRFS DISPATCH
# ------------------------------------------------------------------------------
# plot(combo, type = ...) routes to the right sub-object. See ?plot.hrfs.

# workingHRF types
plot(combo, type = "proportion")
plot(combo, type = "mask")
plot(combo, type = "binary", threshold = 0.2)

# allHRFs types
plot(combo, type = "hrfs")
plot(combo, type = "hrfs", tapered = FALSE)
plot(combo, type = "param_grid")
plot(combo, type = "single_hrf", hrf_idx = 24)
plot(combo, type = "single_hrf", hrf_idx = 40, tapered = FALSE)
plot(combo, type = "multiple_hrf", hrf_idx = c(1, 5, 40),
     colors = c("#2c7fb8", "#d95f02", "#7570b3"), tapered = TRUE)
plot(combo, type = "multiple_hrf", hrf_idx = c(1, 5, 40),
     colors = c("#2c7fb8", "#d95f02", "#7570b3"), tapered = FALSE)

# Design plots — ambiguous so disambiguate with which=
plot(combo, type = "design", which = "working", subject = 1)
plot(combo, type = "design", which = "all",     subject = 1, hrf_idx = 19)
plot(combo, type = "design", which = "all",     subject = 1, hrf_idx = 21)

# regularize types
plot(combo, type = "pop_avg",            param = "a1")
plot(combo, type = "pop_avg_continuous", param = "a1")
plot(combo, type = "mean_all",           param = "a1")
plot(combo, type = "mean",               param = "a1")
plot(combo, type = "param_heatmap")

# hrf_grid preview (independent of combo)
plot_hrf_preview(default_hrf_grid, EVs = session_data$EVs_list[[1]],
                 TR = 0.72, nT = 284, hrf_idx = 80)


# ------------------------------------------------------------------------------
#                                   DEBUG AREA
# ------------------------------------------------------------------------------
file_paths <- attr(combo$fit_allHRFs, "result_paths")
work1 <- qs2::qs_read(file_paths[[1]])
