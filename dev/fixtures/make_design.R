library(here)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z") # hrf-z is Zeshawn's branch with latest hrf package changes.
# devtools::load_all("~/Documents/Github/hrf-HRFcalc-mods") # Temporary while cleaning up hrf-z branch

session_data <- readRDS(here("dev", "fixtures", "session_data_motorlr_4s.rds"))
EVs <- session_data[["EVs_list"]][[1]]

design_1 <- make_design(
  EVs = EVs, nTime = 284, TR = 0.72, dHRF = 0,
  onset = TRUE, offset = TRUE,
  a1 = 3, b1 = 1/2, c = 1/6,
  a2 = 8, b2 = 1/2
)
