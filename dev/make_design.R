library(here)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tictoc)
library(ciftiTools)

devtools::load_all("~/Documents/Github/hrf-z")
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')

task_names <- c("motor", "wm", "gambling")
nT_lookup <- c(284, 405, 253)  # motor, wm, gambling respectively



for (t in 1) {  # 1 = motor, 2 = wm, 3 = gambling

  session_data <- readRDS(here(
    "dev", "fixtures", "session_data_4s",
    paste0("session_data_", task_names[t], "_lr_4s.rds")
  ))

  EVs <- session_data[["EVs_list"]][[1]]
  str(EVs)
  nTime <- nT_lookup[t]

  onsets <- c("left_hand",  "right_hand", "left_foot" , "right_foot", "tongue")
  offsets <-  c("left_hand",  "right_hand", "left_foot" , "right_foot", "tongue")

  # onsets <- names(EVs)
  # offsets <- names(EVs)

  cat("\n****Running make_design for", toupper(task_names[t]),"_LR .... \n")

  # tic("***make_design when not optimized")
  # design_not_optimized <- make_design(
  #   EVs = EVs, nTime = nTime, TR = 0.72, dHRF = 0,
  #   onset = if (t == 3) FALSE else TRUE,
  #   offset = if (t == 3) FALSE else TRUE,
  #   # onset = onsets, offset = offsets,
  #   # optimize = FALSE
  # )
  # toc()

  tic("***make_design when optimized")
  design_optimized <- make_design(
    EVs = EVs, nTime = nTime, TR = 0.72, dHRF = 2,
    onset = if (t == 3) FALSE else TRUE,
    offset = if (t == 3) FALSE else TRUE,
    # onset = onsets, offset = offsets,
    # optimize = TRUE
  )
  toc()

  # plot(design_not_optimized)
  plot(design_optimized)
}









