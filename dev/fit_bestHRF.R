library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z")

session_data <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")
reg_result <- readRDS("dev/fixtures/regularize_result_motorlr_1125s_slim.rds")
reg_result$subject_results <- NULL

result <- fit_bestHRF(
  reg_result,
  BOLD_file = session_data$BOLD_files[1],
  EVs = session_data$EVs_list[[1]],
  nuisance_file = session_data$nuisance_files[1],
  TR = 0.72,
  onsets = TRUE, offsets = TRUE,
  verbose = 1,
)

names(result)
result$df

plot(result$betas, idx = 1, title = "cue betas")
plot(result$contrasts$tstat, idx = 1, title = "cue t-stat")
plot(result$contrasts$pval_adj, idx = 1, title = "cue p-val (BH adjusted)")


# AVG on 1000 subjects
# vs
# AVG on 100 subjects
