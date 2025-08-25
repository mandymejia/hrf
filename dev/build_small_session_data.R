library(here)

subjects <- c("100206", "100307", "100408", "100610")
fixture_root <- here("dev", "fixtures")

# Save each subject's EVs list into its folder
for (i in seq_along(subjects)) {
  subj_dir <- file.path(fixture_root, subjects[i])
  dir.create(subj_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(session_data$EVs_list[[i]], file.path(subj_dir, "EVs_list.rds"))
}

# Build session_data object with proper fixture paths
session_data <- list(
  BOLD_files = file.path(fixture_root, subjects, "tfMRI_MOTOR_LR_Atlas_MSMAll.dtseries.nii"),
  nuisance_files = file.path(fixture_root, subjects, "Movement_Regressors.txt"),
  EVs_list = lapply(subjects, function(subj) {
    readRDS(file.path(fixture_root, subj, "EVs_list.rds"))
  })
)

# Save to fixtures as one combined RDS
saveRDS(session_data, file.path(fixture_root, "session_data_motorlr_4s.rds"))
