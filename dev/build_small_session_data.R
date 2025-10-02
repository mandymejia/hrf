library(here)

# You need to first download 2 files per subject/task/direction
# the .nii file and the movement_regressors
# Then the 3rd file EVs_list.rd is ripped froo the research code's session_data

subjects <- c("100206", "100307", "100408", "100610")
fixture_root <- here("dev", "fixtures", "hcp_data")
save_root <- here("dev", "fixtures", "session_data_4s")

tasks <- c("MOTOR_LR", "WM_LR", "GAMBLING_LR")

for (task in tasks) {
  session_data <- list(
    BOLD_files = file.path(fixture_root, subjects, task, paste0("tfMRI_", task, "_Atlas_MSMAll.dtseries.nii")),
    nuisance_files = file.path(fixture_root, subjects, task, "Movement_Regressors.txt"),
    EVs_list = lapply(subjects, function(subj) {
      readRDS(file.path(fixture_root, subj, task, "EVs_list.rds"))
    })
  )

  saveRDS(session_data, file.path(save_root, paste0("session_data_", tolower(task), "_4s.rds")))
  cat("✓ Saved session_data for", task, "\n")
}
