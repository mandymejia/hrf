library(ciftiTools)
library(here)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

input_file <- here("dev", "fixtures", "hcp_data", "100206", "MOTOR_LR", "tfMRI_MOTOR_LR_Atlas_MSMAll.dtseries.nii")
output_dir <- here("dev", "fixtures", "hcp_data", "100206", "MOTOR_LR", "resampled")

resolutions <- c(2000, 1000, 500, 250, 125, 100)
results <- data.frame(Resolution = integer(), File = character(), Size_MB = numeric())

for (res in resolutions) {
  message(sprintf("Resampling to %d vertices per hemisphere...", res))

  output_name <- sprintf("resampled_%d.dtseries.nii", res)
  output_path <- file.path(output_dir, output_name)

  resampled <- resample_xifti(
    x = input_file,
    resamp_res = res,
    write_dir = output_dir,
    verbose = FALSE
  )

  resampled_file <- resampled[["cifti"]]
  file.rename(resampled_file, output_path)

  size_mb <- file.info(output_path)$size / (1024^2)

  results <- rbind(results, data.frame(
    Resolution = res,
    File = output_name,
    Size_MB = round(size_mb, 2)
  ))

  message(sprintf("Done: %s (%.2f MB)", output_name, size_mb))
}

message("Resampling Summary:")
print(results, row.names = FALSE)



results <- data.frame(Resolution = integer(), File = character(), Time_sec = numeric(), Size_MB = numeric())

for (res in resolutions) {
  file_name <- sprintf("resampled_%d.dtseries.nii", res)
  file_path <- file.path(output_dir, file_name)

  if (!file.exists(file_path)) {
    warning(sprintf("File not found: %s", file_name))
    next
  }

  size_mb <- file.info(file_path)$size / (1024^2)

  message(sprintf("Reading %s...", file_name))
  time_taken <- system.time({
    xifti <- ciftiTools::read_cifti(file_path)
  })[["elapsed"]]

  results <- rbind(results, data.frame(
    Resolution = res,
    File = file_name,
    Time_sec = round(time_taken, 2),
    Size_MB = round(size_mb, 2)
  ))

  message(sprintf("Done: %.2f sec, %.2f MB", time_taken, size_mb))
}

message("Read Performance Summary:")
print(results[order(results$Resolution, decreasing = TRUE), ], row.names = FALSE)
