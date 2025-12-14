library(here)
library(qs)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z")

# Extract rounded params from earlier compute
base_dir <- "/Volumes/LaCie/root/hrf_adaptation"
regularize_allHRFs_result <- readRDS(file.path(base_dir, "validation/regularize_allHRFs/regularize_allHRFs_result_tfMRI_MOTOR_LR_1080s.rds"))
saveRDS(
  regularize_allHRFs_result[["rounded_params"]],
  file = here("dev", "fixtures", "regularized_params_motor_1080s.rds")
)



fit_allHRFs_result <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_4s.rds"))
regularized_params <- readRDS(here("dev", "fixtures", "regularized_params_motor_1080s.rds"))

model <- "rounded_intercept_only"
subjects <- 1080
resolution <- 18792
subject <- 1

chunk_size <- length(regularized_params[[model]]$a1) / subjects

subject_params_a1 <- regularized_params[[model]]$a1[((subject - 1) * chunk_size + 1):(subject * chunk_size)]
subject_params_b1 <- regularized_params[[model]]$b1[((subject - 1) * chunk_size + 1):(subject * chunk_size)]
subject_params_c  <- regularized_params[[model]]$c[((subject - 1) * chunk_size + 1):(subject * chunk_size)]

# Function to find matching HRF index in default_hrf_grid
find_hrf_index <- function(a1, b1, c, hrf_grid = default_hrf_grid) {
  matches <- which(
    hrf_grid$a1 == a1 &
      hrf_grid$b1 == b1 &
      abs(hrf_grid$c - c) < 1e-6  # Use tolerance for floating point comparison
  )
  if (length(matches) == 0) {
    stop("No matching HRF found for a1=", a1, ", b1=", b1, ", c=", c)
  }
  return(matches[1])
}

# Concatenate left and right hemisphere F-stats into one matrix per HRF
n_hrfs <- length(fit_allHRFs_result[["subject_results"]][[subject]][["glm_result"]][["GLMs"]])
fstat_matrix <- matrix(nrow = resolution, ncol = n_hrfs)

for (hrf_idx in 1:n_hrfs) {
  fstat_left <- fit_allHRFs_result[["subject_results"]][[subject]][["glm_result"]][["GLMs"]][[hrf_idx]][["Fstat"]][["data"]][["cortex_left"]][, 1]
  fstat_right <- fit_allHRFs_result[["subject_results"]][[subject]][["glm_result"]][["GLMs"]][[hrf_idx]][["Fstat"]][["data"]][["cortex_right"]][, 1]

  fstat_matrix[, hrf_idx] <- c(fstat_left, fstat_right)
}

# Now extract the best F-stat for each voxel
best_fstat <- numeric(resolution)

for (voxel in 1:resolution) {
  # Get regularized params for this voxel
  a1 <- subject_params_a1[voxel]
  b1 <- subject_params_b1[voxel]
  c <- subject_params_c[voxel]

  # Find corresponding HRF index
  hrf_idx <- find_hrf_index(a1, b1, c)

  # Extract F-stat from the appropriate HRF
  best_fstat[voxel] <- fstat_matrix[voxel, hrf_idx]
}

# Verify
cat("Total voxels with best F-stats:", length(best_fstat), "\n")

template_xii <- fit_allHRFs_result[["subject_results"]][[1]][["glm_result"]][["bestmodel_xii"]]
template_dscalar <- ciftiTools::convert_xifti(template_xii, to = "dscalar")

# Now create the new xifti with best F-stats
best_fstat_xifti <- ciftiTools::newdata_xifti(xifti = template_dscalar,newdata = best_fstat,newnames = "regularized_Fstat")
plot(best_fstat_xifti)


###############################################################################
###### REGRESSORS                                                             #
###############################################################################
# Configuration
model <- "rounded_intercept_only"
subjects <- 1080
resolution <- 18792
subject <- 2

# Extract subject-specific regularized parameters
chunk_size <- length(regularized_params[[model]]$a1) / subjects
subject_params_a1 <- regularized_params[[model]]$a1[((subject - 1) * chunk_size + 1):(subject * chunk_size)]
subject_params_b1 <- regularized_params[[model]]$b1[((subject - 1) * chunk_size + 1):(subject * chunk_size)]
subject_params_c  <- regularized_params[[model]]$c[((subject - 1) * chunk_size + 1):(subject * chunk_size)]

# Function to find matching HRF index in hrf_grid
find_hrf_index <- function(a1, b1, c, hrf_grid) {
  matches <- which(
    hrf_grid$a1 == a1 &
      hrf_grid$b1 == b1 &
      abs(hrf_grid$c - c) < 1e-6  # Use tolerance for floating point comparison
  )
  if (length(matches) == 0) {
    stop("No matching HRF found for a1=", a1, ", b1=", b1, ", c=", c)
  }
  return(matches[1])
}

# Get GLM result and determine number of HRFs and regressors
glm_result <- fit_allHRFs_result[["subject_results"]][[subject]][["glm_result"]]
hrf_grid <- fit_allHRFs_result$hrf_grid
n_hrfs <- length(glm_result[["GLMs"]])

# Determine number of regressors from the first GLM's betas
first_betas_xii <- glm_result[["GLMs"]][[1]][["betas"]]
n_regressors <- ncol(as.matrix(first_betas_xii))

cat("Configuration:\n")
cat("  Subject:", subject, "\n")
cat("  Resolution:", resolution, "\n")
cat("  Number of HRFs:", n_hrfs, "\n")
cat("  Number of regressors:", n_regressors, "\n")

# Get regressor names
regressor_names <- first_betas_xii$meta$cifti$names
cat("  Regressor names:", paste(regressor_names, collapse=", "), "\n\n")

# Build beta matrices: one matrix per HRF, each with dimensions (resolution x n_regressors)
cat("Building beta matrices for all HRFs...\n")
beta_arrays <- array(NA, dim = c(resolution, n_regressors, n_hrfs))

for (hrf_idx in 1:n_hrfs) {
  # Get beta xifti for this HRF
  beta_xii <- glm_result[["GLMs"]][[hrf_idx]][["betas"]]

  # Convert to matrix (resolution x n_regressors)
  # Each column is one regressor
  beta_mat <- as.matrix(beta_xii)

  # Store in array
  beta_arrays[, , hrf_idx] <- beta_mat

  # Progress indicator
  if (hrf_idx %% 20 == 0) {
    cat("  Processed", hrf_idx, "/", n_hrfs, "HRFs\n")
  }
}
cat("  Complete!\n\n")

# Extract the best betas for each voxel based on regularized HRF parameters
cat("Extracting best betas for each voxel...\n")
best_betas <- matrix(NA, nrow = resolution, ncol = n_regressors)

for (voxel in 1:resolution) {
  # Get regularized params for this voxel
  a1 <- subject_params_a1[voxel]
  b1 <- subject_params_b1[voxel]
  c <- subject_params_c[voxel]

  # Find corresponding HRF index
  hrf_idx <- find_hrf_index(a1, b1, c, hrf_grid)

  # Extract betas from the appropriate HRF for all regressors
  best_betas[voxel, ] <- beta_arrays[voxel, , hrf_idx]

  # Progress indicator
  if (voxel %% 5000 == 0) {
    cat("  Processed", voxel, "/", resolution, "voxels\n")
  }
}
cat("  Complete!\n\n")

# Verify dimensions
cat("Results:\n")
cat("  best_betas dimensions:", dim(best_betas), "(voxels x regressors)\n")
cat("  Total betas extracted:", nrow(best_betas) * ncol(best_betas), "\n\n")

# Create xifti object with regularized betas
cat("Creating xifti object...\n")
template_xii <- glm_result[["bestmodel_xii"]]
template_dscalar <- ciftiTools::convert_xifti(template_xii, to = "dscalar")

# Create new xifti with best betas for all regressors
best_betas_xifti <- ciftiTools::newdata_xifti(
  xifti = template_dscalar,
  newdata = best_betas,
  newnames = regressor_names  # Use original regressor names
)

cat("  Complete!\n\n")

# Summary statistics for each regressor
cat("Summary statistics by regressor:\n")
cat(rep("-", 70), "\n", sep="")
for (i in 1:n_regressors) {
  cat(sprintf("%-15s: mean=%7.2f, sd=%7.2f, range=[%7.2f, %7.2f]\n",
              regressor_names[i],
              mean(best_betas[, i], na.rm=TRUE),
              sd(best_betas[, i], na.rm=TRUE),
              min(best_betas[, i], na.rm=TRUE),
              max(best_betas[, i], na.rm=TRUE)))
}
cat(rep("-", 70), "\n", sep="")

# Plot examples
plot(best_betas_xifti, idx = 1, title = "cue betas (best hrfs)", fname=here("dev/htmlplots/cue_betas_s2.png"),
     zlim=c(-1.85,1.85))
plot(best_betas_xifti, idx = 6, title = "tongue betas (best hrfs)", fname=here("dev/htmlplots/tongue_betas_s2.png"),
     zlim=c(-1.49, 1.49))

plot(best_betas_xifti, idx = 1, title = "cue betas (best hrfs)", fname=here("dev/fuck.png"))
plot(fit_allHRFs_result[["subject_results"]][[1]][["glm_result"]][["GLMs"]][[24]][["betas"]], idx = 1,
     fname=here("dev/htmlplots/cue_s2.png"), title = "cue betas (canonical)",
     zlim=c(-1.85,1.85))
plot(fit_allHRFs_result[["subject_results"]][[1]][["glm_result"]][["GLMs"]][[24]][["betas"]], idx = 6,
     fname=here("dev/htmlplots/tongue_s2.png"), title = "tongue betas (canonical)",
     zlim=c(-1.49, 1.49))

# Optionally save the result
# saveRDS(best_betas_xifti, here("dev", "fixtures", "regularized_betas_subject1.rds"))
# cat("\nSaved to:", here("dev", "fixtures", "regularized_betas_subject1.rds"), "\n")

cat("\nExtraction complete!\n")
cat("\nThe object 'best_betas_xifti' contains all", n_regressors, "regressors.\n")
cat("Each voxel's betas come from its regularized HRF model.\n")

###### Verification
