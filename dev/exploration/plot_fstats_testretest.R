library(here)
library(ciftiTools)
devtools::load_all("~/Documents/Github/hrf-z")

ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')

plot_dir <- "~/Documents/Github/hrf-z/dev/exploration/plots/fstats_testretest"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Load session data (test/retest)
session_data_test <- readRDS("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures/session_data_tfMRI_GAMBLING_LR.rds")
session_data_retest <- readRDS("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures/retest_session_data_tfMRI_GAMBLING_LR.rds")

# Function to remap paths from HPC to mounted directory
remap_path <- function(hpc_path) {
  # Handle both hrf_adaptation and hcp_dcwan paths
  hpc_path <- gsub("^/N/project/hrf_adaptation", "~/mnt/hrf_adaptation", hpc_path)
  hpc_path <- gsub("^/N/project/hcp_dcwan", "~/mnt/hcp_dcwan", hpc_path)
  return(hpc_path)
}

# Remap all file paths in session data
session_data_test$BOLD_files <- lapply(session_data_test$BOLD_files, remap_path)
session_data_test$nuisance_files <- lapply(session_data_test$nuisance_files, remap_path)
session_data_retest$BOLD_files <- lapply(session_data_retest$BOLD_files, remap_path)
session_data_retest$nuisance_files <- lapply(session_data_retest$nuisance_files, remap_path)

# Extract subject IDs to match test and retest
extract_subject_id <- function(file_path) {
  match <- regmatches(file_path, regexpr("/[0-9]{6}/", file_path))
  if (length(match) > 0) return(gsub("/", "", match))
  else return(NA)
}

subjects_test <- sapply(session_data_test$BOLD_files, extract_subject_id)
subjects_retest <- sapply(session_data_retest$BOLD_files, extract_subject_id)

cat("Test data has", length(subjects_test), "subjects\n")
cat("Retest data has", length(subjects_retest), "subjects\n")
cat("First 5 retest subjects:", head(subjects_retest, 5), "\n\n")

# Match subjects present in both test and retest
matching_subjects <- intersect(subjects_test, subjects_retest)
matching_subjects <- matching_subjects[!is.na(matching_subjects)]
n_subjects <- length(matching_subjects)

if (n_subjects == 0) {
  stop("No overlapping subjects between test and retest data.")
}

cat("Matched", n_subjects, "subjects between test and retest datasets\n\n")

# Helper to run canonical multiGLM and return F-stat xifti
run_glm_fstat <- function(session_data, idx) {
  BOLD <- ciftiTools::read_cifti(
    session_data[["BOLD_files"]][[idx]],
    brainstructures = c("left", "right"),
    resamp_res = 10000
  )

  design <- make_design(
    EVs = session_data[["EVs_list"]][[idx]],
    nTime = ncol(BOLD),
    TR = 0.72,
    onset = FALSE,  # FALSE for GAMBLING
    offset = FALSE
  )

  # Convert to 3D array for multiGLM (needs at least 2 models)
  design_array <- array(design$design, dim = c(dim(design$design), 2))
  design_array[, , 2] <- design_array[, , 2] + rnorm(length(design_array[, , 2]))

  glm_result <- multiGLM(
    BOLD = BOLD,
    brainstructures = c("left", "right"),
    design = design_array,
    design_canonical = design_array[, , 1],
    nuisance = as.matrix(read.table(session_data[["nuisance_files"]][[idx]], header = FALSE)),
    TR = 0.72,
    hpf = 0.01,
    resamp_res = NULL
  )

  return(glm_result$Fstat_xii)
}

F_test_all <- NULL
F_retest_all <- NULL
template_xii <- NULL

for (i in seq_along(matching_subjects)) {
  subj_id <- matching_subjects[i]
  retest_idx <- which(subjects_retest == subj_id)
  test_idx <- which(subjects_test == subj_id)

  if (length(test_idx) == 0 || length(retest_idx) == 0) {
    cat("WARNING: Subject", subj_id, "not present in both datasets. Skipping.\n")
    next
  }

  cat("\n=== Processing subject", subj_id, "(", i, "/", n_subjects, ") ===\n")
  cat("  Test index:", test_idx, "| Retest index:", retest_idx, "\n")

  fstat_test <- run_glm_fstat(session_data_test, test_idx)
  fstat_retest <- run_glm_fstat(session_data_retest, retest_idx)

  if (is.null(template_xii)) {
    template_xii <- fstat_test
    n_vox <- nrow(as.matrix(fstat_test))
    F_test_all <- matrix(NA_real_, nrow = n_vox, ncol = n_subjects)
    F_retest_all <- matrix(NA_real_, nrow = n_vox, ncol = n_subjects)
  }

  F_test_all[, i] <- as.vector(as.matrix(fstat_test))
  F_retest_all[, i] <- as.vector(as.matrix(fstat_retest))

  rm(fstat_test, fstat_retest)
  gc()
}

F_test_avg_mat <- matrix(rowMeans(F_test_all, na.rm = TRUE), ncol = 1)
F_retest_avg_mat <- matrix(rowMeans(F_retest_all, na.rm = TRUE), ncol = 1)

F_test_avg_xii <- ciftiTools::newdata_xifti(template_xii, F_test_avg_mat)
F_retest_avg_xii <- ciftiTools::newdata_xifti(template_xii, F_retest_avg_mat)

# Correlation between averaged maps
cor_avg <- cor(
  as.vector(F_test_avg_mat),
  as.vector(F_retest_avg_mat),
  use = "pairwise.complete.obs"
)

plot(
  F_test_avg_xii,
  title = paste0("GAMBLING LR - Test Avg (", n_subjects, " subjects) - F-statistics"),
  color_mode = "sequential",
  zlim = c(0, 12),
  fname = file.path(plot_dir, paste0("GAMBLING_LR_test_avg_", n_subjects, "_subjects_Fstat.png"))
)

plot(
  F_retest_avg_xii,
  title = paste0("GAMBLING LR - Retest Avg (", n_subjects, " subjects) - F-statistics"),
  color_mode = "sequential",
  zlim = c(0, 12),
  fname = file.path(plot_dir, paste0("GAMBLING_LR_retest_avg_", n_subjects, "_subjects_Fstat.png"))
)

cat("\nF-statistic correlation (test vs retest averages):", round(cor_avg, 3), "\n")
cat("All plots saved to:", plot_dir, "\n")
