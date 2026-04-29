# Create test fixtures matching the 45 retest subjects
# This extracts only the test data for subjects that have retest data

library(dplyr)

# Helper to extract subject ID from file path
extract_subject_id <- function(file_path) {
  match <- regmatches(file_path, regexpr("/[0-9]{6}/", file_path))
  if (length(match) > 0) return(gsub("/", "", match))
  else return(NA)
}

# Output directory
output_dir <- "~/Documents/Github/hrf-z/dev/exploration/fixtures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Process GAMBLING LR and RL
for (run in c("LR", "RL")) {
  sess <- paste0("tfMRI_GAMBLING_", run)

  cat("\n========================================\n")
  cat("Processing", sess, "\n")
  cat("========================================\n")

  # Load full test session data
  test_file <- paste0("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures/session_data_", sess, ".rds")
  session_data_test <- readRDS(test_file)

  # Load retest session data
  retest_file <- paste0("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures/retest_session_data_", sess, ".rds")
  session_data_retest <- readRDS(retest_file)

  # Extract subject IDs
  subjects_test <- sapply(session_data_test$BOLD_files, extract_subject_id)
  subjects_retest <- sapply(session_data_retest$BOLD_files, extract_subject_id)

  cat("Full test data has", length(subjects_test), "subjects\n")
  cat("Retest data has", length(subjects_retest), "subjects\n")

  # Find matching subjects in test data
  match_indices <- which(subjects_test %in% subjects_retest)

  cat("Found", length(match_indices), "matching subjects in test data\n")

  # Verify all retest subjects are found
  stopifnot(all(subjects_retest %in% subjects_test))
  cat("✓ All retest subjects found in test data\n")

  # Create subset session data
  session_data_test_matched <- list(
    BOLD_files = session_data_test$BOLD_files[match_indices],
    EVs_list = session_data_test$EVs_list[match_indices],
    nuisance_files = session_data_test$nuisance_files[match_indices]
  )

  # Verify subject IDs match
  subjects_test_matched <- sapply(session_data_test_matched$BOLD_files, extract_subject_id)

  # Sort both to compare (they might be in different order)
  cat("\nTest subjects (first 10):", head(sort(subjects_test_matched), 10), "\n")
  cat("Retest subjects (first 10):", head(sort(subjects_retest), 10), "\n")

  # Check they have same subjects (regardless of order)
  if (setequal(subjects_test_matched, subjects_retest)) {
    cat("✓ Test and retest fixtures have matching subject sets\n")
  } else {
    cat("✗ WARNING: Subject sets do not match!\n")
    cat("Missing from test:", setdiff(subjects_retest, subjects_test_matched), "\n")
    cat("Missing from retest:", setdiff(subjects_test_matched, subjects_retest), "\n")
  }

  # Save matched test fixture
  output_file <- file.path(output_dir, paste0("test_matched_session_data_", sess, ".rds"))
  saveRDS(session_data_test_matched, output_file)
  cat("\n✓ Saved matched test fixture to:", output_file, "\n")
  cat("  Contains", length(session_data_test_matched$BOLD_files), "subjects\n")
}

cat("\n========================================\n")
cat("All matched test fixtures created!\n")
cat("Location:", output_dir, "\n")
cat("========================================\n")
