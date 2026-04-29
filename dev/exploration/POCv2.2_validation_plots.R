library(dplyr)
library(ggplot2)
library(tidyr)

extract_subject_id <- function(file_path) {
  match <- regmatches(file_path, regexpr("/[0-9]{6}/", file_path))
  if (length(match) > 0) return(gsub("/", "", match))
  else return(NA)
}

# Process both LR and RL
run_names <- c('LR')
# Change to c('LR', 'RL') to run both

# Local paths (after rsync from HPC)
results_dir <- "/Volumes/LaCie/root/hrf_adaptation/validation_bonferroni_0.05_smoothed_RSS/POCv2.2"
fixture_dir <- "~/Documents/Github/hrf-z/dev/exploration/fixtures"

for (run in run_names) {
  cat("\n========================================\n")
  cat("Processing", run, "run\n")
  cat("========================================\n\n")

  sess <- paste0("tfMRI_GAMBLING_", run)

  # Load session data to get actual subject IDs (45 matched subjects)
  session_data_test_matched <- readRDS(file.path(fixture_dir, paste0('test_matched_session_data_', sess, '.rds')))
  session_data_retest <- readRDS(file.path('~/Documents/Github/HRF-Adaptation-paper/Code/fixtures', paste0('retest_session_data_', sess, '.rds')))

  # Extract subject IDs from file paths (in order)
  subjects_test_matched <- sapply(session_data_test_matched[["BOLD_files"]], extract_subject_id)
  subjects_retest <- sapply(session_data_retest[["BOLD_files"]], extract_subject_id)

  cat("Test_matched", run, "has", length(subjects_test_matched), "subjects\n")
  cat("Retest", run, "has", length(subjects_retest), "subjects\n\n")

  # Load POCv2.2 results
  results_test_matched <- readRDS(file.path(results_dir, paste0("results_POCv2.2_test_matched_gambling_", run, "_weighted.rds")))
  results_retest <- readRDS(file.path(results_dir, paste0("results_POCv2.2_retest_gambling_", run, "_weighted.rds")))

  cat("Test_matched", run, "results has", length(unique(results_test_matched$subject_id)), "unique subject_ids\n")
  cat("Retest", run, "results has", length(unique(results_retest$subject_id)), "unique subject_ids\n\n")

  # Add actual subject IDs (matching by index position)
  results_test_matched$actual_subject_id <- subjects_test_matched[results_test_matched$subject_id]
  results_retest$actual_subject_id <- subjects_retest[results_retest$subject_id]

  # Test: verify no NAs in actual_subject_id
  stopifnot(sum(is.na(results_test_matched$actual_subject_id)) == 0)
  stopifnot(sum(is.na(results_retest$actual_subject_id)) == 0)
  cat("All subjects successfully mapped!\n\n")

  # Find best candidate for each subject using WEIGHTED RSS
  best_test_matched <- results_test_matched %>%
    group_by(actual_subject_id) %>%
    slice_min(weighted_RSS, n = 1) %>%
    ungroup()

  best_retest <- results_retest %>%
    group_by(actual_subject_id) %>%
    slice_min(weighted_RSS, n = 1) %>%
    ungroup()

  # Merge test and retest by actual subject ID
  best_test_compare <- best_test_matched %>%
    select(actual_subject_id,
           candidate_id_test = candidate_id,
           a1_offset_test = a1_offset,
           b1_offset_test = b1_offset,
           t2p_test = t2p,
           fwhm_test = fwhm)

  best_retest_compare <- best_retest %>%
    select(actual_subject_id,
           candidate_id_retest = candidate_id,
           a1_offset_retest = a1_offset,
           b1_offset_retest = b1_offset,
           t2p_retest = t2p,
           fwhm_retest = fwhm)

  comparison <- inner_join(best_test_compare, best_retest_compare, by = "actual_subject_id")

  cat("Number of subjects with both test and retest data:", nrow(comparison), "\n\n")

  # Combined best tables (for histogram/contingency)
  best_test_matched$run <- "Test"
  best_retest$run <- "Retest"
  best_combined <- bind_rows(best_test_matched, best_retest)

  # Contingency table with offsets by run
  contingency_with_offsets <- best_combined %>%
    count(candidate_id, a1_offset, b1_offset, run) %>%
    pivot_wider(names_from = run, values_from = n, values_fill = 0) %>%
    arrange(candidate_id)

  cat("Contingency (candidate_id x offsets x run):\n")
  print(knitr::kable(contingency_with_offsets))
  cat("\n")

  # Pearson correlations
  cor_candidate <- cor(comparison$candidate_id_test, comparison$candidate_id_retest)
  cor_a1 <- cor(comparison$a1_offset_test, comparison$a1_offset_retest)
  cor_b1 <- cor(comparison$b1_offset_test, comparison$b1_offset_retest)
  cor_t2p <- cor(comparison$t2p_test, comparison$t2p_retest)
  cor_fwhm <- cor(comparison$fwhm_test, comparison$fwhm_retest)

  cat("Candidate ID correlation:", cor_candidate, "\n")
  cat("a1_offset correlation:", cor_a1, "\n")
  cat("b1_offset correlation:", cor_b1, "\n")
  cat("t2p correlation:", cor_t2p, "\n")
  cat("fwhm correlation:", cor_fwhm, "\n\n")

  # Create output directory
  plot_dir <- "~/Documents/Github/hrf-z/dev/exploration/plots/POCv2.2_validation"
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  # Histogram (overlapping)
  p_hist <- ggplot(best_combined, aes(x = candidate_id, fill = run)) +
    geom_histogram(binwidth = 1, position = "identity", alpha = 0.6, color = "#e9ecef", linewidth = 0.5) +
    scale_fill_manual(values = c("Test" = "#404080", "Retest" = "#69b3a2")) +
    scale_x_continuous(breaks = 1:25) +
    labs(title = paste0("POCv2.2: Distribution of Best Candidate Maps (", run, ")"),
         x = "Candidate ID",
         y = "Number of Subjects",
         fill = "") +
    theme_minimal()

  ggsave(file.path(plot_dir, paste0("POCv2.2_candidate_hist_", run, ".png")), p_hist, width = 6, height = 4, dpi = 300)

  # Plot 1: Scatter plot - Candidate ID
  p_candidate <- ggplot(comparison, aes(x = candidate_id_test, y = candidate_id_retest)) +
    geom_point(alpha = 0.5, size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste0("POCv2.2 Test vs Retest: Candidate ID (", run, ")"),
         subtitle = paste0("r = ", round(cor_candidate, 3), " (n = ", nrow(comparison), ")"),
         x = "Test Candidate ID", y = "Retest Candidate ID") +
    theme_minimal()

  ggsave(file.path(plot_dir, paste0("POCv2.2_candidate_id_", run, ".png")), p_candidate, width = 6, height = 6, dpi = 300)

  # Plot 2: Scatter plot - a1 offset
  p_a1 <- ggplot(comparison, aes(x = a1_offset_test, y = a1_offset_retest)) +
    geom_jitter(alpha = 0.5, size = 3, width = 0.1, height = 0.1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste0("POCv2.2 Test vs Retest: a1 Offset (", run, ")"),
         subtitle = paste0("r = ", round(cor_a1, 3), " (n = ", nrow(comparison), ")"),
         x = "Test a1 Offset", y = "Retest a1 Offset") +
    theme_minimal()

  ggsave(file.path(plot_dir, paste0("POCv2.2_a1_offset_", run, ".png")), p_a1, width = 6, height = 6, dpi = 300)

  # Plot 3: Scatter plot - b1 offset
  p_b1 <- ggplot(comparison, aes(x = b1_offset_test, y = b1_offset_retest)) +
    geom_jitter(alpha = 0.5, size = 3, width = 0.02, height = 0.02) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste0("POCv2.2 Test vs Retest: b1 Offset (", run, ")"),
         subtitle = paste0("r = ", round(cor_b1, 3), " (n = ", nrow(comparison), ")"),
         x = "Test b1 Offset", y = "Retest b1 Offset") +
    theme_minimal()

  ggsave(file.path(plot_dir, paste0("POCv2.2_b1_offset_", run, ".png")), p_b1, width = 6, height = 6, dpi = 300)

  # Plot 4: Scatter plot - t2p (time to peak)
  p_t2p <- ggplot(comparison, aes(x = t2p_test, y = t2p_retest)) +
    geom_jitter(alpha = 0.5, size = 3, width = 0.05, height = 0.05) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste0("POCv2.2 Test vs Retest: Time to Peak (", run, ")"),
         subtitle = paste0("r = ", round(cor_t2p, 3), " (n = ", nrow(comparison), ")"),
         x = "Test t2p (sec)", y = "Retest t2p (sec)") +
    theme_minimal()

  ggsave(file.path(plot_dir, paste0("POCv2.2_t2p_", run, ".png")), p_t2p, width = 6, height = 6, dpi = 300)

  # Plot 5: Scatter plot - fwhm
  p_fwhm <- ggplot(comparison, aes(x = fwhm_test, y = fwhm_retest)) +
    geom_jitter(alpha = 0.5, size = 3, width = 0.05, height = 0.05) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste0("POCv2.2 Test vs Retest: FWHM (", run, ")"),
         subtitle = paste0("r = ", round(cor_fwhm, 3), " (n = ", nrow(comparison), ")"),
         x = "Test FWHM (sec)", y = "Retest FWHM (sec)") +
    theme_minimal()

  ggsave(file.path(plot_dir, paste0("POCv2.2_fwhm_", run, ".png")), p_fwhm, width = 6, height = 6, dpi = 300)

  cat("\nPlots saved to:", plot_dir, "\n")

  # Show the data
  cat("\nFirst 10 matched subjects:\n")
  print(head(comparison, 10))

  cat("\n")
}

cat("\n========================================\n")
cat("All POCv2.2 plots completed!\n")
cat("========================================\n")
