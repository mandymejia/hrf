library(dplyr)
library(ggplot2)
library(tidyr)

extract_subject_id <- function(file_path) {
  match <- regmatches(file_path, regexpr("/[0-9]{6}/", file_path))
  if (length(match) > 0) return(gsub("/", "", match))
  else return(NA)
}

# Load session data to get actual subject IDs
session_data_LR_original <- readRDS("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures/session_data_tfMRI_GAMBLING_LR.rds")
session_data_LR_retest <- readRDS("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures/retest_session_data_tfMRI_GAMBLING_LR.rds")

# Extract subject IDs from file paths (in order)
subjects_LR_original <- sapply(session_data_LR_original[["BOLD_files"]], extract_subject_id)
subjects_LR_retest <- sapply(session_data_LR_retest[["BOLD_files"]], extract_subject_id)

cat("Original LR has", length(subjects_LR_original), "subjects\n")
cat("Retest LR has", length(subjects_LR_retest), "subjects\n")
cat("First 10 retest subjects:", head(subjects_LR_retest, 10), "\n\n")

# Load results
lr_original <- readRDS("~/Downloads/cm_poc/results_cm_poc_gambling_LR.rds")
lr_retest <- readRDS("/Volumes/LaCie/root/hrf_adaptation/validation_bonferroni_0.05/cm_poc_retest/results_cm_poc_retest_gambling_LR.rds")

cat("Original LR results has", length(unique(lr_original$subject_id)), "unique subject_ids\n")
cat("Retest LR results has", length(unique(lr_retest$subject_id)), "unique subject_ids\n\n")

# Add actual subject IDs (matching by index position)
lr_original$actual_subject_id <- subjects_LR_original[lr_original$subject_id]
lr_retest$actual_subject_id <- subjects_LR_retest[lr_retest$subject_id]

# Test: verify no NAs in actual_subject_id
stopifnot(sum(is.na(lr_original$actual_subject_id)) == 0)
stopifnot(sum(is.na(lr_retest$actual_subject_id)) == 0)
cat("All subjects successfully mapped!\n\n")

# Find best candidate for each subject
best_lr_original <- lr_original %>%
  group_by(actual_subject_id) %>%
  slice_min(total_RSS, n = 1) %>%
  ungroup()

best_lr_retest <- lr_retest %>%
  group_by(actual_subject_id) %>%
  slice_min(total_RSS, n = 1) %>%
  ungroup()

# Merge test and retest by actual subject ID
best_lr_original_compare <- best_lr_original %>%
  select(actual_subject_id, candidate_id_test = candidate_id, a1_offset_test = a1_offset, b1_offset_test = b1_offset)

best_lr_retest_compare <- best_lr_retest %>%
  select(actual_subject_id, candidate_id_retest = candidate_id, a1_offset_retest = a1_offset, b1_offset_retest = b1_offset)

test_retest_comparison <- inner_join(best_lr_original_compare, best_lr_retest_compare, by = "actual_subject_id")

cat("Number of subjects with both test and retest data:", nrow(test_retest_comparison), "\n\n")

# Combined best tables (for histogram/contingency)
best_lr_original$run <- "Test"
best_lr_retest$run <- "Retest"
best_combined <- bind_rows(best_lr_original, best_lr_retest)

# Contingency table with offsets by run
contingency_with_offsets <- best_combined %>%
  count(candidate_id, a1_offset, b1_offset, run) %>%
  pivot_wider(names_from = run, values_from = n, values_fill = 0) %>%
  arrange(candidate_id)

cat("Contingency (candidate_id x offsets x run):\n")
print(knitr::kable(contingency_with_offsets))
cat("\n")

# Histogram of best candidate IDs by run (overlaid)
p_hist <- ggplot(best_combined, aes(x = candidate_id, fill = run)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.5) +
  scale_fill_manual(values = c("Test" = "blue", "Retest" = "red")) +
  labs(title = "Distribution of Best Candidate Maps: Test vs Retest (LR)",
       x = "Candidate ID (1-25)",
       y = "Number of Subjects") +
  theme_minimal()

# Pearson correlations
cor_candidate <- cor(test_retest_comparison$candidate_id_test, test_retest_comparison$candidate_id_retest)
cor_a1 <- cor(test_retest_comparison$a1_offset_test, test_retest_comparison$a1_offset_retest)
cor_b1 <- cor(test_retest_comparison$b1_offset_test, test_retest_comparison$b1_offset_retest)

cat("Candidate ID correlation:", cor_candidate, "\n")
cat("a1_offset correlation:", cor_a1, "\n")
cat("b1_offset correlation:", cor_b1, "\n\n")

# Create output directory
plot_dir <- "~/Documents/Github/hrf-z/dev/exploration/plots/poc_retest"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Save histogram
ggsave(file.path(plot_dir, "test_retest_candidate_hist_LR.png"), p_hist, width = 6, height = 4, dpi = 300)

# Scatter plot - Candidate ID
p_candidate <- ggplot(test_retest_comparison, aes(x = candidate_id_test, y = candidate_id_retest)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Test vs Retest: Best Candidate ID per Subject (LR)",
       subtitle = paste0("r = ", round(cor_candidate, 3), " (n = ", nrow(test_retest_comparison), ")"),
       x = "Test Candidate ID", y = "Retest Candidate ID") +
  theme_minimal()

ggsave(file.path(plot_dir, "test_retest_candidate_id_LR.png"), p_candidate, width = 6, height = 6, dpi = 300)

# Scatter plot - a1 offset
p_a1 <- ggplot(test_retest_comparison, aes(x = a1_offset_test, y = a1_offset_retest)) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.1, height = 0.1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Test vs Retest: Best a1 Offset per Subject (LR)",
       subtitle = paste0("r = ", round(cor_a1, 3), " (n = ", nrow(test_retest_comparison), ")"),
       x = "Test a1 Offset", y = "Retest a1 Offset") +
  theme_minimal()

ggsave(file.path(plot_dir, "test_retest_a1_offset_LR.png"), p_a1, width = 6, height = 6, dpi = 300)

# Scatter plot - b1 offset
p_b1 <- ggplot(test_retest_comparison, aes(x = b1_offset_test, y = b1_offset_retest)) +
  geom_jitter(alpha = 0.5, size = 3, width = 0.02, height = 0.02) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "Test vs Retest: Best b1 Offset per Subject (LR)",
       subtitle = paste0("r = ", round(cor_b1, 3), " (n = ", nrow(test_retest_comparison), ")"),
       x = "Test b1 Offset", y = "Retest b1 Offset") +
  theme_minimal()

ggsave(file.path(plot_dir, "test_retest_b1_offset_LR.png"), p_b1, width = 6, height = 6, dpi = 300)

cat("\nPlots saved to:", plot_dir, "\n")

# Show the data
cat("\nFirst 10 matched subjects:\n")
print(head(test_retest_comparison, 10))
