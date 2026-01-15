library(dplyr)
library(ggplot2)
library(tidyr)

extract_subject_id <- function(file_path) {
  match <- regmatches(file_path, regexpr("/[0-9]{6}/", file_path))
  if (length(match) > 0) return(gsub("/", "", match))
  else return(NA)
}

# Load session data to get actual subject IDs
session_data_LR <- readRDS("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures/session_data_tfMRI_GAMBLING_LR.rds")
session_data_RL <- readRDS("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures/session_data_tfMRI_GAMBLING_RL.rds")

# Extract subject IDs from file paths (in order)
subjects_LR <- sapply(session_data_LR[["BOLD_files"]], extract_subject_id)
subjects_RL <- sapply(session_data_RL[["BOLD_files"]], extract_subject_id)

# Load results
lr <- readRDS("~/Downloads/cm_poc/results_cm_poc_gambling_LR.rds")
rl <- readRDS("~/Downloads/cm_poc/results_cm_poc_gambling_RL.rds")

# Test: verify number of subjects matches
cat("LR: session_data has", length(subjects_LR), "subjects, results has", length(unique(lr$subject_id)), "subjects\n")
cat("RL: session_data has", length(subjects_RL), "subjects, results has", length(unique(rl$subject_id)), "subjects\n")
stopifnot(length(subjects_LR) == length(unique(lr$subject_id)))
stopifnot(length(subjects_RL) == length(unique(rl$subject_id)))
cat("Subject counts match!\n\n")

# Add actual subject IDs (matching by index position)
lr$actual_subject_id <- subjects_LR[lr$subject_id]
rl$actual_subject_id <- subjects_RL[rl$subject_id]

# Test: verify no NAs in actual_subject_id
stopifnot(sum(is.na(lr$actual_subject_id)) == 0)
stopifnot(sum(is.na(rl$actual_subject_id)) == 0)
cat("All subjects successfully mapped!\n\n")

# Find best candidate for each subject
best_lr <- lr %>%
  group_by(actual_subject_id) %>%
  slice_min(total_RSS, n = 1) %>%
  ungroup()

best_rl <- rl %>%
  group_by(actual_subject_id) %>%
  slice_min(total_RSS, n = 1) %>%
  ungroup()

best_lr$run <- "LR"
best_rl$run <- "RL"
best_combined <- rbind(best_lr, best_rl)

# Histogram (side-by-side)
ggplot(best_combined, aes(x = candidate_id, fill = run)) +
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
  scale_fill_manual(values = c("LR" = "blue", "RL" = "red")) +
  labs(title = "Distribution of Best Candidate Maps: LR vs RL",
       x = "Candidate ID (1-25)",
       y = "Number of Subjects") +
  theme_minimal()

################################################################################

# Contingency table with offsets
contingency_with_offsets <- best_combined %>%
  count(candidate_id, a1_offset, b1_offset, run) %>%
  pivot_wider(names_from = run, values_from = n, values_fill = 0) %>%
  arrange(candidate_id)

knitr::kable(contingency_with_offsets)

################################################################################

# Merge LR and RL by actual subject ID
best_lr_compare <- best_lr %>%
  select(actual_subject_id, candidate_id_lr = candidate_id, a1_offset_lr = a1_offset, b1_offset_lr = b1_offset)

best_rl_compare <- best_rl %>%
  select(actual_subject_id, candidate_id_rl = candidate_id, a1_offset_rl = a1_offset, b1_offset_rl = b1_offset)

lr_rl_comparison <- inner_join(best_lr_compare, best_rl_compare, by = "actual_subject_id")

cat("Number of subjects with both LR and RL data:", nrow(lr_rl_comparison), "\n\n")

# Pearson correlations
cor_candidate <- cor(lr_rl_comparison$candidate_id_lr, lr_rl_comparison$candidate_id_rl)
cor_a1 <- cor(lr_rl_comparison$a1_offset_lr, lr_rl_comparison$a1_offset_rl)
cor_b1 <- cor(lr_rl_comparison$b1_offset_lr, lr_rl_comparison$b1_offset_rl)

cat("Candidate ID correlation:", cor_candidate, "\n")
cat("a1_offset correlation:", cor_a1, "\n")
cat("b1_offset correlation:", cor_b1, "\n\n")

# Scatter plots
ggplot(lr_rl_comparison, aes(x = candidate_id_lr, y = candidate_id_rl)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "LR vs RL: Best Candidate ID per Subject",
       subtitle = paste0("r = ", round(cor_candidate, 3), " (n = ", nrow(lr_rl_comparison), ")"),
       x = "LR Candidate ID", y = "RL Candidate ID") +
  theme_minimal()

ggplot(lr_rl_comparison, aes(x = a1_offset_lr, y = a1_offset_rl)) +
  geom_jitter(alpha = 0.5, width = 0.1, height = 0.1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "LR vs RL: Best a1 Offset per Subject",
       subtitle = paste0("r = ", round(cor_a1, 3), " (n = ", nrow(lr_rl_comparison), ")"),
       x = "LR a1 Offset", y = "RL a1 Offset") +
  theme_minimal()

ggplot(lr_rl_comparison, aes(x = b1_offset_lr, y = b1_offset_rl)) +
  geom_jitter(alpha = 0.5, width = 0.02, height = 0.02) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(title = "LR vs RL: Best b1 Offset per Subject",
       subtitle = paste0("r = ", round(cor_b1, 3), " (n = ", nrow(lr_rl_comparison), ")"),
       x = "LR b1 Offset", y = "RL b1 Offset") +
  theme_minimal()
