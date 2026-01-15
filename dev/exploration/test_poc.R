library(dplyr)

extract_subject_id <- function(file_path) {
  match <- regmatches(file_path, regexpr("/[0-9]{6}/", file_path))
  if (length(match) > 0) return(gsub("/", "", match))
  else return(NA)
}

cat("*** LOADING SESSION DATA ***\n")
session_data_LR <- readRDS("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures/session_data_tfMRI_GAMBLING_LR.rds")
session_data_RL <- readRDS("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures/session_data_tfMRI_GAMBLING_RL.rds")

# Extract subject IDs from file paths (in order)
subjects_LR <- sapply(session_data_LR[["BOLD_files"]], extract_subject_id)
subjects_RL <- sapply(session_data_RL[["BOLD_files"]], extract_subject_id)

cat("LR has", length(subjects_LR), "subjects\n")
cat("RL has", length(subjects_RL), "subjects\n\n")

cat("*** FIRST 10 SUBJECTS IN EACH ***\n")
cat("LR:", head(subjects_LR, 10), "\n")
cat("RL:", head(subjects_RL, 10), "\n\n")

cat("*** CHECKING OVERLAP ***\n")
subjects_both <- intersect(subjects_LR, subjects_RL)
subjects_LR_only <- setdiff(subjects_LR, subjects_RL)
subjects_RL_only <- setdiff(subjects_RL, subjects_LR)

cat("Subjects in both LR and RL:", length(subjects_both), "\n")
cat("Subjects in LR only:", length(subjects_LR_only), "\n")
cat("Subjects in RL only:", length(subjects_RL_only), "\n\n")

if (length(subjects_LR_only) > 0) {
  cat("LR-only subjects:", head(subjects_LR_only, 10), "...\n")
}
if (length(subjects_RL_only) > 0) {
  cat("RL-only subjects:", head(subjects_RL_only, 10), "...\n")
}
cat("\n")

cat("*** CHECKING ORDER DIFFERENCES ***\n")
# Find a few subjects that exist in both and check their positions
test_subjects <- head(subjects_both, 5)
for (subj in test_subjects) {
  pos_LR <- which(subjects_LR == subj)
  pos_RL <- which(subjects_RL == subj)
  cat(sprintf("Subject %s: position %d in LR, position %d in RL\n", subj, pos_LR, pos_RL))
}
cat("\n")

cat("*** LOADING POC RESULTS ***\n")
lr <- readRDS("~/Downloads/cm_poc/results_cm_poc_gambling_LR.rds")
rl <- readRDS("~/Downloads/cm_poc/results_cm_poc_gambling_RL.rds")

cat("LR results has", nrow(lr), "rows,", length(unique(lr$subject_id)), "unique subject_ids\n")
cat("RL results has", nrow(rl), "rows,", length(unique(rl$subject_id)), "unique subject_ids\n")
cat("subject_id ranges: LR [", min(lr$subject_id), "-", max(lr$subject_id), "], RL [", min(rl$subject_id), "-", max(rl$subject_id), "]\n\n")

cat("*** MAPPING subject_id TO actual_subject_id ***\n")
lr$actual_subject_id <- subjects_LR[lr$subject_id]
rl$actual_subject_id <- subjects_RL[rl$subject_id]

cat("LR: NAs after mapping =", sum(is.na(lr$actual_subject_id)), "\n")
cat("RL: NAs after mapping =", sum(is.na(rl$actual_subject_id)), "\n\n")

cat("*** SHOWING EXAMPLE MAPPINGS ***\n")
# Show first few rows of each to verify mapping
cat("LR sample:\n")
print(lr %>% select(subject_id, actual_subject_id, candidate_id, total_RSS) %>% head(10))
cat("\nRL sample:\n")
print(rl %>% select(subject_id, actual_subject_id, candidate_id, total_RSS) %>% head(10))
cat("\n")

cat("*** TESTING SUBJECT MATCHING ***\n")
# Find best candidate for each subject
best_lr <- lr %>%
  group_by(actual_subject_id) %>%
  slice_min(total_RSS, n = 1) %>%
  ungroup() %>%
  select(actual_subject_id, candidate_id_lr = candidate_id, a1_offset_lr = a1_offset, b1_offset_lr = b1_offset)

best_rl <- rl %>%
  group_by(actual_subject_id) %>%
  slice_min(total_RSS, n = 1) %>%
  ungroup() %>%
  select(actual_subject_id, candidate_id_rl = candidate_id, a1_offset_rl = a1_offset, b1_offset_rl = b1_offset)

cat("best_lr has", nrow(best_lr), "subjects\n")
cat("best_rl has", nrow(best_rl), "subjects\n\n")

# Join to see matched subjects
lr_rl_matched <- inner_join(best_lr, best_rl, by = "actual_subject_id")
cat("After inner_join:", nrow(lr_rl_matched), "subjects with both LR and RL data\n")
cat("This should match subjects_both:", length(subjects_both), "\n")
stopifnot(nrow(lr_rl_matched) == length(subjects_both))
cat("Match confirmed!\n\n")

cat("*** EXAMPLE MATCHED SUBJECTS ***\n")
print(head(lr_rl_matched, 10))
cat("\n")

cat("*** VERIFYING SPECIFIC SUBJECT ***\n")
# Pick a subject that exists in both and verify end-to-end
test_subj <- subjects_both[1]
cat("Testing subject:", test_subj, "\n")

# Find positions in session_data
pos_LR <- which(subjects_LR == test_subj)
pos_RL <- which(subjects_RL == test_subj)
cat(sprintf("Position in session_data: LR=%d, RL=%d\n", pos_LR, pos_RL))

# Find in results (should match positions)
lr_rows <- lr %>% filter(subject_id == pos_LR)
rl_rows <- rl %>% filter(subject_id == pos_RL)
cat(sprintf("Rows in results: LR=%d (subject_id=%d), RL=%d (subject_id=%d)\n",
            nrow(lr_rows), unique(lr_rows$subject_id),
            nrow(rl_rows), unique(rl_rows$subject_id)))

# Verify actual_subject_id matches
cat(sprintf("actual_subject_id in results: LR=%s, RL=%s\n",
            unique(lr_rows$actual_subject_id),
            unique(rl_rows$actual_subject_id)))
stopifnot(unique(lr_rows$actual_subject_id) == test_subj)
stopifnot(unique(rl_rows$actual_subject_id) == test_subj)
cat("Mapping verified for", test_subj, "\n\n")

cat("*** ALL TESTS PASSED! ***\n")
