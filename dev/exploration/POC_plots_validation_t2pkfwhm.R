library(dplyr)
library(ggplot2)
library(tidyr)
devtools::load_all("~/Documents/Github/hrf-z")

extract_subject_id <- function(file_path) {
  match <- regmatches(file_path, regexpr("/[0-9]{6}/", file_path))
  if (length(match) > 0) return(gsub("/", "", match))
  else return(NA)
}

get_hrf_metrics <- function(a1, b1, c, tapered = TRUE) {
  n <- length(a1)
  time_to_peak <- numeric(n)
  FWHM <- numeric(n)

  inds <- seq(1/100, 20, 1/100) * 2

  for (i in 1:n) {
    a1_i <- a1[i]
    b1_i <- b1[i]
    c_i <- c[i]
    a2_i <- (16 / sqrt(6)) * sqrt(a1_i) * sqrt(b1_i)
    b2_i <- b1_i

    time_to_peak[i] <- a1_i - b1_i

    hrf_raw <- hrf::HRF_calc(t = inds, deriv = 0, a1 = a1_i, b1 = b1_i,
                              a2 = a2_i, b2 = b2_i, c = c_i)

    if (tapered && c_i > 0) {
      peak2_idx <- which.min(hrf_raw)
      peak2_time <- inds[peak2_idx]
      taper_start <- min(peak2_time, 25)

      if (abs(hrf_raw[which.min(abs(inds - 30))]) > 0.01) {
        hrf_use <- hrf::HRF_calc(t = inds, deriv = 0, a1 = a1_i, b1 = b1_i,
                                  a2 = a2_i, b2 = b2_i, c = c_i,
                                  taper_start = taper_start, taper_end = 30, taper_power = 1)
      } else {
        hrf_use <- hrf_raw
      }
    } else {
      hrf_use <- hrf_raw
    }

    peak_val <- max(hrf_use)
    peak_idx <- which.max(hrf_use)
    vals_left <- hrf_use[1:peak_idx]
    vals_right <- hrf_use[peak_idx:length(hrf_use)]

    x1 <- min(which(vals_left > 0.5 * peak_val))
    x2 <- max(which(vals_right > 0.5 * peak_val))
    time1 <- inds[x1]
    time2 <- inds[peak_idx + x2 - 1]
    FWHM[i] <- time2 - time1
  }

  return(data.frame(a1 = a1, b1 = b1, c = c, t2p = time_to_peak, fwhm = FWHM))
}

# Load regularize result
cat("Loading regularize result...\n")
reg_result <- readRDS("/Volumes/LaCie/root/hrf_adaptation/validation_bonferroni_0.05_smoothed/regularize_allHRFs/regularize_allHRFs_result_tfMRI_GAMBLING_LR_1084s.rds")
hrf_grid <- reg_result$hrf_grid

# snap_to_grid function from POC_validation.R
snap_to_grid <- function(param_df, hrf_grid, a1_col = "a1_pop", b1_col = "b1_pop", c_col = "c_pop") {
  grid_combos <- hrf_grid %>% select(a1, b1, c) %>% distinct()
  a_range <- max(grid_combos$a1) - min(grid_combos$a1)
  b_range <- max(grid_combos$b1) - min(grid_combos$b1)

  snapped_df <- param_df %>%
    rowwise() %>%
    mutate(
      grid_idx = {
        a1_val <- .data[[a1_col]]
        b1_val <- .data[[b1_col]]
        c_val <- .data[[c_col]]
        c_match <- abs(grid_combos$c - c_val) < 1e-6
        if (sum(c_match) == 0) {
          c_dists <- abs(grid_combos$c - c_val)
          c_match <- c_dists == min(c_dists)
        }
        candidates <- grid_combos[c_match, ]
        distances <- sqrt(((a1_val - candidates$a1) / a_range)^2 +
                          ((b1_val - candidates$b1) / b_range)^2)
        which(c_match)[which.min(distances)]
      }
    ) %>%
    ungroup() %>%
    mutate(
      a1_snapped = grid_combos$a1[grid_idx],
      b1_snapped = grid_combos$b1[grid_idx],
      c_snapped = grid_combos$c[grid_idx]
    ) %>%
    select(-grid_idx)
  return(snapped_df)
}

# Extract continuous population averages
pop_params <- reg_result$regularized_params
a1_pop <- pop_params$a1$results_OLS %>% select(voxel, a1_pop = param_mean0) %>% distinct()
b1_pop <- pop_params$b1$results_OLS %>% select(voxel, b1_pop = param_mean0) %>% distinct()
c_pop <- pop_params$c$results_OLS %>% select(voxel, c_pop = param_mean0) %>% distinct()

group_map <- a1_pop %>%
  left_join(b1_pop, by = "voxel") %>%
  left_join(c_pop, by = "voxel")

# Snap to grid
cat("Snapping population averages to grid...\n")
group_map_snapped <- snap_to_grid(group_map, hrf_grid)

voxel_base <- group_map_snapped %>%
  select(voxel, a1_base = a1_snapped, b1_base = b1_snapped, c_base = c_snapped)

cat("Loaded base params for", nrow(voxel_base), "voxels\n")

# Precompute t2p/fwhm for all possible (base + offset) combinations
a1_offsets <- c(-2, -1, 0, 1, 2)
b1_offsets <- c(-0.5, -0.25, 0, 0.25, 0.5)

unique_base <- voxel_base %>%
  select(a1_base, b1_base, c_base) %>%
  distinct()

cat("Found", nrow(unique_base), "unique snapped (a1, b1, c) tuples\n")

# Apply all offset combinations to each unique base tuple
offset_grid <- expand.grid(a1_offset = a1_offsets, b1_offset = b1_offsets)

all_combos <- do.call(rbind, lapply(1:nrow(offset_grid), function(i) {
  unique_base %>%
    mutate(
      actual_a1 = a1_base + offset_grid$a1_offset[i],
      actual_b1 = b1_base + offset_grid$b1_offset[i],
      actual_c = c_base
    ) %>%
    select(actual_a1, actual_b1, actual_c)
})) %>%
  distinct()

cat("Precomputing metrics for", nrow(all_combos), "unique (a1, b1, c) combinations...\n")
precomputed <- get_hrf_metrics(all_combos$actual_a1, all_combos$actual_b1, all_combos$actual_c, tapered = TRUE)
all_combos$t2p <- precomputed$t2p
all_combos$fwhm <- precomputed$fwhm
cat("Done.\n\n")

run_names <- c('LR')

for (run in run_names) {
  cat("\n========================================\n")
  cat("Processing", run, "run\n")
  cat("========================================\n\n")

  # Load session data to get actual subject IDs
  session_data_test_matched <- readRDS(paste0("~/Documents/Github/hrf-z/dev/exploration/fixtures/test_matched_session_data_tfMRI_GAMBLING_", run, ".rds"))
  session_data_retest <- readRDS(paste0("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures/retest_session_data_tfMRI_GAMBLING_", run, ".rds"))

  subjects_test_matched <- sapply(session_data_test_matched[["BOLD_files"]], extract_subject_id)
  subjects_retest <- sapply(session_data_retest[["BOLD_files"]], extract_subject_id)

  # Load POC results
  results_test_matched <- readRDS(paste0("~/Documents/Github/hrf-z/dev/exploration/poc_results/weighted_smoothed/results_cm_poc_test_matched_gambling_", run, "_weighted.rds"))
  results_retest <- readRDS(paste0("~/Documents/Github/hrf-z/dev/exploration/poc_results/weighted_smoothed/results_cm_poc_retest_gambling_", run, "_weighted.rds"))

  # Add actual subject IDs
  results_test_matched$actual_subject_id <- subjects_test_matched[results_test_matched$subject_id]
  results_retest$actual_subject_id <- subjects_retest[results_retest$subject_id]

  # Find best candidate for each subject
  best_test_matched <- results_test_matched %>%
    group_by(actual_subject_id) %>%
    slice_min(weighted_RSS, n = 1) %>%
    ungroup()

  best_retest <- results_retest %>%
    group_by(actual_subject_id) %>%
    slice_min(weighted_RSS, n = 1) %>%
    ungroup()

  # For each subject, compute average t2p and fwhm across all voxels
  # by applying their best offset to the per-voxel base values (using precomputed lookup)

  compute_subject_metrics <- function(a1_offset, b1_offset, voxel_base, lookup) {
    # Apply offset to all voxels
    voxel_params <- voxel_base %>%
      mutate(
        actual_a1 = a1_base + a1_offset,
        actual_b1 = b1_base + b1_offset,
        actual_c = c_base
      ) %>%
      left_join(lookup, by = c("actual_a1", "actual_b1", "actual_c"))

    # Return average across voxels
    return(c(t2p = mean(voxel_params$t2p, na.rm = TRUE),
             fwhm = mean(voxel_params$fwhm, na.rm = TRUE)))
  }

  cat("Computing per-subject metrics (test)...\n")
  test_metrics <- t(sapply(1:nrow(best_test_matched), function(i) {
    compute_subject_metrics(best_test_matched$a1_offset[i],
                            best_test_matched$b1_offset[i],
                            voxel_base, all_combos)
  }))
  best_test_matched$t2p <- test_metrics[, "t2p"]
  best_test_matched$fwhm <- test_metrics[, "fwhm"]

  cat("Computing per-subject metrics (retest)...\n")
  retest_metrics <- t(sapply(1:nrow(best_retest), function(i) {
    compute_subject_metrics(best_retest$a1_offset[i],
                            best_retest$b1_offset[i],
                            voxel_base, all_combos)
  }))
  best_retest$t2p <- retest_metrics[, "t2p"]
  best_retest$fwhm <- retest_metrics[, "fwhm"]

  # Merge for comparison
  best_test_compare <- best_test_matched %>%
    select(actual_subject_id, t2p_test = t2p, fwhm_test = fwhm)

  best_retest_compare <- best_retest %>%
    select(actual_subject_id, t2p_retest = t2p, fwhm_retest = fwhm)

  comparison <- inner_join(best_test_compare, best_retest_compare, by = "actual_subject_id")

  cat("Matched subjects:", nrow(comparison), "\n\n")

  # Correlations
  cor_t2p <- cor(comparison$t2p_test, comparison$t2p_retest, use = "complete.obs")
  cor_fwhm <- cor(comparison$fwhm_test, comparison$fwhm_retest, use = "complete.obs")

  cat("t2p correlation:", cor_t2p, "\n")
  cat("fwhm correlation:", cor_fwhm, "\n\n")

  # Create output directory
  plot_dir <- "~/Documents/Github/hrf-z/dev/exploration/plots/poc_validation_t2pkfwhm"
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  # Scatter plot - t2p
  p_t2p <- ggplot(comparison, aes(x = t2p_test, y = t2p_retest)) +
    geom_jitter(alpha = 0.5, size = 3, width = 0.1, height = 0.1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste0("Test vs Retest: Time to Peak per Subject (", run, ")"),
         subtitle = paste0("r = ", round(cor_t2p, 3), " (n = ", nrow(comparison), ")"),
         x = "Test t2p (sec)", y = "Retest t2p (sec)") +
    theme_minimal()

  ggsave(file.path(plot_dir, paste0("test_retest_t2p_", run, ".png")), p_t2p, width = 6, height = 6, dpi = 300)

  # Scatter plot - fwhm
  p_fwhm <- ggplot(comparison, aes(x = fwhm_test, y = fwhm_retest)) +
    geom_jitter(alpha = 0.5, size = 3, width = 0.05, height = 0.05) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    labs(title = paste0("Test vs Retest: FWHM per Subject (", run, ")"),
         subtitle = paste0("r = ", round(cor_fwhm, 3), " (n = ", nrow(comparison), ")"),
         x = "Test FWHM (sec)", y = "Retest FWHM (sec)") +
    theme_minimal()

  ggsave(file.path(plot_dir, paste0("test_retest_fwhm_", run, ".png")), p_fwhm, width = 6, height = 6, dpi = 300)

  cat("Plots saved to:", plot_dir, "\n")
}

cat("\n========================================\n")
cat("All plots completed!\n")
cat("========================================\n")
