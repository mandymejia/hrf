library(dplyr)
devtools::load_all("~/Documents/Github/hrf-z")

########################################
##### BUILD DENSE LOOKUP TABLE #########
########################################

build_dense_hrf_lookup <- function(c_val = 1/6,
                                    a1_range = c(4, 11),
                                    b1_range = c(0.5, 2.0),
                                    decimal_places = 2) {

  step <- 10^(-decimal_places)

  a1_vals <- seq(a1_range[1], a1_range[2], by = step)
  b1_vals <- seq(b1_range[1], b1_range[2], by = step)

  cat("Building dense lookup table:\n")
  cat("  a1:", length(a1_vals), "values from", a1_range[1], "to", a1_range[2], "step", step, "\n")
  cat("  b1:", length(b1_vals), "values from", b1_range[1], "to", b1_range[2], "step", step, "\n")
  cat("  c:", c_val, "\n")

  dense <- expand.grid(a1 = a1_vals, b1 = b1_vals)
  dense$c <- c_val
  cat("  Total combos:", nrow(dense), "\n\n")

  # Compute t2p (analytical)
  dense$t2p <- dense$a1 - dense$b1

  # Compute fwhm (numerical, with tapering like fit_allHRFs)
  cat("Computing FWHM for", nrow(dense), "grid points (with tapering)...\n")
  inds <- seq(0.01, 40, 0.01)

  dense$fwhm <- NA
  pb_interval <- max(1, floor(nrow(dense) / 20))

  for (i in 1:nrow(dense)) {
    if (i %% pb_interval == 0) cat("  Progress:", round(i/nrow(dense)*100), "%\n")

    a1 <- dense$a1[i]
    b1 <- dense$b1[i]
    a2 <- (16 / sqrt(6)) * sqrt(a1) * sqrt(b1)

    # Raw HRF first
    hrf_raw <- hrf::HRF_calc(t = inds, deriv = 0, a1 = a1, b1 = b1,
                              a2 = a2, b2 = b1, c = c_val)

    # Apply tapering if c > 0 (same logic as fit_allHRFs)
    if (c_val > 0) {
      peak2_idx <- which.min(hrf_raw)
      peak2_time <- inds[peak2_idx]
      taper_start <- min(peak2_time, 25)

      # Only taper if HRF hasn't resolved by 30s
      if (abs(hrf_raw[which.min(abs(inds - 30))]) > 0.01) {
        hrf <- hrf::HRF_calc(t = inds, deriv = 0, a1 = a1, b1 = b1,
                              a2 = a2, b2 = b1, c = c_val,
                              taper_start = taper_start, taper_end = 30, taper_power = 1)
      } else {
        hrf <- hrf_raw
      }
    } else {
      hrf <- hrf_raw
    }

    peak_val <- max(hrf)
    peak_idx <- which.max(hrf)
    vals_left <- hrf[1:peak_idx]
    vals_right <- hrf[peak_idx:length(hrf)]

    x1 <- min(which(vals_left > 0.5 * peak_val))
    x2 <- max(which(vals_right > 0.5 * peak_val))

    dense$fwhm[i] <- inds[peak_idx + x2 - 1] - inds[x1]
  }

  cat("Done!\n\n")

  # Precompute normalized ranges for distance calculation
  t2p_range <- max(dense$t2p) - min(dense$t2p)
  fwhm_range <- max(dense$fwhm) - min(dense$fwhm)

  cat("t2p range:", range(dense$t2p), "\n")
  cat("fwhm range:", range(dense$fwhm), "\n")

  return(list(
    grid = dense,
    c_val = c_val,
    t2p_range = t2p_range,
    fwhm_range = fwhm_range
  ))
}

# Inverse lookup: given (t2p, fwhm) -> find (a1, b1)
inverse_lookup <- function(target_t2p, target_fwhm, lookup) {
  dense <- lookup$grid

  # Normalized distance
  dist <- sqrt(((dense$t2p - target_t2p) / lookup$t2p_range)^2 +
               ((dense$fwhm - target_fwhm) / lookup$fwhm_range)^2)

  best_idx <- which.min(dist)

  c(a1 = dense$a1[best_idx], b1 = dense$b1[best_idx])
}

# Vectorized inverse lookup for many voxels
inverse_lookup_batch <- function(target_t2p_vec, target_fwhm_vec, lookup) {
  n <- length(target_t2p_vec)
  results <- matrix(NA, nrow = n, ncol = 2)
  colnames(results) <- c("a1", "b1")

  for (i in 1:n) {
    results[i,] <- inverse_lookup(target_t2p_vec[i], target_fwhm_vec[i], lookup)
  }

  return(as.data.frame(results))
}

########################################
##### BUILD AND SAVE LOOKUP ############
########################################

fixture_dir <- "~/Documents/Github/hrf-z/dev/exploration/fixtures"

cat("========================================\n")
cat("Building 1-decimal dense lookup table\n")
cat("========================================\n\n")

lookup_1dec <- build_dense_hrf_lookup(c_val = 1/6, decimal_places = 1)

# Save to fixtures
lookup_file <- file.path(fixture_dir, "dense_hrf_lookup_c16_1dec.rds")
saveRDS(lookup_1dec, lookup_file)
cat("\nSaved lookup table to:", lookup_file, "\n")
cat("File size:", round(file.size(lookup_file) / 1024 / 1024, 2), "MB\n")

########################################
##### TEST IT ##########################
########################################

cat("\n========================================\n")
cat("Testing dense lookup table\n")
cat("========================================\n\n")

# Speed test with realistic inputs from the valid region
set.seed(42)
n_test <- 5000

# Generate test inputs from VALID region (use actual grid to get valid t2p/fwhm pairs)
valid_indices <- sample(1:nrow(lookup_1dec$grid), n_test, replace = TRUE)
test_t2p <- lookup_1dec$grid$t2p[valid_indices] + runif(n_test, -0.1, 0.1)  # small perturbation
test_fwhm <- lookup_1dec$grid$fwhm[valid_indices] + runif(n_test, -0.1, 0.1)

cat("--- Speed test:", n_test, "voxels ---\n")
start_time <- Sys.time()
results <- inverse_lookup_batch(test_t2p, test_fwhm, lookup_1dec)
end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("Time:", round(elapsed, 3), "seconds\n")
cat("Per voxel:", round(elapsed/n_test * 1000, 4), "ms\n")

# Verify results
results$t2p_in <- test_t2p
results$fwhm_in <- test_fwhm
results$t2p_reconstructed <- results$a1 - results$b1
results$t2p_error <- abs(results$t2p_reconstructed - results$t2p_in)

cat("\nSample results:\n")
print(head(results, 10))

cat("\nt2p reconstruction error:\n")
cat("  Max:", round(max(results$t2p_error), 4), "\n")
cat("  Mean:", round(mean(results$t2p_error), 4), "\n")

cat("\n========================================\n")
cat("Step 0 complete! Lookup table ready.\n")
cat("========================================\n")
