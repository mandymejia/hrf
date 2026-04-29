library(dplyr)
devtools::load_all("~/Documents/Github/hrf-z")

fit_results_root_dir <- "/N/project/hrf_adaptation/validation_bonferroni_0.05_smoothed_RSS"

get_hrf_metrics <- function(a1, b1, c, tapered = TRUE) {
  n <- length(a1)
  time_to_peak <- numeric(n)
  FWHM <- numeric(n)

  inds <- seq(1/100, 20, 1/100) * 2  # TR=2, sr_factor=100, max_time=20

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

  return(data.frame(a1 = a1, b1 = b1, c = c, time_to_peak = time_to_peak, FWHM = FWHM))
}

# Reverse: snap t2p/fwhm back to nearest (a1, b1) grid point with fixed c
snap_to_grid_t2p_fwhm <- function(t2p, fwhm, hrf_grid, fixed_c) {
  # Build lookup table with t2p/fwhm for grid points with matching c
  grid_subset <- hrf_grid %>% filter(abs(c - fixed_c) < 1e-6)
  metrics <- get_hrf_metrics(grid_subset$a1, grid_subset$b1, grid_subset$c, tapered = TRUE)

  lookup <- data.frame(
    a1 = grid_subset$a1,
    b1 = grid_subset$b1,
    c = grid_subset$c,
    t2p = metrics$time_to_peak,
    fwhm = metrics$FWHM
  )

  # Normalize for distance calculation
  t2p_range <- max(lookup$t2p) - min(lookup$t2p)
  fwhm_range <- max(lookup$fwhm) - min(lookup$fwhm)

  # For each input (t2p, fwhm), find nearest grid point
  n <- length(t2p)
  a1_out <- numeric(n)
  b1_out <- numeric(n)

  for (i in 1:n) {
    dist <- sqrt(((t2p[i] - lookup$t2p) / t2p_range)^2 +
                 ((fwhm[i] - lookup$fwhm) / fwhm_range)^2)
    nearest_idx <- which.min(dist)
    a1_out[i] <- lookup$a1[nearest_idx]
    b1_out[i] <- lookup$b1[nearest_idx]
  }

  return(data.frame(a1 = a1_out, b1 = b1_out, c = fixed_c))
}

create_population_avg <- function(root_dir, sess_name) {

  # Step 1: Load results
  cat("Loading fit_workingHRF result...\n")
  fit_workingHRF_result <- readRDS(file.path(root_dir, "fit_workingHRF", paste0("fit_workingHRF_result_tfMRI_", sess_name, ".rds")))

  cat("Loading fit_allHRFs result...\n")
  fit_allHRFs_result <- readRDS(file.path(root_dir, "fit_allHRFs", paste0("fit_allHRFs_result_tfMRI_", sess_name, ".rds")))

  # Step 2: Extract masks
  # Population mask: which voxels have enough activation across subjects
  mask_prop_NA <- fit_workingHRF_result[["activation_masks"]][["mask_prop_NA"]]
  pop_mask_voxels <- which(!is.na(mask_prop_NA))
  cat("Population mask: ", length(pop_mask_voxels), " voxels pass threshold\n")

  # Subject masks: (voxels x subjects) matrix of TRUE/FALSE
  subject_masks <- fit_workingHRF_result[["activation_masks"]][["masks"]]
  cat("Subject masks: ", nrow(subject_masks), " voxels x ", ncol(subject_masks), " subjects\n")

  # Step 3: Get HRF grid and pre-compute t2p/fwhm for all grid points
  hrf_grid <- fit_allHRFs_result[["hrf_grid"]]
  c0_indices <- which(hrf_grid$c == 0)
  c1_indices <- which(hrf_grid$c != 0)
  cat("HRF grid: ", nrow(hrf_grid), " models (", length(c0_indices), " c=0, ", length(c1_indices), " c=1/6)\n")

  # Pre-compute t2p/fwhm for all grid points (do this ONCE)
  cat("Pre-computing t2p/fwhm for all grid points...\n")
  grid_metrics <- get_hrf_metrics(hrf_grid$a1, hrf_grid$b1, hrf_grid$c, tapered = TRUE)
  hrf_grid$t2p <- grid_metrics$time_to_peak
  hrf_grid$fwhm <- grid_metrics$FWHM
  cat("Done. Grid now has t2p/fwhm columns.\n")

  n_subjects <- length(fit_allHRFs_result[["subject_results"]])
  n_voxels <- nrow(subject_masks)
  cat("Processing ", n_subjects, " subjects...\n")

  # Step 4: Build best_params_df by extracting RSS and computing best model per winning c
  # First pass: determine winning c using subject masks
  cat("\nStep 4a: Counting c values to determine winner...\n")
  c_votes <- c(c0 = 0, c1 = 0)

  for (i in 1:n_subjects) {
    # Get bestmodel (original, all models)
    bestmodel_L <- fit_allHRFs_result[["subject_results"]][[i]][["glm_result"]][["mGLM0s"]][["cortexL"]][["bestmodel"]]
    bestmodel_R <- fit_allHRFs_result[["subject_results"]][[i]][["glm_result"]][["mGLM0s"]][["cortexR"]][["bestmodel"]]
    bestmodel <- c(bestmodel_L, bestmodel_R)

    # Get subject mask for this subject
    subj_mask <- subject_masks[, i]

    # Count c values where BOTH masks are TRUE
    both_mask <- subj_mask & (1:n_voxels %in% pop_mask_voxels)
    best_c <- hrf_grid$c[bestmodel[both_mask]]
    c_votes["c0"] <- c_votes["c0"] + sum(best_c == 0)
    c_votes["c1"] <- c_votes["c1"] + sum(best_c != 0)
  }

  winning_c <- if (c_votes["c0"] > c_votes["c1"]) 0 else 1/6
  winning_c_indices <- if (winning_c == 0) c0_indices else c1_indices
  cat("C votes: c=0 has ", c_votes["c0"], ", c=1/6 has ", c_votes["c1"], "\n")
  cat("Winning c: ", winning_c, "\n\n")

  # Step 4b: For each subject, compute best model among winning c models using RSS
  cat("Step 4b: Computing best model per voxel among winning c models...\n")

  best_params_list <- list()

  for (i in 1:n_subjects) {
    if (i %% 100 == 0) cat("  Subject ", i, "/", n_subjects, "\n")

    # Get RSS matrices
    RSS_L <- fit_allHRFs_result[["subject_results"]][[i]][["glm_result"]][["mGLM0s"]][["cortexL"]][["RSS"]]
    RSS_R <- fit_allHRFs_result[["subject_results"]][[i]][["glm_result"]][["mGLM0s"]][["cortexR"]][["RSS"]]
    RSS <- rbind(RSS_L, RSS_R)

    # Find best model among winning c models only
    RSS_winning_c <- RSS[, winning_c_indices]
    bestmodel_local <- apply(RSS_winning_c, 1, which.min)
    bestmodel_global <- winning_c_indices[bestmodel_local]

    # Get subject mask
    subj_mask <- subject_masks[, i]

    # Build df for this subject - lookup t2p/fwhm from pre-computed grid
    best_params_list[[i]] <- data.frame(
      voxel = 1:n_voxels,
      subject = i,
      a1 = hrf_grid$a1[bestmodel_global],
      b1 = hrf_grid$b1[bestmodel_global],
      c = hrf_grid$c[bestmodel_global],
      t2p = hrf_grid$t2p[bestmodel_global],
      fwhm = hrf_grid$fwhm[bestmodel_global],
      mask = subj_mask
    )
  }

  best_params_df <- bind_rows(best_params_list)
  cat("Built best_params_df: ", nrow(best_params_df), " rows\n\n")

  # Step 5: Compute population average using both masks
  cat("Step 5: Computing population averages...\n")

  pop_avg <- best_params_df %>%
    filter(
      mask == TRUE,                    # Subject mask
      voxel %in% pop_mask_voxels       # Population mask
    ) %>%
    group_by(voxel) %>%
    summarize(
      t2p_mean = mean(t2p, na.rm = TRUE),
      fwhm_mean = mean(fwhm, na.rm = TRUE),
      n_subjects = n(),
      .groups = "drop"
    )

  cat("Population averages computed for ", nrow(pop_avg), " voxels\n")

  # Step 6: Snap averaged t2p/fwhm back to a1/b1 grid with fixed winning_c
  cat("Step 6: Snapping t2p/fwhm back to a1/b1 grid...\n")
  snapped <- snap_to_grid_t2p_fwhm(pop_avg$t2p_mean, pop_avg$fwhm_mean, hrf_grid, winning_c)
  pop_avg$a1_snapped <- snapped$a1
  pop_avg$b1_snapped <- snapped$b1
  pop_avg$c_snapped <- snapped$c

  cat("Snapped to grid: ", length(unique(paste(pop_avg$a1_snapped, pop_avg$b1_snapped))), " unique (a1, b1) combos\n")

  return(list(
    best_params_df = best_params_df,
    pop_avg = pop_avg,
    winning_c = winning_c,
    hrf_grid = hrf_grid,
    mask_prop_NA = mask_prop_NA
  ))
}

# Run it and save
result <- create_population_avg(fit_results_root_dir, "GAMBLING_LR_1084s")
saveRDS(result, "~/Documents/Github/hrf-z/dev/exploration/fixtures/pop_avg_GAMBLING_LR_1084s.rds")

########################################
##### POC BEGINS #######################
########################################

library(ciftiTools)
source("~/Documents/Github/HRF-Adaptation-paper/Code/0_analysis.R")
ciftiTools.setOption('wb_path', if (Sys.info()[['sysname']] == 'Darwin') '/Applications/workbench/bin_macosxub/wb_command' else '~/workbench/bin_rh_linux64/wb_command')

# snap_to_grid: snap (a1, b1, c) to nearest grid point
snap_to_grid <- function(param_df, hrf_grid, a1_col = "a1", b1_col = "b1", c_col = "c") {
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

# create_candidate_maps_v2: uses pop_avg from create_population_avg
create_candidate_maps_v2 <- function(pop_avg_result) {
  pop_avg <- pop_avg_result$pop_avg
  hrf_grid <- pop_avg_result$hrf_grid

  # Use snapped values as base
  group_map_snapped <- pop_avg %>%
    select(voxel, a1_snapped, b1_snapped, c_snapped)

  # Offsets
  a1_offsets <- c(-2, -1, 0, 1, 2)
  b1_offsets <- c(-0.5, -0.25, 0, 0.25, 0.5)
  offset_combos <- expand.grid(a1_offset = a1_offsets, b1_offset = b1_offsets)
  cat("Total candidate maps:", nrow(offset_combos), "\n")

  candidate_maps <- list()

  for (i in 1:nrow(offset_combos)) {
    a1_off <- offset_combos$a1_offset[i]
    b1_off <- offset_combos$b1_offset[i]

    # Apply offset, truncate to valid ranges
    candidate_map_raw <- group_map_snapped %>%
      mutate(
        a1_offset_applied = pmin(pmax(a1_snapped + a1_off, 3), 12),
        b1_offset_applied = pmin(pmax(b1_snapped + b1_off, 0.5), 2.0),
        c_offset_applied = c_snapped
      )

    # Snap to grid
    candidate_map <- snap_to_grid(
      candidate_map_raw, hrf_grid,
      a1_col = "a1_offset_applied",
      b1_col = "b1_offset_applied",
      c_col = "c_offset_applied"
    ) %>%
      select(voxel, a1_candidate = a1_snapped, b1_candidate = b1_snapped, c_candidate = c_snapped)

    candidate_maps[[i]] <- candidate_map
    candidate_maps[[i]]$offset_id <- i
    candidate_maps[[i]]$a1_offset <- a1_off
    candidate_maps[[i]]$b1_offset <- b1_off
  }

  # Create full HRF grid with model indices
  full_hrf_grid <- do.call(rbind, lapply(candidate_maps, function(cm) {
    cm %>% select(a1_candidate, b1_candidate, c_candidate)
  })) %>%
    distinct() %>%
    mutate(model_idx = row_number())

  # Add t2p/fwhm via lookup (hrf_grid already has these from create_population_avg)
  grid_lookup <- hrf_grid %>% select(a1, b1, c, t2p, fwhm) %>% distinct()
  full_hrf_grid <- full_hrf_grid %>%
    left_join(grid_lookup, by = c("a1_candidate" = "a1", "b1_candidate" = "b1", "c_candidate" = "c"))

  cat("Total unique HRF combos:", nrow(full_hrf_grid), "\n")

  # Add model_idx to each candidate map
  for (i in 1:length(candidate_maps)) {
    candidate_maps[[i]] <- candidate_maps[[i]] %>%
      left_join(full_hrf_grid, by = c("a1_candidate", "b1_candidate", "c_candidate"))
  }

  return(list(candidate_maps = candidate_maps, full_hrf_grid = full_hrf_grid))
}

# fit_maps_v2: returns t2p/fwhm in results
fit_maps_v2 <- function(candidate_maps, full_hrf_grid, session_data, subject_idx) {
  BOLD_xii <- ciftiTools::read_cifti(
    session_data[["BOLD_files"]][[subject_idx]],
    brainstructures = c("left", "right"),
    resamp_res = 10000
  )
  nuisance <- as.matrix(read.table(session_data[["nuisance_files"]][[subject_idx]], header = FALSE))

  n_models <- nrow(full_hrf_grid)
  n_timepoints <- ncol(BOLD_xii)
  n_regressors <- ncol(make_design(
    EVs = session_data[["EVs_list"]][[subject_idx]],
    nTime = n_timepoints, TR = 0.72, dHRF = 0, onset = FALSE, offset = FALSE
  )$design)

  design_canonical <- make_design(
    EVs = session_data[["EVs_list"]][[subject_idx]],
    nTime = n_timepoints, TR = 0.72, dHRF = 0, onset = FALSE, offset = FALSE
  )$design

  design_3D <- array(NA, dim = c(n_timepoints, n_regressors, n_models))
  for (j in 1:n_models) {
    design_3D[,,j] <- make_design(
      EVs = session_data[["EVs_list"]][[subject_idx]],
      nTime = n_timepoints, TR = 0.72, dHRF = 0, onset = FALSE, offset = FALSE,
      a1 = full_hrf_grid$a1_candidate[j],
      b1 = full_hrf_grid$b1_candidate[j],
      c  = full_hrf_grid$c_candidate[j],
      a2 = (16 / sqrt(6)) * sqrt(full_hrf_grid$a1_candidate[j]) * sqrt(full_hrf_grid$b1_candidate[j]),
      b2 = full_hrf_grid$b1_candidate[j]
    )$design
  }

  glm_result <- multiGLM(
    BOLD = BOLD_xii, brainstructures = c("left", "right"),
    design = design_3D, design_canonical = design_canonical,
    scrub = NULL, nuisance = nuisance, TR = 0.72, hpf = 0.01, resamp_res = NULL
  )

  RSS_matrix <- rbind(glm_result$mGLM0s$cortexL$RSS, glm_result$mGLM0s$cortexR$RSS)
  Fstat_vector <- c(glm_result$mGLM0s$cortexL$Fstat, glm_result$mGLM0s$cortexR$Fstat)

  total_RSS_per_candidate <- numeric(length(candidate_maps))
  weighted_RSS_per_candidate <- numeric(length(candidate_maps))

  for (i in 1:length(candidate_maps)) {
    model_indices <- candidate_maps[[i]]$model_idx
    voxel_rss <- RSS_matrix[cbind(1:nrow(RSS_matrix), model_indices)]
    total_RSS_per_candidate[i] <- sum(voxel_rss, na.rm = TRUE)
    # V2.2: Use Fstat^1 weighting (milder than ^2)
    weighted_RSS_per_candidate[i] <- sum(voxel_rss * Fstat_vector^1, na.rm = TRUE) /
      sum(Fstat_vector^1, na.rm = TRUE)
  }

  # Compute mean t2p/fwhm per candidate (averaged across voxels)
  t2p_per_candidate <- numeric(length(candidate_maps))
  fwhm_per_candidate <- numeric(length(candidate_maps))

  for (i in 1:length(candidate_maps)) {
    t2p_per_candidate[i] <- mean(candidate_maps[[i]]$t2p, na.rm = TRUE)
    fwhm_per_candidate[i] <- mean(candidate_maps[[i]]$fwhm, na.rm = TRUE)
  }

  # Build results with t2p/fwhm
  results_df <- data.frame(
    candidate_id = 1:length(candidate_maps),
    a1_offset = sapply(candidate_maps, function(x) x$a1_offset[1]),
    b1_offset = sapply(candidate_maps, function(x) x$b1_offset[1]),
    t2p = t2p_per_candidate,
    fwhm = fwhm_per_candidate,
    total_RSS = total_RSS_per_candidate,
    weighted_RSS = weighted_RSS_per_candidate
  ) %>%
    arrange(weighted_RSS)

  return(results_df)
}

########################################
##### MAIN LOOP ########################
########################################

run_names <- c('LR','RL')
# Change to 1:1 to run only LR, or 2:2 to run only RL
runs_to_process <- 1:1

# HPC paths
output_dir <- "/N/project/hrf_adaptation/validation_bonferroni_0.05_smoothed_RSS/POCv2.2"
fixture_dir <- "~/Documents/Github/hrf-z/dev/exploration/fixtures"

for(r in runs_to_process){
  sess <- paste("tfMRI_GAMBLING", run_names[r], sep='_')
  sess_name <- paste0(sess, "_1084s")  # Using 1084 subjects for population

  cat("\n********** Processing", sess, "**********\n"); flush.console()

  # Load pop_avg from fixtures (this is the v2 difference: uses t2p/fwhm averaging)
  pop_avg_file <- file.path(fixture_dir, paste0("pop_avg_", sess_name, ".rds"))
  cat("Loading pop_avg from", pop_avg_file, "\n"); flush.console()
  pop_avg_result <- readRDS(pop_avg_file)

  # Create candidate maps (v2: based on pop_avg instead of regularize_result)
  cat("Creating candidate maps...\n"); flush.console()
  cm_result <- create_candidate_maps_v2(pop_avg_result)

  # Process both test_matched and retest for this run
  data_types <- c("test_matched", "retest")

  for(data_type in data_types){
    cat("\n********** Processing", sess, toupper(data_type), "**********\n"); flush.console()

    # Load appropriate session data
    if(data_type == "test_matched"){
      session_data <- readRDS(file.path(fixture_dir, paste0('test_matched_session_data_', sess, '.rds')))
    } else {
      session_data <- readRDS(file.path('~/Documents/Github/HRF-Adaptation-paper/Code/fixtures', paste0('retest_session_data_', sess, '.rds')))
    }

    # Collect results
    results_list <- list()
    n_subjects <- length(session_data$BOLD_files)
    cat("Testing on", n_subjects, "subjects\n"); flush.console()

    for (subject_idx in 1:n_subjects) {
      cat(toupper(data_type), "subject", subject_idx, "/", n_subjects, "\n"); flush.console()

      results_list[[subject_idx]] <- fit_maps_v2(
        cm_result$candidate_maps,
        cm_result$full_hrf_grid,
        session_data,
        subject_idx
      ) %>%
        mutate(subject_id = subject_idx)
    }

    # Combine all results
    results_combined <- bind_rows(results_list)

    # Save complete results with proper naming
    output_file <- safe_path(output_dir, paste0("results_POCv2.2_", data_type, "_gambling_", run_names[r], "_weighted.rds"))
    saveRDS(results_combined, output_file)
    cat("\n*** Saved", data_type, "results to:", output_file, "***\n"); flush.console()

    rm(results_list, results_combined, session_data); gc()
  }

  rm(pop_avg_result, cm_result); gc()
}

cat("\n========================================\n")
cat("POCv2.2 complete!\n")
cat("========================================")
