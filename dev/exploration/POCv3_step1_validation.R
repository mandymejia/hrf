library(dplyr)
devtools::load_all("~/Documents/Github/hrf-z")

# HPC path
fit_results_root_dir <- "/N/project/hrf_adaptation/validation_bonferroni_0.05_smoothed_RSS"
fixture_dir <- "~/Documents/Github/hrf-z/dev/exploration/fixtures"

########################################
##### INVERSE LOOKUP (from step0) ######
########################################

# Inverse lookup: given (t2p, fwhm) -> find continuous (a1, b1)
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
##### HRF METRICS ######################
########################################

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

# NOTE: snap_to_grid_t2p_fwhm removed in v3 - using inverse_lookup instead

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

  # Step 6: Convert averaged t2p/fwhm back to CONTINUOUS a1/b1 (NO SNAPPING in v3)
  cat("Step 6: Converting t2p/fwhm to continuous a1/b1 using dense lookup...\n")
  dense_lookup <- readRDS(file.path(fixture_dir, "dense_hrf_lookup_c16_1dec.rds"))
  continuous <- inverse_lookup_batch(pop_avg$t2p_mean, pop_avg$fwhm_mean, dense_lookup)
  pop_avg$a1_continuous <- continuous$a1
  pop_avg$b1_continuous <- continuous$b1
  pop_avg$c_continuous <- winning_c

  cat("Continuous a1/b1 computed (no grid snapping). Unique a1 values: ",
      length(unique(pop_avg$a1_continuous)), "\n")

  return(list(
    best_params_df = best_params_df,
    pop_avg = pop_avg,
    winning_c = winning_c,
    hrf_grid = hrf_grid,
    mask_prop_NA = mask_prop_NA
  ))
}

# Run it and save (v3 uses continuous a1/b1)
result <- create_population_avg(fit_results_root_dir, "GAMBLING_LR_1084s")
saveRDS(result, file.path(fixture_dir, "pop_avg_v3.1_tfMRI_GAMBLING_LR_1084s.rds"))

########################################
##### POC BEGINS #######################
########################################

library(ciftiTools)
source("~/Documents/Github/HRF-Adaptation-paper/Code/0_analysis.R")
ciftiTools.setOption('wb_path', if (Sys.info()[['sysname']] == 'Darwin') '/Applications/workbench/bin_macosxub/wb_command' else '~/workbench/bin_rh_linux64/wb_command')

# NOTE: snap_to_grid removed in v3 - no longer snapping after offsets

# create_candidate_maps_v3: uses CONTINUOUS a1/b1 (NO SNAPPING)
create_candidate_maps_v3 <- function(pop_avg_result) {
  pop_avg <- pop_avg_result$pop_avg
  winning_c <- pop_avg_result$winning_c

  # Use CONTINUOUS values as base (not snapped)
  group_map <- pop_avg %>%
    select(voxel, a1_continuous, b1_continuous, c_continuous)

  # Offsets (more granular in v3)
  a1_offsets <- c(-2, -1, 0, 1, 2)
  b1_offsets <- c(-0.5, -0.25, 0, 0.25, 0.5)
  offset_combos <- expand.grid(a1_offset = a1_offsets, b1_offset = b1_offsets)
  cat("Total candidate maps:", nrow(offset_combos), "\n")

  # Load dense lookup for t2p/fwhm computation
  dense_lookup <- readRDS(file.path(fixture_dir, "dense_hrf_lookup_c16_1dec.rds"))

  candidate_maps <- list()

  for (i in 1:nrow(offset_combos)) {
    a1_off <- offset_combos$a1_offset[i]
    b1_off <- offset_combos$b1_offset[i]

    # Apply offset, truncate to valid ranges (NO SNAPPING)
    candidate_map <- group_map %>%
      mutate(
        a1_candidate = pmin(pmax(a1_continuous + a1_off, 4), 11),
        b1_candidate = pmin(pmax(b1_continuous + b1_off, 0.5), 2.0),
        c_candidate = c_continuous
      ) %>%
      select(voxel, a1_candidate, b1_candidate, c_candidate)

    candidate_maps[[i]] <- candidate_map
    candidate_maps[[i]]$offset_id <- i
    candidate_maps[[i]]$a1_offset <- a1_off
    candidate_maps[[i]]$b1_offset <- b1_off
  }

  # Create full HRF grid with model indices (all unique continuous combos)
  full_hrf_grid <- do.call(rbind, lapply(candidate_maps, function(cm) {
    cm %>% select(a1_candidate, b1_candidate, c_candidate)
  })) %>%
    distinct() %>%
    mutate(model_idx = row_number())

  # Compute t2p/fwhm for continuous values using dense lookup
  cat("Computing t2p/fwhm for", nrow(full_hrf_grid), "unique HRF combos...\n")
  full_hrf_grid$t2p <- full_hrf_grid$a1_candidate - full_hrf_grid$b1_candidate

  # For fwhm, find nearest match in dense lookup
  fwhm_vals <- numeric(nrow(full_hrf_grid))
  for (j in 1:nrow(full_hrf_grid)) {
    # Find closest (a1, b1) in dense lookup
    dist <- sqrt(((dense_lookup$grid$a1 - full_hrf_grid$a1_candidate[j]))^2 +
                 ((dense_lookup$grid$b1 - full_hrf_grid$b1_candidate[j]))^2)
    nearest_idx <- which.min(dist)
    fwhm_vals[j] <- dense_lookup$grid$fwhm[nearest_idx]
  }
  full_hrf_grid$fwhm <- fwhm_vals

  cat("Total unique HRF combos:", nrow(full_hrf_grid), "\n")

  # Add model_idx and t2p/fwhm to each candidate map
  for (i in 1:length(candidate_maps)) {
    candidate_maps[[i]] <- candidate_maps[[i]] %>%
      left_join(full_hrf_grid, by = c("a1_candidate", "b1_candidate", "c_candidate"))
  }

  return(list(candidate_maps = candidate_maps, full_hrf_grid = full_hrf_grid))
}

# fit_maps_v3: returns t2p/fwhm in results (same as v2, works with continuous values)
fit_maps_v3 <- function(candidate_maps, full_hrf_grid, session_data, subject_idx) {
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
    weighted_RSS_per_candidate[i] <- sum(voxel_rss * Fstat_vector^2, na.rm = TRUE) /
      sum(Fstat_vector^2, na.rm = TRUE)
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
output_dir <- "/N/project/hrf_adaptation/validation_bonferroni_0.05_smoothed_RSS/POCv3"

for(r in runs_to_process){
  sess <- paste("tfMRI_GAMBLING", run_names[r], sep='_')
  sess_name <- paste0(sess, "_1084s")  # Using 1084 subjects for population

  cat("\n********** Processing", sess, "**********\n"); flush.console()

  # Load pop_avg from fixtures (v3: uses CONTINUOUS a1/b1)
  pop_avg_file <- file.path(fixture_dir, paste0("pop_avg_v3.1_", sess_name, ".rds"))
  cat("Loading pop_avg from", pop_avg_file, "\n"); flush.console()
  pop_avg_result <- readRDS(pop_avg_file)

  # Create candidate maps (v3: NO SNAPPING, continuous values)
  cat("Creating candidate maps (v3: no snapping)...\n"); flush.console()
  cm_result <- create_candidate_maps_v3(pop_avg_result)

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

      results_list[[subject_idx]] <- fit_maps_v3(
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
    output_file <- safe_path(output_dir, paste0("results_POCv3_", data_type, "_gambling_", run_names[r], "_weighted.rds"))
    saveRDS(results_combined, output_file)
    cat("\n*** Saved", data_type, "results to:", output_file, "***\n"); flush.console()

    rm(results_list, results_combined, session_data); gc()
  }

  rm(pop_avg_result, cm_result); gc()
}

cat("\n========================================\n")
cat("POCv3 complete!\n")
cat("========================================")
