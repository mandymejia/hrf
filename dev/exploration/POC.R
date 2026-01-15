devtools::load_all("~/Documents/Github/hrf-z")
library(dplyr)
library(here)
library(ciftiTools)
ciftiTools.setOption('wb_path', if (Sys.info()[['sysname']] == 'Darwin') '/Applications/workbench/bin_macosxub/wb_command' else '~/workbench/bin_rh_linux64/wb_command')

snap_to_grid <- function(param_df, hrf_grid,
                           a1_col = "a1_pop",
                           b1_col = "b1_pop",
                           c_col = "c_pop") {

  # Get unique grid combinations
  grid_combos <- hrf_grid %>%
    select(a1, b1, c) %>%
    distinct()

  # Calculate ranges for normalization
  a_range <- max(grid_combos$a1) - min(grid_combos$a1)
  b_range <- max(grid_combos$b1) - min(grid_combos$b1)

  # For each row, find nearest grid point
  snapped_df <- param_df %>%
    rowwise() %>%
    mutate(
      grid_idx = {
        # Get current parameter values
        a1_val <- .data[[a1_col]]
        b1_val <- .data[[b1_col]]
        c_val <- .data[[c_col]]

        # Find grid points with matching c
        c_match <- abs(grid_combos$c - c_val) < 1e-6

        if (sum(c_match) == 0) {
          # No exact match, find closest c
          c_dists <- abs(grid_combos$c - c_val)
          c_match <- c_dists == min(c_dists)
        }

        # Among matching c, find closest (a1, b1) using normalized distance
        candidates <- grid_combos[c_match, ]
        distances <- sqrt(
          ((a1_val - candidates$a1) / a_range)^2 +
            ((b1_val - candidates$b1) / b_range)^2
        )

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

create_candidate_maps <- function(regularize_result) {
  a1_pop <- regularize_result[["regularized_params"]][["a1"]][["results_OLS"]] %>%
    select(voxel, a1_pop = param_mean0) %>%
    distinct()

  b1_pop <- regularize_result[["regularized_params"]][["b1"]][["results_OLS"]] %>%
    select(voxel, b1_pop = param_mean0) %>%
    distinct()

  c_pop <- regularize_result[["regularized_params"]][["c"]][["results_OLS"]] %>%
    select(voxel, c_pop = param_mean0) %>%
    distinct()

  group_map <- a1_pop %>%
    left_join(b1_pop, by = "voxel") %>%
    left_join(c_pop, by = "voxel")

  hrf_grid <- regularize_result[["hrf_grid"]]
  group_map_snapped <- snap_to_grid(group_map, hrf_grid,
                                    a1_col = "a1_pop",
                                    b1_col = "b1_pop",
                                    c_col = "c_pop")
  # Step 1 Create offsets
  a1_offsets <- c(-2, -1, 0, 1, 2)
  b1_offsets <- c(-0.5, -0.25, 0, 0.25, 0.5)
  # c is fixed (no offset, per Mandy's suggestion)

  # Create all combinations
  offset_combos <- expand.grid(a1_offset = a1_offsets, b1_offset = b1_offsets)
  cat("Total candidate maps:", nrow(offset_combos), "\n")

  # Step 2 create candidate maps
  # Store all candidate maps in a list
  candidate_maps <- list()

  for (i in 1:nrow(offset_combos)) {
    a1_off <- offset_combos$a1_offset[i]
    b1_off <- offset_combos$b1_offset[i]

    # Apply global offset to all voxels
    candidate_map_raw <- group_map_snapped %>%
      mutate(
        a1_offset_applied = a1_snapped + a1_off,
        b1_offset_applied = b1_snapped + b1_off,
        c_offset_applied = c_snapped  # Fixed, no offset
      ) %>%
      # Truncate to valid ranges FIRST
      mutate(
        a1_offset_applied = pmin(pmax(a1_offset_applied, 3), 12),
        b1_offset_applied = pmin(pmax(b1_offset_applied, 0.5), 2.0)
      )

    # NOW snap to grid
    candidate_map <- snap_to_grid(
      candidate_map_raw,
      hrf_grid,
      a1_col = "a1_offset_applied",
      b1_col = "b1_offset_applied",
      c_col = "c_offset_applied"
    ) %>%
      select(voxel, a1_candidate = a1_snapped, b1_candidate = b1_snapped, c_candidate = c_snapped)

    # Store it
    candidate_maps[[i]] <- candidate_map
    candidate_maps[[i]]$offset_id <- i
    candidate_maps[[i]]$a1_offset <- a1_off
    candidate_maps[[i]]$b1_offset <- b1_off
  }

  # Step 3: Create full HRF grid and add global model indices to candidate maps

  # Extract ALL unique combos across all candidate maps
  full_hrf_grid <- do.call(rbind, lapply(candidate_maps, function(cm) {
    cm %>% select(a1_candidate, b1_candidate, c_candidate)
  })) %>%
    distinct() %>%
    mutate(model_idx = row_number())

  cat("Total unique HRF combos across all candidate maps:", nrow(full_hrf_grid), "\n")

  # For each candidate map, add model_idx by matching to full_hrf_grid
  for (i in 1:length(candidate_maps)) {
    candidate_maps[[i]] <- candidate_maps[[i]] %>%
      left_join(full_hrf_grid, by = c("a1_candidate", "b1_candidate", "c_candidate"))
  }

  return(list(
    candidate_maps = candidate_maps,
    full_hrf_grid = full_hrf_grid
  ))

}


fit_maps <- function(candidate_maps, full_hrf_grid, session_data, subject_idx) {
  # Step 4 make design matrix and run multiGLM
  # Load BOLD and nuisance
  BOLD_xii <- ciftiTools::read_cifti(
    session_data[["BOLD_files"]][[subject_idx]],
    brainstructures = c("left", "right"),
    resamp_res = 10000
  )
  nuisance <- as.matrix(read.table(session_data[["nuisance_files"]][[subject_idx]], header = FALSE))

  # Get dimensions for design_3D
  n_models <- nrow(full_hrf_grid)
  n_timepoints <- ncol(BOLD_xii)
  n_regressors <- ncol(make_design(
    EVs = session_data[["EVs_list"]][[subject_idx]],
    nTime = n_timepoints,
    TR = 0.72,
    dHRF = 0,
    onset = FALSE, # FALSE FOR GAMBLING
    offset = FALSE
  )$design)

  # Build design_3D
  design_3D <- array(NA, dim = c(n_timepoints, n_regressors, n_models))

  for (j in 1:n_models) {
    design_3D[,,j] <- make_design(
      EVs = session_data[["EVs_list"]][[subject_idx]],
      nTime = n_timepoints,
      TR = 0.72,
      dHRF = 0,
      onset = FALSE,
      offset = FALSE,
      a1 = full_hrf_grid$a1_candidate[j],
      b1 = full_hrf_grid$b1_candidate[j],
      c  = full_hrf_grid$c_candidate[j],
      a2 = (16 / sqrt(6)) * sqrt(full_hrf_grid$a1_candidate[j]) * sqrt(full_hrf_grid$b1_candidate[j]),
      b2 = full_hrf_grid$b1_candidate[j]
    )$design
  }

  # Run multiGLM once
  glm_result <- multiGLM(
    BOLD = BOLD_xii,
    brainstructures = c("left", "right"),
    design = design_3D,
    scrub = NULL,
    nuisance = nuisance,
    TR = 0.72,
    hpf = 0.01,
    resamp_res = NULL
  )

  # Extract RSS matrix
  RSS_left <- glm_result$mGLM0s$cortexL$RSS
  RSS_right <- glm_result$mGLM0s$cortexR$RSS
  RSS_matrix <- rbind(RSS_left, RSS_right)

  cat("RSS matrix dimensions:", dim(RSS_matrix), "\n")


  # Step 5: For each candidate map, compute total RSS
  total_RSS_per_candidate <- numeric(length(candidate_maps))

  for (i in 1:length(candidate_maps)) {
    # Get model indices for this candidate map's voxels
    model_indices <- candidate_maps[[i]]$model_idx

    # Extract RSS for each voxel's assigned model
    # Use matrix indexing: RSS_matrix[row, col] where row=voxel, col=model
    voxel_rss <- RSS_matrix[cbind(1:nrow(RSS_matrix), model_indices)]

    # Sum across all voxels
    total_RSS_per_candidate[i] <- sum(voxel_rss, na.rm = TRUE)

    cat("Candidate", i, "(a1=", candidate_maps[[i]]$a1_offset[1],
        ", b1=", candidate_maps[[i]]$b1_offset[1], "): RSS =",
        total_RSS_per_candidate[i], "\n")
  }

  # Find best offset
  best_idx <- which.min(total_RSS_per_candidate)

  cat("\n***BEST OFFSET***\n")
  cat("Candidate Map ID:", best_idx, "\n")
  cat("Offset ID:", candidate_maps[[best_idx]]$offset_id[1], "\n")
  cat("a1 offset:", candidate_maps[[best_idx]]$a1_offset[1], "\n")
  cat("b1 offset:", candidate_maps[[best_idx]]$b1_offset[1], "\n")
  cat("Minimum total RSS:", total_RSS_per_candidate[best_idx], "\n")

  # Show all results sorted
  results_df <- data.frame(
    candidate_id = 1:length(candidate_maps),
    offset_id = sapply(candidate_maps, function(x) x$offset_id[1]),
    a1_offset = sapply(candidate_maps, function(x) x$a1_offset[1]),
    b1_offset = sapply(candidate_maps, function(x) x$b1_offset[1]),
    total_RSS = total_RSS_per_candidate
  ) %>%
    arrange(total_RSS)

  return(results_df)
}


regularize_path <- "/Volumes/LaCie/root/hrf_adaptation/validation/regularize_allHRFs/regularize_allHRFs_result_tfMRI_GAMBLING_LR_1084s.rds"
session_data <- readRDS(here("dev", "fixtures", "session_data_4s", paste0("session_data_", "motor", "_lr_4s.rds")))

regularize_result <- readRDS(regularize_path)

cm_result <- create_candidate_maps(regularize_result)
results <- fit_maps(cm_result$candidate_maps, cm_result$full_hrf_grid, session_data, 1)

head(results)



