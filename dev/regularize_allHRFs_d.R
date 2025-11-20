library(here)

ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z") # hrf-z is Zeshawn's branch with latest hrf package changes.
# devtools::load_all("~/Documents/Github/hrf-HRFcalc-mods") # Temporary while cleaning up hrf-z branch

workingHRF_results <- readRDS(here("dev", "fixtures", "fit_workingHRF_result_motorlr_4s.rds"))
allHRF_results <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_4s.rds"))

# workingHRF_results <- readRDS("/Users/zeshawnzahid/Downloads/19s_fit_workingHRF_result.rds")
# allHRF_results <- readRDS("/Users/zeshawnzahid/Downloads/19s_fit_allHRFs_result.rds")

workingHRF_results <- readRDS(here("dev", "fixtures", "fit_workingHRF_result_motorlr_500s.rds"))
allHRF_results <- readRDS(here("dev", "fixtures", "fit_allHRFs_result_motorlr_500s.rds"))

allHRF_results <- readRDS("~/Downloads/fit_allHRFs_result_tfMRI_MOTOR_LR_1000s.rds")
workingHRF_results <- readRDS("~/Downloads/fit_workingHRF_result_tfMRI_MOTOR_LR_1000s.rds")

allHRF_results <- readRDS("/Users/zeshawnzahid/sshfs/hrf_adaptation/validation/fit_allHRFs/fit_allHRFs_result_tfMRI_MOTOR_LR_25s.rds")
workingHRF_results <- readRDS("/Users/zeshawnzahid/sshfs/hrf_adaptation/validation/fit_workingHRF/fit_workingHRF_result_tfMRI_MOTOR_LR_25s.rds")

regularize_allHRFs_result <- regularize_allHRFs(
  workingHRF_results,
  allHRF_results
)

regularize_allHRFs_result <- regularize_allHRFs(
  fit_workingHRF_result,
  fit_allHRFs_result
)

saveRDS(
  regularize_allHRFs_result,
  file = here("dev", "fixtures", "regularize_allHRFs_result_motorlr_500s.rds")
)

regularize_allHRFs_result <- readRDS(here("dev", "fixtures", "regularize_allHRFs_result_motorlr_4s.rds"))
base_dir <- "/Volumes/LaCie/root/hrf_adaptation"
regularize_allHRFs_result <- readRDS(file.path(base_dir, "validation/regularize_allHRFs/regularize_allHRFs_result_tfMRI_WM_LR_1083s.rds"))




plot(regularize_allHRFs_result, type = "variance", param = "a1", model = "IS", method = "OLS", fname = here("dev", "test_plots", "regularize_allHRFs", "OLS_A_a1"))
plot(regularize_allHRFs_result, type = "variance", param = "a1", model = "IS", method = "WLS", fname = here("dev", "test_plots", "regularize_allHRFs", "WLS_A_a1"))

plot(regularize_allHRFs_result, type = "wls_weights", param = "a1", model = "IS", fname = here("dev", "test_plots", "regularize_allHRFs", "OLS_A_a1_precision"))

plot(regularize_allHRFs_result, type = "mean", param = "b1", fname = here("dev", "test_plots", "regularize_allHRFs", "mean_b1.png"))
plot(regularize_allHRFs_result, type = "mean", param = "a1", fname = here("dev", "test_plots", "regularize_allHRFs", "mean_a1.png"))
plot(regularize_allHRFs_result, type = "mean", param = "c", fname = here("dev", "test_plots", "regularize_allHRFs", "mean_c.png"))


plot(regularize_allHRFs_result,  type = "mean_all",  param = "b1",  fname = here("dev", "test_plots", "regularize_allHRFs", "mean_b1_filtered.png"))


a <- plot(regularize_allHRFs_result,  type = "param_heatmap")
ggsave("dev/test_plots/regularize_allHRFs/heatmap_IO.png", a, width = 14, height = 7)
plot(regularize_allHRFs_result,  type = "param_heatmap", model = "intercept_slope")
plot(regularize_allHRFs_result,  type = "param_heatmap", model = "intercept_only")
plot(regularize_allHRFs_result,  type = "param_heatmap", model = "best_params")
plot(regularize_allHRFs_result,  type = "param_heatmap", model = "average")

###############################################################################
# DEBUG AREA                                                                  #
###############################################################################

var_range_A <- range(regularize_allHRFs_result[["regularized_params"]][["a1"]][["results_OLS"]][["variance_A"]], na.rm = TRUE)
var_range_B <- range(regularize_allHRFs_result[["regularized_params"]][["a1"]][["results_OLS"]][["variance_B"]], na.rm = TRUE)
range(c(var_range_A, var_range_B))


plot(regularize_allHRFs_result[["regularized_params"]][["a1"]][["residual_variance_xii_A_OLS"]],
     title = 'Variance of Model A Predictions for a1 using OLS method',
     color_mode = 'sequential',
     zlim = range(c(regularize_allHRFs_result[["regularized_params"]][["a1"]][["residual_variance_xii_A_OLS"]]$data,
                    regularize_allHRFs_result[["regularized_params"]][["a1"]][["residual_variance_xii_B_OLS"]]$data),
                  na.rm = TRUE),
     fname = NULL,
     legend_fname = NULL)





snap_to_grid <- function(pop_avg_df, hrf_grid) {

  # Get unique grid combinations
  grid_combos <- hrf_grid %>%
    dplyr::select(a1, b1, c) %>%
    dplyr::distinct()

  # Calculate ranges for normalization
  a_range <- max(grid_combos$a1) - min(grid_combos$a1)
  b_range <- max(grid_combos$b1) - min(grid_combos$b1)

  # For each voxel, find nearest grid point
  snapped_df <- pop_avg_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      grid_idx = {
        # Find grid points with same c value (or closest c)
        c_match <- abs(grid_combos$c - c) < 1e-6

        if (sum(c_match) == 0) {
          # If no exact c match, find closest c
          c_dists <- abs(grid_combos$c - c)
          c_match <- c_dists == min(c_dists)
        }

        # Among matching c values, find closest (a1, b1) using normalized Euclidean distance
        candidates <- grid_combos[c_match, ]
        distances <- sqrt(
          ((a1 - candidates$a1) / a_range)^2 +
            ((b1 - candidates$b1) / b_range)^2
        )

        # Return index of nearest point
        which(c_match)[which.min(distances)]
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      a1 = grid_combos$a1[grid_idx],
      b1 = grid_combos$b1[grid_idx],
      c = grid_combos$c[grid_idx]
    ) %>%
    dplyr::select(-grid_idx)

  return(snapped_df)
}


# Updated build function that includes snapping
build_pop_avg_df <- function(regularize_result, snap = TRUE) {

  # Extract param_mean0 for each parameter
  a1_df <- regularize_result[["regularized_params"]][["a1"]][["results_OLS"]] %>%
    dplyr::select(voxel, a1 = param_mean0) %>%
    dplyr::distinct()

  b1_df <- regularize_result[["regularized_params"]][["b1"]][["results_OLS"]] %>%
    dplyr::select(voxel, b1 = param_mean0) %>%
    dplyr::distinct()

  c_df <- regularize_result[["regularized_params"]][["c"]][["results_OLS"]] %>%
    dplyr::select(voxel, c = param_mean0) %>%
    dplyr::distinct()

  # Merge by voxel
  pop_avg_df <- a1_df %>%
    dplyr::left_join(b1_df, by = "voxel") %>%
    dplyr::left_join(c_df, by = "voxel")

  # Optionally snap to grid
  if (snap) {
    hrf_grid <- regularize_result[["hrf_grid"]]
    if (is.null(hrf_grid)) {
      warning("hrf_grid not found in regularize_result, cannot snap to grid")
    } else {
      pop_avg_df <- snap_to_grid(pop_avg_df, hrf_grid)
    }
  }

  return(pop_avg_df)
}

pop_avg_df <- build_pop_avg_df(regularize_allHRFs_result)



