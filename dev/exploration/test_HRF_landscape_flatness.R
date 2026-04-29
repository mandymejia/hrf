library(here)
library(dplyr)
library(ciftiTools)

ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z")

# Load session data (mount hcp_dcwan first: sshfs quartz:/N/project/hcp_dcwan ~/mnt/hcp_dcwan)
sess <- "tfMRI_GAMBLING_LR"
session_data <- readRDS(file.path('~/Documents/Github/HRF-Adaptation-paper/Code/fixtures', paste0('session_data_', sess, '.rds')))

# Remap paths from HPC to mounted directory
remap_path <- function(hpc_path) {
  hpc_path <- gsub("^/N/project/hrf_adaptation", "~/mnt/hrf_adaptation", hpc_path)
  hpc_path <- gsub("^/N/project/hcp_dcwan", "~/mnt/hcp_dcwan", hpc_path)
  return(hpc_path)
}

session_data$BOLD_files <- lapply(session_data$BOLD_files, remap_path)
session_data$nuisance_files <- lapply(session_data$nuisance_files, remap_path)

subject_idx <- 1  # First subject

cat("Loading subject", subject_idx, "data...\n")
BOLD_xii <- ciftiTools::read_cifti(
  session_data[["BOLD_files"]][[subject_idx]],
  brainstructures = c("left", "right"),
  resamp_res = 10000
)
nuisance <- as.matrix(read.table(session_data[["nuisance_files"]][[subject_idx]], header = FALSE))

# Load HRF grid from package
data("default_hrf_grid")
hrf_grid <- default_hrf_grid

n_models <- nrow(hrf_grid)
n_timepoints <- ncol(BOLD_xii)
n_regressors <- ncol(make_design(
  EVs = session_data[["EVs_list"]][[subject_idx]],
  nTime = n_timepoints,
  TR = 0.72,
  dHRF = 0,
  onset = FALSE,
  offset = FALSE
)$design)

cat("Building design matrix for", n_models, "HRF models...\n")

# Create canonical design
design_canonical <- make_design(
  EVs = session_data[["EVs_list"]][[subject_idx]],
  nTime = n_timepoints,
  TR = 0.72,
  dHRF = 0,
  onset = FALSE,
  offset = FALSE
)$design

# Build design_3D with all HRFs from grid
design_3D <- array(NA, dim = c(n_timepoints, n_regressors, n_models))

for (j in 1:n_models) {
  design_3D[,,j] <- make_design(
    EVs = session_data[["EVs_list"]][[subject_idx]],
    nTime = n_timepoints,
    TR = 0.72,
    dHRF = 0,
    onset = FALSE,
    offset = FALSE,
    a1 = hrf_grid$a1[j],
    b1 = hrf_grid$b1[j],
    c  = hrf_grid$c[j],
    a2 = hrf_grid$a2[j],
    b2 = hrf_grid$b2[j]
  )$design
}

cat("Running multiGLM...\n")
glm_result <- multiGLM(
  BOLD = BOLD_xii,
  brainstructures = c("left", "right"),
  design = design_3D,
  design_canonical = design_canonical,
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

cat("\n========================================\n")
cat("RSS LANDSCAPE ANALYSIS (fit_allHRFs stage)\n")
cat("========================================\n\n")

# For each voxel, analyze RSS spread across all HRF models
voxel_stats <- data.frame(
  voxel = 1:nrow(RSS_matrix),
  min_rss = apply(RSS_matrix, 1, min, na.rm = TRUE),
  max_rss = apply(RSS_matrix, 1, max, na.rm = TRUE),
  mean_rss = apply(RSS_matrix, 1, mean, na.rm = TRUE)
) %>%
  mutate(
    rss_range = max_rss - min_rss,
    rss_range_pct = 100 * rss_range / min_rss
  )

cat("Across all", nrow(RSS_matrix), "voxels:\n\n")
cat("Mean RSS range (best to worst HRF):", round(mean(voxel_stats$rss_range), 6), "\n")
cat("Median RSS range:", round(median(voxel_stats$rss_range), 6), "\n")
cat("\n")
cat("Mean RSS range as % of minimum:", round(mean(voxel_stats$rss_range_pct), 2), "%\n")
cat("Median RSS range as % of minimum:", round(median(voxel_stats$rss_range_pct), 2), "%\n")
cat("\n")

# Show distribution
cat("Distribution of RSS range %:\n")
cat("  < 5%:", sum(voxel_stats$rss_range_pct < 5), "voxels\n")
cat("  5-10%:", sum(voxel_stats$rss_range_pct >= 5 & voxel_stats$rss_range_pct < 10), "voxels\n")
cat("  10-20%:", sum(voxel_stats$rss_range_pct >= 10 & voxel_stats$rss_range_pct < 20), "voxels\n")
cat("  > 20%:", sum(voxel_stats$rss_range_pct >= 20), "voxels\n")
cat("\n")

# Sample a few voxels
cat("Sample of 5 random voxels:\n")
set.seed(42)
sample_voxels <- sample(1:nrow(RSS_matrix), 5)
for (v in sample_voxels) {
  cat(sprintf("Voxel %d: min RSS = %.4f, max RSS = %.4f, range = %.4f (%.1f%%)\n",
              v, voxel_stats$min_rss[v], voxel_stats$max_rss[v],
              voxel_stats$rss_range[v], voxel_stats$rss_range_pct[v]))
}

cat("\n")
cat("INTERPRETATION:\n")
cat("===============\n")
cat("If most voxels have range < 10%, the RSS landscape is FLAT at fit_allHRFs.\n")
cat("This means the BOLD signal lacks sufficient information to distinguish HRF shapes.\n")
cat("If this is true, NO amount of regularization or POC tricks can fix the problem -\n")
cat("the limitation is fundamental to the SNR of fMRI data.\n")
cat("\n========================================\n")
