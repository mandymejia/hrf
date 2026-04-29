# HPC sample: WM 2bk vs 0bk contrast using the 500-subject population HRF map.
# Tested on subject 501 (held out from the 500-subject population).
#
# Contrast: -0.25 each on the four 0bk tasks, +0.25 each on the four 2bk tasks.
#
# See dev/hpc/README.md for setup (RStudio launch, dependency install).

library(ciftiTools)
devtools::load_all("~/Documents/Github/hrf-z")

ciftiTools.setOption('wb_path', '~/workbench/bin_rh_linux64/wb_command')

# ** Paths **
run_name <- "validation_final"
output_root <- file.path("/N/project/hrf_adaptation", run_name)
sess <- "tfMRI_WM_LR"
n_use <- 500  # population built from the first 500 subjects

regularize_path <- file.path(output_root, "regularize_allHRFs",
  paste0("regularize_allHRFs_result_", sess, "_", n_use, "s.rds"))
session_path <- file.path("~/Documents/Github/HRF-Adaptation-paper/Code/fixtures",
  paste0("session_data_", sess, ".rds"))

plot_dir <- file.path("~/Documents/Github/hrf-z/dev/hpc/plots", sess, "2bk_vs_0bk")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ** Load **
cat("Loading regularize result...\n")
regularize_result <- readRDS(regularize_path)

cat("Loading session_data...\n")
session_data <- readRDS(session_path)

# ** Build the contrast matrix **
# WM EV order (HCP convention): four 0bk tasks first, then four 2bk tasks.
ev_names <- names(session_data$EVs_list[[1]])
cat("EV names:", paste(ev_names, collapse = ", "), "\n")

n_main <- length(ev_names)
n_task <- n_main + 2  # main EVs + 1 shared onset + 1 shared offset
A <- matrix(0, nrow = 1, ncol = n_task)
colnames(A) <- c(ev_names, "onset", "offset")
rownames(A) <- "2bk_minus_0bk"

A[1, grep("^0bk", ev_names, value = TRUE)] <- -0.25
A[1, grep("^2bk", ev_names, value = TRUE)] <-  0.25

cat("Contrast matrix:\n"); print(A)

# ** Run on held-out subjects (not in the 500-subject population) **
subjects <- c(501, 550)
results <- list()

for (s in subjects) {
  cat("\n--- Subject", s, "---\n")
  results[[as.character(s)]] <- fit_bestHRF(
    regularize_result = regularize_result,
    BOLD_file     = session_data$BOLD_files[s],
    EVs           = session_data$EVs_list[[s]],
    nuisance_file = session_data$nuisance_files[s],
    TR = 0.72,
    use = "population",
    contrasts = A,
    onsets  = TRUE,
    offsets = TRUE,
    verbose = 1
  )
}

# ** Common zlim across all subjects + both HRF modes (so plots are comparable) **
all_t <- unlist(lapply(results, function(r) {
  c(as.matrix(r$working$contrasts$tstat),
    as.matrix(r$population$contrasts$tstat))
}))
zmax <- max(abs(all_t), na.rm = TRUE)
zlim <- c(-zmax, zmax)
cat("\nCommon zlim:", zlim, "\n")

# ** Plot **
cat("Saving plots to", plot_dir, "\n")
for (s in subjects) {
  r <- results[[as.character(s)]]
  plot(r$working$contrasts$tstat,    idx = 1, zlim = zlim,
       fname = file.path(plot_dir, paste0("subject", s, "_working_2bk_minus_0bk_tstat.png")),
       title = paste0("Subject ", s, " - Working HRF - 2bk - 0bk (t-stat)"))
  plot(r$population$contrasts$tstat, idx = 1, zlim = zlim,
       fname = file.path(plot_dir, paste0("subject", s, "_population_2bk_minus_0bk_tstat.png")),
       title = paste0("Subject ", s, " - Population HRF - 2bk - 0bk (t-stat)"))
}

cat("Done.\n")
