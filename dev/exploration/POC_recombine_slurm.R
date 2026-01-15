library(dplyr)

# Recombine chunked results from SLURM array jobs

output_dir <- "/N/project/hrf_adaptation/validation_bonferroni_0.05/cm_poc/chunks"
run_names <- c('LR', 'RL')
n_chunks <- 10

for(r in 1:2) {
  cat("\n********** Recombining", run_names[r], "**********\n")

  # Load all chunks
  all_chunks <- list()
  for (chunk_id in 1:n_chunks) {
    chunk_file <- file.path(output_dir, paste0("results_cm_poc_gambling_", run_names[r], "_chunk_", chunk_id, ".rds"))

    if (!file.exists(chunk_file)) {
      warning("Chunk file not found: ", chunk_file)
      next
    }

    cat("Loading chunk", chunk_id, "...\n")
    all_chunks[[chunk_id]] <- readRDS(chunk_file)
  }

  if (length(all_chunks) == 0) {
    warning("No chunks found for ", run_names[r])
    next
  }

  # Flatten list of lists and combine
  cat("Combining", length(all_chunks), "chunks...\n")
  results_all <- bind_rows(unlist(all_chunks, recursive = FALSE))

  cat("Total rows:", nrow(results_all), "\n")
  cat("Unique subjects:", length(unique(results_all$subject_id)), "\n")

  # Save combined results
  output_file <- file.path(output_dir, paste0("results_combined_cm_poc_gambling_", run_names[r], ".rds"))
  saveRDS(results_all, output_file)
  cat("*** Saved combined results to:", output_file, "***\n")
}

cat("\n*** Recombination complete! ***\n")