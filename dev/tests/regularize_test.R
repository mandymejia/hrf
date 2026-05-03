# Test: regularize_allHRFs in refit mode (save_rss = FALSE).
# Re-runs regularize on the 4-subject MOTOR_LR fixture and compares to the
# saved baseline. Used as the regression test before/after parallelizing
# fit_candidate_maps_refit.
#
# Run from repo root:
#   Rscript dev/tests/regularize_test.R

library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)

# ** Inputs **
fwh <- readRDS("dev/fixtures/fit_workingHRF_result_motorlr_4s.rds")
fah <- readRDS("dev/fixtures/fit_allHRFs_result_motorlr_4s_norss.rds")
sd  <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")

# ** Saved baseline (sequential refit, recorded BEFORE the refit refactor) **
expected <- readRDS("dev/fixtures/regularize_allHRFs_result_motorlr_4s_norss.rds")

# ** Assertion harness **
fail <- 0
check <- function(name, ok, detail = "") {
  if (ok) {
    cat("  [PASS]", name, "\n")
  } else {
    cat("  [FAIL]", name, if (nzchar(detail)) paste0(" — ", detail) else "", "\n")
    fail <<- fail + 1
  }
}

compare_against_baseline <- function(got, label) {
  cat("\n--- ", label, " ---\n", sep = "")
  check("class", identical(class(got), class(expected)))
  check("names match", identical(sort(names(got)), sort(names(expected))))
  check("winning_c", identical(got$winning_c, expected$winning_c),
        sprintf("got=%s expected=%s", got$winning_c, expected$winning_c))
  check("c_votes", identical(got$c_votes, expected$c_votes))
  check("subject_results n_rows",
        nrow(got$subject_results) == nrow(expected$subject_results))
  check("winning_candidate_id per subject",
        identical(got$subject_results$winning_candidate_id,
                  expected$subject_results$winning_candidate_id),
        sprintf("got=[%s] expected=[%s]",
                paste(got$subject_results$winning_candidate_id, collapse=","),
                paste(expected$subject_results$winning_candidate_id, collapse=",")))
  check("winning_weighted_RSS (within 1e-8)",
        all(abs(got$subject_results$winning_weighted_RSS -
                expected$subject_results$winning_weighted_RSS) < 1e-8))
  check("pop_avg identical", isTRUE(all.equal(got$pop_avg, expected$pop_avg)))
  check("candidate_maps length",
        length(got$candidate_maps) == length(expected$candidate_maps))
  check("candidate_maps content identical",
        isTRUE(all.equal(got$candidate_maps, expected$candidate_maps)))
  check("mask_prop_NA identical",
        isTRUE(all.equal(got$mask_prop_NA, expected$mask_prop_NA)))
  check("hrf_grid identical", isTRUE(all.equal(got$hrf_grid, expected$hrf_grid)))
}

run_refit <- function(n_cores) {
  regularize_allHRFs(
    workingHRF_results = fwh,
    allHRF_results     = fah,
    BOLD     = sd$BOLD_files,
    EVs      = sd$EVs_list,
    nuisance = sd$nuisance_files,
    TR = 0.72,
    onsets = TRUE, offsets = TRUE,
    seffects = TRUE,
    n_cores = n_cores,
    log_dir = "dev/logs",
    verbose = 0
  )
}

# === Test 1: sequential (n_cores = 1) — must match baseline ===
cat("Running regularize_allHRFs in refit mode (n_cores = 1)...\n")
t0 <- Sys.time()
got_seq <- run_refit(n_cores = 1)
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n")
compare_against_baseline(got_seq, "Sequential vs baseline")

# === Test 2: parallel (n_cores = 4) — must also match baseline ===
cat("\nRunning regularize_allHRFs in refit mode (n_cores = 4)...\n")
t0 <- Sys.time()
got_par <- run_refit(n_cores = 4)
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n")
compare_against_baseline(got_par, "Parallel vs baseline")

cat("\n")
if (fail == 0) {
  cat("==== All checks passed ====\n")
} else {
  cat("==== ", fail, " check(s) FAILED ====\n", sep = "")
  quit(status = 1)
}
