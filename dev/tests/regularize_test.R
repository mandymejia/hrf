# Test: regularize_allHRFs after Phase 2 strip of candidate fitting.
#
# regularize is now aggregation only: workingHRF + allHRFs -> pop_avg,
# best_params_df, winning_c, c_votes, hrf_grid (with t2p/FWHM), mask_prop_NA,
# plus xifti template as an attribute. Candidate scoring + per-subject
# winner selection moved to fit_bestHRF.
#
# Compares a fresh run on the 4-subject MOTOR_LR fixtures to the saved
# baseline.
#
# Run from repo root:
#   Rscript dev/tests/regularize_test.R

library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)

# ** Inputs **
fwh <- readRDS("dev/fixtures/fit_workingHRF_result_motorlr_4s.rds")
fah <- readRDS("dev/fixtures/fit_allHRFs_result_motorlr_4s.rds")

# ** Saved baseline (post-Phase 2 lean output) **
expected <- readRDS("dev/fixtures/regularize_allHRFs_result_motorlr_4s.rds")

# ** Assertion harness **
fail <- 0
check <- function(name, ok, detail = "") {
  if (ok) {
    cat("  [PASS]", name, "\n")
  } else {
    cat("  [FAIL]", name, if (nzchar(detail)) paste0(" - ", detail) else "", "\n")
    fail <<- fail + 1
  }
}

cat("Running regularize_allHRFs...\n")
t0 <- Sys.time()
got <- regularize_allHRFs(
  workingHRF_results = fwh,
  allHRF_results     = fah,
  verbose = 0
)
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n\n")

cat("--- Structural checks ---\n")
check("class regularizeHRFs", identical(class(got), class(expected)))
check("names match",          identical(sort(names(got)), sort(names(expected))))
check("no candidate_maps in result",  !"candidate_maps"  %in% names(got))
check("no subject_results in result", !"subject_results" %in% names(got))
check("no offset_combos in result",   !"offset_combos"   %in% names(got))
check("xii attribute attached",       !is.null(attr(got, "xii")))

cat("\n--- Aggregate outputs vs baseline ---\n")
check("winning_c identical", identical(got$winning_c, expected$winning_c),
      sprintf("got=%s expected=%s", got$winning_c, expected$winning_c))
check("c_votes identical",   identical(got$c_votes,   expected$c_votes))
check("pop_avg identical (cols + numeric)",
      isTRUE(all.equal(got$pop_avg, expected$pop_avg)))
check("best_params_df identical",
      isTRUE(all.equal(got$best_params_df, expected$best_params_df)))
check("hrf_grid identical (with t2p/FWHM)",
      isTRUE(all.equal(got$hrf_grid, expected$hrf_grid)))
check("mask_prop_NA identical",
      isTRUE(all.equal(got$mask_prop_NA, expected$mask_prop_NA)))

cat("\n--- pop_avg shape checks ---\n")
required_cols <- c("voxel", "t2p_mean", "fwhm_mean",
                   "a1_snapped", "b1_snapped", "c_snapped")
check("pop_avg has expected columns",
      all(required_cols %in% names(got$pop_avg)),
      paste("missing:", paste(setdiff(required_cols, names(got$pop_avg)), collapse=",")))
check("pop_avg$c_snapped is single value (== winning_c)",
      length(unique(got$pop_avg$c_snapped)) == 1 &&
        abs(unique(got$pop_avg$c_snapped) - got$winning_c) < 1e-12)

cat("\n")
if (fail == 0) {
  cat("==== All checks passed ====\n")
} else {
  cat("==== ", fail, " check(s) FAILED ====\n", sep = "")
  quit(status = 1)
}
