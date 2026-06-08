# Test: candidate-fitting helpers ported into fit_bestHRF.R produce numerics
# identical to regularize_allHRFs's saved output.
#
# Locks (against dev/fixtures/regularize_allHRFs_result_motorlr_4s_norss.rds):
#   - build_candidate_maps(pop_avg, hrf_grid) -> $candidate_maps, $offset_combos
#   - pick_winning_candidate(scores) -> which.min behavior + ties
#
# These helpers live in BOTH fit_bestHRF.R and regularize_allHRFs.R during the
# migration window. The test only exercises the fit_bestHRF copies; regularize
# still uses its own copies (covered by regularize_test.R).
#
# Run from repo root:
#   Rscript dev/tests/fit_bestHRF_candidate_helpers_test.R

devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)

expected <- readRDS("dev/fixtures/regularize_allHRFs_result_motorlr_4s_norss.rds")

fail <- 0
check <- function(name, ok, detail = "") {
  if (ok) {
    cat("  [PASS]", name, "\n")
  } else {
    cat("  [FAIL]", name, if (nzchar(detail)) paste0(" - ", detail) else "", "\n")
    fail <<- fail + 1
  }
}

# === build_candidate_maps ===
cat("\n--- build_candidate_maps vs regularize candidate_maps ---\n")
got <- build_candidate_maps(
  pop_avg  = expected$pop_avg,
  hrf_grid = expected$hrf_grid,
  verbose  = 0
)
check("returns 25 candidates",     length(got$candidate_maps) == 25)
check("offset_combos identical",   isTRUE(all.equal(got$offset_combos, expected$offset_combos)))
check("candidate_maps identical",  isTRUE(all.equal(got$candidate_maps, expected$candidate_maps)))

# Per-candidate spot checks for clearer failure messages if mismatched
for (j in c(1, 7, 13, 19, 25)) {
  check(sprintf("candidate_maps[[%d]] identical", j),
        isTRUE(all.equal(got$candidate_maps[[j]], expected$candidate_maps[[j]])))
}

# === pick_winning_candidate ===
cat("\n--- pick_winning_candidate ---\n")
check("which.min basic",       pick_winning_candidate(c(10, 3, 7)) == 2L)
check("which.min first tie",   pick_winning_candidate(c(5, 5, 9)) == 1L)
check("matches regularize subject_results winning_candidate_id",
      all(vapply(seq_len(nrow(expected$subject_results)), function(i) {
        scores <- as.numeric(expected$subject_results[i, paste0("weighted_RSS_", seq_along(expected$candidate_maps))])
        pick_winning_candidate(scores) == expected$subject_results$winning_candidate_id[i]
      }, logical(1))))

cat("\n")
if (fail == 0) {
  cat("==== All checks passed ====\n")
} else {
  cat("==== ", fail, " check(s) FAILED ====\n", sep = "")
  quit(status = 1)
}
