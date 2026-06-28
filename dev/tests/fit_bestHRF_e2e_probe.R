# Probe: end-to-end pipeline + edge cases for the post-migration hrf-z package.
# Not a regression test against fixtures -- this is an integration-level smoke
# probe to catch issues the unit/regression suites might miss.
#
# Covers:
#   1. Full pipeline run on subject 1 from raw inputs:
#        regularize (fixtures) -> fit_bestHRF (refit + lookup)
#   2. use = c("personalized", "adapted") -- both modes in one call.
#   3. n_cores = 3 with both modes -- 3-way mode parallelism + personalized.
#   4. Validation errors: invalid subject_idx, missing qs file.
#   5. Lookup vs refit numerically equivalent on training subject (a meaningful
#      sanity check: both paths should pick the same winner with similar scores
#      since they operate on the same smoothed-unscaled BOLD).
#
# Run from repo root:
#   Rscript dev/tests/fit_bestHRF_e2e_probe.R

library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)

session_data <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")
reg_result   <- readRDS("dev/fixtures/regularize_allHRFs_result_motorlr_4s.rds")
allHRF_res   <- readRDS("dev/fixtures/fit_allHRFs_result_motorlr_4s.rds")
subj <- 1L

fail <- 0
check <- function(name, ok, detail = "") {
  if (ok) {
    cat("  [PASS]", name, "\n")
  } else {
    cat("  [FAIL]", name, if (nzchar(detail)) paste0(" - ", detail) else "", "\n")
    fail <<- fail + 1
  }
}
xifti_eq <- function(a, b, tol = 1e-8) {
  am <- as.matrix(a); bm <- as.matrix(b)
  identical(dim(am), dim(bm)) &&
    all(abs(am - bm) < tol | (is.na(am) & is.na(bm)))
}

common <- list(
  regularize_result = reg_result,
  BOLD_file = session_data$BOLD_files[subj],
  EVs       = session_data$EVs_list[[subj]],
  nuisance_file = session_data$nuisance_files[subj],
  TR = 0.72,
  onsets = TRUE, offsets = TRUE,
  verbose = 0
)

# === 1. use = c("personalized", "adapted") both modes ===
cat("\n== 1. Both modes in one call (lookup) ==\n")
t0 <- Sys.time()
got_both <- do.call(fit_bestHRF, c(common, list(
  use = c("personalized", "adapted"),
  allHRF_result = allHRF_res, subject_idx = subj
)))
cat("  Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n")
check("result has $working",      !is.null(got_both$working))
check("result has $adapted",      !is.null(got_both$adapted))
check("result has $personalized", !is.null(got_both$personalized))
check("personalized has winning_candidate_id", !is.null(got_both$personalized$winning_candidate_id))
check("personalized has candidate_scores",     !is.null(got_both$personalized$candidate_scores))
check("adapted has no winning_candidate_id",   is.null(got_both$adapted$winning_candidate_id))
check("class bestHRF",                          identical(class(got_both), "bestHRF"))

# === 2. n_cores = 3 + both modes (parallel + personalized in PSOCK worker) ===
cat("\n== 2. n_cores = 3 + both modes ==\n")
t0 <- Sys.time()
got_par <- do.call(fit_bestHRF, c(common, list(
  use = c("personalized", "adapted"),
  allHRF_result = allHRF_res, subject_idx = subj,
  n_cores = 3
)))
cat("  Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n")
check("parallel betas match sequential (working)",
      xifti_eq(got_par$working$betas, got_both$working$betas))
check("parallel betas match sequential (adapted)",
      xifti_eq(got_par$adapted$betas, got_both$adapted$betas))
check("parallel betas match sequential (personalized)",
      xifti_eq(got_par$personalized$betas, got_both$personalized$betas))
check("personalized winner identical par vs seq",
      identical(got_par$personalized$winning_candidate_id,
                got_both$personalized$winning_candidate_id))

# === 3. Lookup vs refit consistency on training subject ===
cat("\n== 3. Lookup vs refit consistency (training subject) ==\n")
got_refit <- do.call(fit_bestHRF, c(common, list(use = "personalized")))
check("lookup winner == refit winner",
      identical(got_both$personalized$winning_candidate_id,
                got_refit$personalized$winning_candidate_id),
      sprintf("lookup=%d refit=%d",
              got_both$personalized$winning_candidate_id,
              got_refit$personalized$winning_candidate_id))
score_corr <- cor(got_both$personalized$candidate_scores,
                  got_refit$personalized$candidate_scores)
check("lookup vs refit candidate score correlation > 0.99",
      score_corr > 0.99,
      sprintf("cor = %.4f", score_corr))

# === 4. Validation errors ===
cat("\n== 4. Validation: bad inputs error cleanly ==\n")

err_bad_subj_idx <- tryCatch({
  do.call(fit_bestHRF, c(common, list(
    use = "personalized", allHRF_result = allHRF_res, subject_idx = 999
  )))
  "no error"
}, error = function(e) e$message)
check("subject_idx out of range -> error",
      grepl("not found|out of range|empty", err_bad_subj_idx, ignore.case = TRUE),
      err_bad_subj_idx)

err_missing_qs <- tryCatch({
  fake_allHRF <- allHRF_res
  attr(fake_allHRF, "result_paths")[subj] <- "/tmp/nonexistent.qs"
  do.call(fit_bestHRF, c(common, list(
    use = "personalized", allHRF_result = fake_allHRF, subject_idx = subj
  )))
  "no error"
}, error = function(e) e$message)
check("qs file missing -> error",
      grepl("not found|does not exist", err_missing_qs, ignore.case = TRUE),
      err_missing_qs)

err_bad_pop_avg <- tryCatch({
  bad_reg <- reg_result
  bad_reg$pop_avg$a1_snapped <- NULL
  do.call(fit_bestHRF, c(common[setdiff(names(common), "regularize_result")], list(
    regularize_result = bad_reg, use = "adapted"
  )))
  "no error"
}, error = function(e) e$message)
check("pop_avg missing columns -> error",
      grepl("missing column", err_bad_pop_avg, ignore.case = TRUE),
      err_bad_pop_avg)

# === 5. regularize on tiny data: lean output sanity ===
cat("\n== 5. regularize lean output sanity ==\n")
check("regularize result has no candidate_maps", !"candidate_maps" %in% names(reg_result))
check("regularize result has no subject_results", !"subject_results" %in% names(reg_result))
check("regularize result has xii attr",          !is.null(attr(reg_result, "xii")))
check("pop_avg c_snapped is uniform = winning_c",
      length(unique(reg_result$pop_avg$c_snapped)) == 1L &&
        abs(unique(reg_result$pop_avg$c_snapped) - reg_result$winning_c) < 1e-12)

cat("\n")
if (fail == 0) {
  cat("==== All probes passed ====\n")
} else {
  cat("==== ", fail, " probe(s) FAILED ====\n", sep = "")
  quit(status = 1)
}
