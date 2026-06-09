# Test: fit_bestHRF after the candidate-fitting migration.
#
# New surface (per the refactor):
#   - fit_bestHRF takes regularize_result (which supplies pop_avg + hrf_grid)
#     plus an optional lookup-mode pair: allHRF_result + subject_idx. The qs
#     cache carries both RSS and canonical-HRF Fstat. When allHRF_result is
#     NULL, personalized mode refits.
#
# Exercises three configurations on the 4-subject MOTOR_LR fixture (subject 1)
# and compares each to a saved baseline:
#   1. adapted only
#   2. personalized via lookup (allHRF_result + subject_idx)
#   3. personalized via refit (allHRF_result = NULL)
#
# Run from repo root:
#   Rscript dev/tests/fit_bestHRF_test.R

library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)

# ** Inputs **
session_data <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")
reg_result   <- readRDS("dev/fixtures/regularize_allHRFs_result_motorlr_4s.rds")
allHRF_res   <- readRDS("dev/fixtures/fit_allHRFs_result_motorlr_4s.rds")

subj <- 1L

# ** Saved baseline **
expected <- readRDS("dev/fixtures/fit_bestHRF_result_motorlr_4s_pop.rds")

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
xifti_eq <- function(a, b, tol = 1e-8) {
  am <- as.matrix(a); bm <- as.matrix(b)
  identical(dim(am), dim(bm)) &&
    all(abs(am - bm) < tol | (is.na(am) & is.na(bm)))
}

compare_section <- function(label, got_sec, exp_sec) {
  check(paste0(label, "$betas"),
        xifti_eq(got_sec$betas, exp_sec$betas))
  check(paste0(label, "$contrasts$tstat"),
        xifti_eq(got_sec$contrasts$tstat, exp_sec$contrasts$tstat))
  check(paste0(label, "$contrasts$pval"),
        xifti_eq(got_sec$contrasts$pval, exp_sec$contrasts$pval))
  check(paste0(label, "$contrasts$pval_adj"),
        xifti_eq(got_sec$contrasts$pval_adj, exp_sec$contrasts$pval_adj))
  check(paste0(label, "$hrf_assignments"),
        isTRUE(all.equal(got_sec$hrf_assignments, exp_sec$hrf_assignments)))
}

call_fit <- function(use, allHRF_result = NULL, subject_idx = NULL) {
  fit_bestHRF(
    regularize_result = reg_result,
    BOLD_file = session_data$BOLD_files[subj],
    EVs       = session_data$EVs_list[[subj]],
    nuisance_file = session_data$nuisance_files[subj],
    TR = 0.72,
    use = use,
    allHRF_result = allHRF_result,
    subject_idx   = subject_idx,
    onsets = TRUE, offsets = TRUE,
    verbose = 0
  )
}

# === Run 1: adapted only ===
cat("\n== Run 1: use = 'adapted' ==\n")
t0 <- Sys.time()
got_adapted <- call_fit("adapted")
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "sec\n")
check("adapted run: class",          identical(class(got_adapted), "bestHRF"))
check("adapted run: df identical",   identical(got_adapted$df, expected$adapted_run$df))
check("adapted run: contrast_matrix identical",
      isTRUE(all.equal(got_adapted$contrast_matrix, expected$adapted_run$contrast_matrix)))
compare_section("adapted_run$working",  got_adapted$working,  expected$adapted_run$working)
compare_section("adapted_run$adapted",  got_adapted$adapted,  expected$adapted_run$adapted)

# === Run 2: personalized via lookup ===
cat("\n== Run 2: use = 'personalized' + lookup pair ==\n")
t0 <- Sys.time()
got_lookup <- call_fit("personalized",
                       allHRF_result = allHRF_res,
                       subject_idx = subj)
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "sec\n")
compare_section("personalized_lookup_run$working", got_lookup$working,
                expected$personalized_lookup_run$working)
compare_section("personalized_lookup_run$personalized", got_lookup$personalized,
                expected$personalized_lookup_run$personalized)
check("personalized_lookup_run: winning_candidate_id",
      identical(got_lookup$personalized$winning_candidate_id,
                expected$personalized_lookup_run$personalized$winning_candidate_id))
check("personalized_lookup_run: candidate_scores",
      isTRUE(all.equal(got_lookup$personalized$candidate_scores,
                       expected$personalized_lookup_run$personalized$candidate_scores)))

# === Run 3: personalized via refit ===
cat("\n== Run 3: use = 'personalized', allHRF_result = NULL (refit) ==\n")
t0 <- Sys.time()
got_refit <- call_fit("personalized")  # no lookup trio = refit
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "sec\n")
compare_section("personalized_refit_run$working", got_refit$working,
                expected$personalized_refit_run$working)
compare_section("personalized_refit_run$personalized", got_refit$personalized,
                expected$personalized_refit_run$personalized)
check("personalized_refit_run: winning_candidate_id",
      identical(got_refit$personalized$winning_candidate_id,
                expected$personalized_refit_run$personalized$winning_candidate_id))
check("personalized_refit_run: candidate_scores",
      isTRUE(all.equal(got_refit$personalized$candidate_scores,
                       expected$personalized_refit_run$personalized$candidate_scores)))

cat("\n")
if (fail == 0) {
  cat("==== All checks passed ====\n")
} else {
  cat("==== ", fail, " check(s) FAILED ====\n", sep = "")
  quit(status = 1)
}
