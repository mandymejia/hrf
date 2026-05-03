# Test: fit_bestHRF in population mode.
# Re-runs fit_bestHRF on the 4-subject MOTOR_LR fixture (subject 1) and
# compares to the saved baseline.
#
# Run from repo root:
#   Rscript dev/tests/fit_bestHRF_test.R

library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)

# ** Inputs **
session_data <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")
reg_result   <- readRDS("dev/fixtures/regularize_result_motorlr_1125s_slim.rds")
reg_result$subject_results <- NULL  # force population-only regularize

# ** Saved baseline **
expected <- readRDS("dev/fixtures/fit_bestHRF_result_motorlr_4s_pop.rds")

cat("Running fit_bestHRF (population mode, subject 1)...\n")
t0 <- Sys.time()
got <- fit_bestHRF(
  reg_result,
  BOLD_file = session_data$BOLD_files[1],
  EVs       = session_data$EVs_list[[1]],
  nuisance_file = session_data$nuisance_files[1],
  TR = 0.72,
  use = "population",
  onsets = TRUE, offsets = TRUE,
  verbose = 0
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("Elapsed:", round(elapsed, 1), "sec\n\n")

# ** Assertions **
fail <- 0
check <- function(name, ok, detail = "") {
  if (ok) {
    cat("  [PASS]", name, "\n")
  } else {
    cat("  [FAIL]", name, if (nzchar(detail)) paste0(" — ", detail) else "", "\n")
    fail <<- fail + 1
  }
}

# xifti compares via the underlying matrix
xifti_eq <- function(a, b, tol = 1e-8) {
  a_mat <- as.matrix(a); b_mat <- as.matrix(b)
  identical(dim(a_mat), dim(b_mat)) &&
    all(abs(a_mat - b_mat) < tol | (is.na(a_mat) & is.na(b_mat)))
}

cat("Comparing against baseline:\n")

check("class", identical(class(got), class(expected)))
check("names match", identical(sort(names(got)), sort(names(expected))))
check("df identical", identical(got$df, expected$df))
check("contrast_matrix identical", isTRUE(all.equal(got$contrast_matrix, expected$contrast_matrix)))

check("working$betas",            xifti_eq(got$working$betas,            expected$working$betas))
check("working$contrasts$tstat",  xifti_eq(got$working$contrasts$tstat,  expected$working$contrasts$tstat))
check("working$contrasts$pval",   xifti_eq(got$working$contrasts$pval,   expected$working$contrasts$pval))
check("working$contrasts$pval_adj", xifti_eq(got$working$contrasts$pval_adj, expected$working$contrasts$pval_adj))
check("working$hrf_assignments",  isTRUE(all.equal(got$working$hrf_assignments, expected$working$hrf_assignments)))

check("population$betas",         xifti_eq(got$population$betas,         expected$population$betas))
check("population$contrasts$tstat", xifti_eq(got$population$contrasts$tstat, expected$population$contrasts$tstat))
check("population$contrasts$pval", xifti_eq(got$population$contrasts$pval, expected$population$contrasts$pval))
check("population$contrasts$pval_adj", xifti_eq(got$population$contrasts$pval_adj, expected$population$contrasts$pval_adj))
check("population$hrf_assignments", isTRUE(all.equal(got$population$hrf_assignments, expected$population$hrf_assignments)))

cat("\n")
if (fail == 0) {
  cat("==== All checks passed ====\n")
} else {
  cat("==== ", fail, " check(s) FAILED ====\n", sep = "")
  quit(status = 1)
}
