# Test: fit_allHRFs on the 4-subject MOTOR_LR fixture (save_rss = TRUE).
# Re-runs the fit and compares to the saved baseline. Used as a regression
# net before package-wide refactors (e.g., dependency cleanup).
#
# NOTE: fit_allHRFs writes per-subject .qs files to a fresh
# "work_<timestamp>_<rand>" subdir each run, so the result_paths attr will
# differ across runs. That's expected — we don't assert on it.
#
# Run from repo root:
#   Rscript dev/tests/fit_allHRFs_test.R

library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)

# ** Inputs **
sd <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")
my_hrf_grid <- generate_hrf_grid()

# ** Saved baseline **
expected <- readRDS("dev/fixtures/fit_allHRFs_result_motorlr_4s.rds")

cat("Running fit_allHRFs on 4-subject MOTOR_LR fixture (save_rss = TRUE)...\n")
t0 <- Sys.time()
got <- fit_allHRFs(
  BOLD = sd$BOLD_files,
  EVs = sd$EVs_list,
  TR = 0.72,
  hrf_grid = my_hrf_grid,
  brainstructures = c("left", "right"),
  resamp_res = 10000,
  hpf = 0.01,
  nuisance = sd$nuisance_files,
  onsets = TRUE, offsets = TRUE,
  scrub = NULL,
  save_rss = TRUE,
  verbose = 0,
  n_cores = 1,
  log_dir = "dev/logs",
  work_dir = "dev/work/"
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

cat("Comparing against baseline:\n")

check("class", identical(class(got), class(expected)))
check("names match", identical(sort(names(got)), sort(names(expected))))

# subject_results count
check("subject_results length",
      length(got$subject_results) == length(expected$subject_results))

# best_params_results data.frame — main scientific output
bp_got <- got$best_params_results
bp_exp <- expected$best_params_results
check("best_params_results dims identical",
      identical(dim(bp_got), dim(bp_exp)),
      sprintf("got=%s expected=%s",
              paste(dim(bp_got), collapse="x"),
              paste(dim(bp_exp), collapse="x")))
check("best_params_results column names identical",
      identical(sort(colnames(bp_got)), sort(colnames(bp_exp))))
# numeric columns within tolerance
for (col in intersect(colnames(bp_got), colnames(bp_exp))) {
  if (is.numeric(bp_got[[col]])) {
    check(paste0("best_params_results$", col, " (within 1e-8)"),
          length(bp_got[[col]]) == length(bp_exp[[col]]) &&
            all(abs(bp_got[[col]] - bp_exp[[col]]) < 1e-8 |
                (is.na(bp_got[[col]]) & is.na(bp_exp[[col]]))))
  } else {
    check(paste0("best_params_results$", col, " identical"),
          identical(bp_got[[col]], bp_exp[[col]]))
  }
}

# hrf_grid
check("hrf_grid identical", isTRUE(all.equal(got$hrf_grid, expected$hrf_grid)))

# session_info
check("session_info identical",
      identical(got$session_info, expected$session_info))

# call_info bits that are stable across runs
check("call_info$n_subjects",
      identical(got$call_info$n_subjects, expected$call_info$n_subjects))
check("call_info$n_hrf_models",
      identical(got$call_info$n_hrf_models, expected$call_info$n_hrf_models))
check("call_info$save_rss",
      identical(got$call_info$save_rss, expected$call_info$save_rss))

cat("\n")
if (fail == 0) {
  cat("==== All checks passed ====\n")
} else {
  cat("==== ", fail, " check(s) FAILED ====\n", sep = "")
  quit(status = 1)
}
