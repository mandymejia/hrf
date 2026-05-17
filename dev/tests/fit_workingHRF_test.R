# Test: fit_workingHRF on the 4-subject MOTOR_LR fixture.
# Re-runs the fit and compares to the saved baseline. Used as a regression
# net before package-wide refactors (e.g., dependency cleanup).
#
# Run from repo root:
#   Rscript dev/tests/fit_workingHRF_test.R

library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)

# ** Inputs **
sd <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")

# ** Saved baseline **
expected <- readRDS("dev/fixtures/fit_workingHRF_result_motorlr_4s.rds")

cat("Running fit_workingHRF on 4-subject MOTOR_LR fixture...\n")
t0 <- Sys.time()
got <- fit_workingHRF(
  BOLD = sd$BOLD_files,
  EVs = sd$EVs_list,
  TR = 0.72,
  brainstructures = c("left", "right"),
  resamp_res = 10000,
  hpf = 0.01,
  nuisance = sd$nuisance_files,
  hrf_params = list(a1 = 6, b1 = 1, c = 1/6, a2 = 16, b2 = 1),
  derivatives = FALSE,
  onsets = TRUE, offsets = TRUE,
  scrub = NULL,
  alpha = 0.01,
  min_active_subjects = 20,
  verbose = 0,
  n_cores = 1,
  log_dir = "dev/logs"
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

# subject statuses
got_status <- sapply(got$subject_results, function(s) s$status)
exp_status <- sapply(expected$subject_results, function(s) s$status)
check("subject statuses identical", identical(got_status, exp_status),
      sprintf("got=[%s] expected=[%s]",
              paste(got_status, collapse=","),
              paste(exp_status, collapse=",")))

# activation masks
am_got <- got$activation_masks
am_exp <- expected$activation_masks
check("activation_masks names",
      identical(sort(names(am_got)), sort(names(am_exp))))
check("activation_masks$alpha",
      identical(am_got$alpha, am_exp$alpha))
check("activation_masks$n_subjects",
      identical(am_got$n_subjects, am_exp$n_subjects))
check("activation_masks$prop (within 1e-8)",
      length(am_got$prop) == length(am_exp$prop) &&
        all(abs(am_got$prop - am_exp$prop) < 1e-8 |
            (is.na(am_got$prop) & is.na(am_exp$prop))))
check("activation_masks$masks identical",
      identical(am_got$masks, am_exp$masks))

# hrf_params
check("hrf_params identical", identical(got$hrf_params, expected$hrf_params))

# session_info (everything except call timestamps)
check("session_info identical",
      identical(got$session_info, expected$session_info))

# call_info bits that should be stable
check("call_info$n_subjects",
      identical(got$call_info$n_subjects, expected$call_info$n_subjects))
check("call_info$alpha",
      identical(got$call_info$alpha, expected$call_info$alpha))
check("call_info$min_active_subjects",
      identical(got$call_info$min_active_subjects, expected$call_info$min_active_subjects))

cat("\n")
if (fail == 0) {
  cat("==== All checks passed ====\n")
} else {
  cat("==== ", fail, " check(s) FAILED ====\n", sep = "")
  quit(status = 1)
}
