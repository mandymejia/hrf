# Test: fit_bestHRF with n_cores > 1 (PSOCK cluster, mode-level parallelism).
# 1) Verifies parallel run produces numerically identical results vs n_cores = 1.
# 2) Verifies log_file captures per-mode timestamp lines from both workers.
#
# Run from repo root:
#   Rscript dev/tests/fit_bestHRF_parallel_test.R

library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)

session_data <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")
reg_result   <- readRDS("dev/fixtures/regularize_result_motorlr_1125s_slim.rds")
reg_result$subject_results <- NULL  # adapted-only

log_path <- tempfile(pattern = "fit_bestHRF_parlog_", fileext = ".txt")
cat("Log file:", log_path, "\n\n")

cat("== Run 1: n_cores=2, log_file set ==\n")
t0 <- Sys.time()
got_par <- fit_bestHRF(
  reg_result,
  BOLD_file = session_data$BOLD_files[1],
  EVs       = session_data$EVs_list[[1]],
  nuisance_file = session_data$nuisance_files[1],
  TR = 0.72,
  use = "adapted",
  onsets = TRUE, offsets = TRUE,
  log_file = log_path,
  n_cores = 2,
  verbose = 1
)
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n\n")

cat("== Run 2: n_cores=1 (sequential reference) ==\n")
t0 <- Sys.time()
got_seq <- fit_bestHRF(
  reg_result,
  BOLD_file = session_data$BOLD_files[1],
  EVs       = session_data$EVs_list[[1]],
  nuisance_file = session_data$nuisance_files[1],
  TR = 0.72,
  use = "adapted",
  onsets = TRUE, offsets = TRUE,
  n_cores = 1,
  verbose = 0
)
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n\n")

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

cat("Parallel vs sequential numerical equivalence:\n")
check("working$betas",   xifti_eq(got_par$working$betas,   got_seq$working$betas))
check("working$tstat",   xifti_eq(got_par$working$contrasts$tstat, got_seq$working$contrasts$tstat))
check("working$pval",    xifti_eq(got_par$working$contrasts$pval,  got_seq$working$contrasts$pval))
check("adapted$betas",   xifti_eq(got_par$adapted$betas,   got_seq$adapted$betas))
check("adapted$tstat",   xifti_eq(got_par$adapted$contrasts$tstat, got_seq$adapted$contrasts$tstat))
check("adapted$pval",    xifti_eq(got_par$adapted$contrasts$pval,  got_seq$adapted$contrasts$pval))
check("df equal",        identical(got_par$df, got_seq$df))

cat("\nLog file checks:\n")
if (!file.exists(log_path)) {
  check("log file exists", FALSE)
} else {
  log_lines <- readLines(log_path)
  ts_lines  <- grep("\\] \\[(working|adapted)\\] \\[pid [0-9]+\\]", log_lines, value = TRUE)
  pids      <- unique(sub(".*\\[pid ([0-9]+)\\].*", "\\1", ts_lines))
  labels    <- unique(sub(".*\\] \\[(working|adapted)\\] \\[.*", "\\1", ts_lines))
  check("log file exists",                        TRUE)
  check("contains per-mode timestamp lines",      length(ts_lines) >= 4,
        paste("got", length(ts_lines), "lines"))
  check("both mode labels present",               all(c("working", "adapted") %in% labels))
  check(">=2 distinct pids (parallel workers)",   length(pids) >= 2,
        paste("pids =", paste(pids, collapse = ",")))
  cat("\n  First timestamp lines:\n")
  for (l in head(ts_lines, 6)) cat("   ", l, "\n")
}

cat("\n")
if (fail == 0) {
  cat("==== All parallel checks passed ====\n")
} else {
  cat("==== ", fail, " check(s) FAILED ====\n", sep = "")
  quit(status = 1)
}
