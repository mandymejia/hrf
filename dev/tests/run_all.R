# Run all dev regression tests and print a summary at the end.
# Each test runs in its own Rscript subprocess so a quit(status=1) in one
# doesn't abort the rest.
#
# Run from repo root:
#   Rscript dev/tests/run_all.R

tests <- c(
  "dev/tests/fit_workingHRF_test.R",
  "dev/tests/fit_allHRFs_test.R",
  "dev/tests/regularize_test.R",
  "dev/tests/fit_bestHRF_test.R"
)

results <- list()
for (t in tests) {
  cat(sprintf("Running %s ... ", basename(t)))
  utils::flush.console()
  t0 <- Sys.time()
  out <- suppressWarnings(system2("Rscript", c(t), stdout = TRUE, stderr = TRUE))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  status <- attr(out, "status")
  passed <- is.null(status) || status == 0

  pass_count <- length(grep("\\[PASS\\]", out))
  fail_count <- length(grep("\\[FAIL\\]", out))
  total <- pass_count + fail_count

  label <- if (passed && fail_count == 0) "PASS" else "FAIL"
  cat(sprintf("[%s] %d/%d checks  %.1fs\n", label, pass_count, total, elapsed))

  if (!passed || fail_count > 0) {
    cat("--- failed test output (last 30 lines) ---\n")
    cat(tail(out, 30), sep = "\n"); cat("\n")
  }

  results[[t]] <- list(
    passed = passed,
    pass_count = pass_count,
    fail_count = fail_count,
    total = total,
    elapsed = elapsed
  )
}

cat("\n\n==================== SUMMARY ====================\n")
total_pass <- 0; total_fail <- 0; total_time <- 0
for (t in names(results)) {
  r <- results[[t]]
  label <- if (r$passed && r$fail_count == 0) "PASS" else "FAIL"
  cat(sprintf("  [%s] %-40s %d/%d checks  %5.1fs\n",
              label, basename(t), r$pass_count, r$total, r$elapsed))
  total_pass <- total_pass + r$pass_count
  total_fail <- total_fail + r$fail_count
  total_time <- total_time + r$elapsed
}
cat("=================================================\n")
cat(sprintf("  TOTAL: %d/%d checks passed  (%.1fs)\n",
            total_pass, total_pass + total_fail, total_time))

if (total_fail > 0 || any(!sapply(results, `[[`, "passed"))) {
  cat("\nFailures detected.\n")
  quit(status = 1)
}
