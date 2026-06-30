# Test: fit_allHRFs on the 4-subject MOTOR_LR fixture.
#
# Asserts the combined-return contract:
#   class(combo)                          = "hrfs"
#   class(combo$fit_workingHRF)           = "workingHRF"
#   class(combo$fit_allHRFs)              = "allHRFs"
#   class(combo$regularize_allHRFs)       = "regularizeHRFs"
#   attr(combo$fit_allHRFs, "result_paths") = character vector (qs paths)
#
# Each sub-object numerically equivalent to the standalone-fn fixture:
#   combo$fit_workingHRF      vs  dev/fixtures/fit_workingHRF_result_motorlr_4s.rds
#   combo$fit_allHRFs         vs  dev/fixtures/fit_allHRFs_result_motorlr_4s.rds (best_per_c)
#   combo$regularize_allHRFs  vs  dev/fixtures/regularize_allHRFs_result_motorlr_4s.rds
#
# Run from repo root:
#   Rscript dev/tests/fit_allHRFs_test.R

library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)

sd <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")

expected_working <- readRDS("dev/fixtures/fit_workingHRF_result_motorlr_4s.rds")
expected_allHRF  <- readRDS("dev/fixtures/fit_allHRFs_result_motorlr_4s.rds")
expected_reg     <- readRDS("dev/fixtures/regularize_allHRFs_result_motorlr_4s.rds")

# Match the baseline fit_workingHRF call: derivatives = FALSE, alpha = 0.01.
cat("Running fit_allHRFs on 4-subject MOTOR_LR fixture...\n")
t0 <- Sys.time()
combo <- fit_allHRFs(
  BOLD = sd$BOLD_files, EVs = sd$EVs_list, TR = 0.72,
  brainstructures = c("left", "right"), resamp_res = 10000,
  hpf = 0.01, nuisance = sd$nuisance_files,
  onsets = TRUE, offsets = TRUE,
  scrub = NULL, smoothing = TRUE, surf_FWHM = 5,
  n_cores = 1, save_rss = FALSE, log_dir = "dev/logs", work_dir = NULL,
  working_hrf = list(a1 = 6, b1 = 1, c = 1/6, a2 = 16, b2 = 1),
  derivatives = FALSE,
  alpha = 0.01,
  min_active_subjects = 20,
  verbose = 0
)
cat("Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "sec\n\n")

fail <- 0
check <- function(name, ok, detail = "") {
  if (ok) cat("  [PASS]", name, "\n")
  else { cat("  [FAIL]", name, if (nzchar(detail)) paste0(" -- ", detail) else "", "\n"); fail <<- fail + 1 }
}

cat("== Top-level shape ==\n")
check("class(combo) == 'hrfs'", identical(class(combo), "hrfs"))
check("combo has fit_workingHRF",   "fit_workingHRF"   %in% names(combo))
check("combo has fit_allHRFs",      "fit_allHRFs"      %in% names(combo))
check("combo has regularize_allHRFs", "regularize_allHRFs" %in% names(combo))
check("top-level has no stray legacy fields",
      all(!c("workingHRF_results","best_params_results","subject_results") %in% names(combo)),
      sprintf("names = %s", paste(names(combo), collapse = ",")))

cat("\n== fit_workingHRF sub-object equivalence ==\n")
w_got <- combo$fit_workingHRF
w_exp <- expected_working
check("class(combo$fit_workingHRF) == 'workingHRF'", identical(class(w_got), "workingHRF"))
got_status <- sapply(w_got$subject_results, function(s) s$status)
exp_status <- sapply(w_exp$subject_results, function(s) s$status)
check("subject statuses identical", identical(got_status, exp_status))
am_got <- w_got$activation_masks
am_exp <- w_exp$activation_masks
check("activation_masks names",      identical(sort(names(am_got)), sort(names(am_exp))))
check("activation_masks$alpha",      identical(am_got$alpha, am_exp$alpha))
check("activation_masks$n_subjects", identical(am_got$n_subjects, am_exp$n_subjects))
check("activation_masks$prop (1e-8)",
      length(am_got$prop) == length(am_exp$prop) &&
        all(abs(am_got$prop - am_exp$prop) < 1e-8 | (is.na(am_got$prop) & is.na(am_exp$prop))))
check("activation_masks$masks identical",       identical(am_got$masks, am_exp$masks))
check("activation_masks$mask_prop_NA identical", identical(am_got$mask_prop_NA, am_exp$mask_prop_NA))
check("hrf_params identical", identical(w_got$hrf_params, w_exp$hrf_params))
shared <- intersect(names(combo$session_info), names(w_exp$session_info))
check("top-level session_info shared fields identical to baseline",
      identical(combo$session_info[shared], w_exp$session_info[shared]))
check("no session_info on fit_workingHRF sub-object", is.null(w_got$session_info))
check("no session_info on fit_allHRFs sub-object",   is.null(combo$fit_allHRFs$session_info))

cat("\n== fit_allHRFs sub-object equivalence ==\n")
a_got <- combo$fit_allHRFs
a_exp <- expected_allHRF
check("class(combo$fit_allHRFs) == 'allHRFs'", identical(class(a_got), "allHRFs"))
check("attr(combo$fit_allHRFs, 'result_paths') is character",
      is.character(attr(a_got, "result_paths")))
check("result_paths length matches subject count",
      length(attr(a_got, "result_paths")) == length(sd$BOLD_files))
check("hrf_grid identical to baseline", identical(a_got$hrf_grid, a_exp$hrf_grid))
check("subject_results length matches baseline",
      length(a_got$subject_results) == length(a_exp$subject_results))
for (i in seq_along(a_got$subject_results)) {
  g <- a_got$subject_results[[i]]; e <- a_exp$subject_results[[i]]
  for (hemi in c("cortexL", "cortexR")) {
    g_best <- g$glm_result$mGLM0s[[hemi]]$best_per_c
    e_best <- e$glm_result$mGLM0s[[hemi]]$best_per_c
    if (is.null(g_best) || is.null(e_best)) next
    for (c_key in intersect(names(g_best), names(e_best))) {
      g_df <- g_best[[c_key]]; e_df <- e_best[[c_key]]
      check(sprintf("subj %d %s c=%s a1 identical", i, hemi, c_key),  identical(g_df$a1, e_df$a1))
      check(sprintf("subj %d %s c=%s b1 identical", i, hemi, c_key),  identical(g_df$b1, e_df$b1))
      check(sprintf("subj %d %s c=%s rss (1e-6)", i, hemi, c_key),
            length(g_df$rss) == length(e_df$rss) &&
              all(abs(g_df$rss - e_df$rss) < 1e-6 | (is.na(g_df$rss) & is.na(e_df$rss))))
    }
  }
}

cat("\n== regularize_allHRFs sub-object equivalence (folded inside fit_allHRFs) ==\n")
r_got <- combo$regularize_allHRFs
r_exp <- expected_reg
check("class(combo$regularize_allHRFs) == 'regularizeHRFs'", identical(class(r_got), "regularizeHRFs"))
check("hrf_grid identical", identical(r_got$hrf_grid, r_exp$hrf_grid))
check("winning_c identical", identical(r_got$winning_c, r_exp$winning_c))
check("c_votes identical",   identical(r_got$c_votes, r_exp$c_votes))
check("pop_avg t2p_mean (1e-6)",
      length(r_got$pop_avg$t2p_mean) == length(r_exp$pop_avg$t2p_mean) &&
        all(abs(r_got$pop_avg$t2p_mean - r_exp$pop_avg$t2p_mean) < 1e-6 |
            (is.na(r_got$pop_avg$t2p_mean) & is.na(r_exp$pop_avg$t2p_mean))))
check("pop_avg fwhm_mean (1e-6)",
      length(r_got$pop_avg$fwhm_mean) == length(r_exp$pop_avg$fwhm_mean) &&
        all(abs(r_got$pop_avg$fwhm_mean - r_exp$pop_avg$fwhm_mean) < 1e-6 |
            (is.na(r_got$pop_avg$fwhm_mean) & is.na(r_exp$pop_avg$fwhm_mean))))
check("pop_avg a1_snapped identical", identical(r_got$pop_avg$a1_snapped, r_exp$pop_avg$a1_snapped))
check("pop_avg b1_snapped identical", identical(r_got$pop_avg$b1_snapped, r_exp$pop_avg$b1_snapped))
check("best_params_df nrow identical",
      nrow(r_got$best_params_df) == nrow(r_exp$best_params_df))
check("mask_prop_NA identical", identical(r_got$mask_prop_NA, r_exp$mask_prop_NA))

cat("\n")
if (fail == 0) {
  cat("==== All checks passed ====\n")
} else {
  cat("==== ", fail, " check(s) FAILED ====\n", sep = "")
  quit(status = 1)
}
