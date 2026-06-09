# Regenerates dev/fixtures/fit_bestHRF_result_motorlr_4s_pop.rds with the
# three-configuration baseline used by dev/tests/fit_bestHRF_test.R:
#   - adapted_run
#   - personalized_lookup_run
#   - personalized_refit_run
#
# Run from repo root:
#   Rscript dev/regenerate_fit_bestHRF_baseline.R

library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)

session_data <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")
reg_result   <- readRDS("dev/fixtures/regularize_allHRFs_result_motorlr_4s.rds")
allHRF_res   <- readRDS("dev/fixtures/fit_allHRFs_result_motorlr_4s.rds")

subj <- 1L

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
    verbose = 1
  )
}

cat("\n== adapted only ==\n")
adapted_run <- call_fit("adapted")

cat("\n== personalized via lookup ==\n")
personalized_lookup_run <- call_fit("personalized",
                                    allHRF_result = allHRF_res,
                                    subject_idx = subj)

cat("\n== personalized via refit ==\n")
personalized_refit_run <- call_fit("personalized")  # no allHRF_result = refit

baseline <- list(
  adapted_run             = adapted_run,
  personalized_lookup_run = personalized_lookup_run,
  personalized_refit_run  = personalized_refit_run
)

saveRDS(baseline, "dev/fixtures/fit_bestHRF_result_motorlr_4s_pop.rds")
cat("\nBaseline saved.\n")
cat("Personalized lookup winner:", personalized_lookup_run$personalized$winning_candidate_id, "\n")
cat("Personalized refit  winner:", personalized_refit_run$personalized$winning_candidate_id, "\n")
