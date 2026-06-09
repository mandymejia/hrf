# Test: candidate-fitting helpers ported into fit_bestHRF.R produce numerics
# identical to regularize_allHRFs's saved output.
#
# Locks:
#   - build_candidate_maps(pop_avg, hrf_grid) -> $candidate_maps, $offset_combos
#       vs regularize_allHRFs_result_motorlr_4s_norss.rds (refit-mode fixture;
#       candidate_maps geometry doesn't depend on save_rss).
#   - pick_winning_candidate(scores) -> which.min + saved winners (refit
#       fixture's subject_results).
#   - score_candidates_lookup(candidate_maps, RSS, fstat) -> matches
#       regularize_allHRFs_result_motorlr_4s.rds (save_rss=TRUE / lookup mode)
#       subject_results$weighted_RSS_1..25 columns. Exercised under FOUR
#       conditions: {qs vs in-memory RSS} x {sequential vs parallel}.
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

# === score_candidates_lookup ===
# Pull the lookup-mode regularize fixture, plus the inputs it was built from.
cat("\n--- score_candidates_lookup vs lookup-mode regularize ---\n")
expected_lookup <- readRDS("dev/fixtures/regularize_allHRFs_result_motorlr_4s.rds")
allHRF_res      <- readRDS("dev/fixtures/fit_allHRFs_result_motorlr_4s.rds")
working_res     <- readRDS("dev/fixtures/fit_workingHRF_result_motorlr_4s.rds")
cands           <- expected_lookup$candidate_maps
n_subjects      <- nrow(expected_lookup$subject_results)
expected_scores <- as.matrix(expected_lookup$subject_results[, paste0("weighted_RSS_", seq_along(cands))])

# Helper: load one subject's RSS + Fstat (qs read + workingHRF Fstat).
load_subject_data <- function(i) {
  qs_path <- attr(allHRF_res, "result_paths")[i]
  qs_obj  <- qs2::qs_read(qs_path)
  RSS <- rbind(qs_obj$glm_result$mGLM0s$cortexL$RSS,
               qs_obj$glm_result$mGLM0s$cortexR$RSS)
  Fstat <- c(working_res$subject_results[[i]]$glm_results$mGLM0s$cortexL$Fstat,
             working_res$subject_results[[i]]$glm_results$mGLM0s$cortexR$Fstat)
  list(RSS = RSS, Fstat = Fstat)
}

scores_to_mat <- function(scores_list) {
  do.call(rbind, scores_list)
}

# Condition 1: with qs, sequential -- each iteration qs_reads.
got_qs_seq <- scores_to_mat(lapply(seq_len(n_subjects), function(i) {
  d <- load_subject_data(i)
  score_candidates_lookup(cands, d$RSS, d$Fstat)
}))
check("qs + sequential matches regularize lookup",
      isTRUE(all.equal(got_qs_seq, expected_scores, check.attributes = FALSE)))

# Condition 2: without qs, sequential -- preload RSS+Fstat once, score in loop.
preloaded <- lapply(seq_len(n_subjects), load_subject_data)
got_mem_seq <- scores_to_mat(lapply(seq_len(n_subjects), function(i) {
  score_candidates_lookup(cands, preloaded[[i]]$RSS, preloaded[[i]]$Fstat)
}))
check("in-memory + sequential matches regularize lookup",
      isTRUE(all.equal(got_mem_seq, expected_scores, check.attributes = FALSE)))

# Condition 3: with qs, parallel -- each worker qs_reads its own subject.
cl <- parallel::makeCluster(2)
parallel::clusterEvalQ(cl, {
  devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)
  library(qs2)
})
parallel::clusterExport(cl, c("allHRF_res", "working_res", "cands",
                              "load_subject_data"))
got_qs_par <- scores_to_mat(parallel::parLapply(cl, seq_len(n_subjects), function(i) {
  d <- load_subject_data(i)
  score_candidates_lookup(cands, d$RSS, d$Fstat)
}))
parallel::stopCluster(cl)
check("qs + parallel matches regularize lookup",
      isTRUE(all.equal(got_qs_par, expected_scores, check.attributes = FALSE)))

# Condition 4: without qs, parallel -- preloaded RSS shipped via clusterExport.
cl <- parallel::makeCluster(2)
parallel::clusterEvalQ(cl, devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE))
parallel::clusterExport(cl, c("preloaded", "cands"))
got_mem_par <- scores_to_mat(parallel::parLapply(cl, seq_len(n_subjects), function(i) {
  score_candidates_lookup(cands, preloaded[[i]]$RSS, preloaded[[i]]$Fstat)
}))
parallel::stopCluster(cl)
check("in-memory + parallel matches regularize lookup",
      isTRUE(all.equal(got_mem_par, expected_scores, check.attributes = FALSE)))

# Sanity: pick_winning_candidate over the lookup scores reproduces the
# lookup fixture's winning_candidate_id.
got_winners <- apply(got_qs_seq, 1, pick_winning_candidate)
check("pick_winning_candidate on lookup scores matches saved winners",
      identical(as.integer(got_winners),
                as.integer(expected_lookup$subject_results$winning_candidate_id)))

# === score_candidates_refit ===
# Refit-mode regularize fixture (no qs RSS, scores come from a fresh multiGLM
# refit of all unique candidate designs). Lock single-subject helper output
# against the per-subject rows of that fixture, sequentially and in parallel.
cat("\n--- score_candidates_refit vs refit-mode regularize ---\n")
library(ciftiTools)
ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
sd          <- readRDS("dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds")
cands_refit <- expected$candidate_maps  # same offset grid; refit fixture
expected_refit_scores <- as.matrix(
  expected$subject_results[, paste0("weighted_RSS_", seq_along(cands_refit))]
)

# Helper: load subject i's BOLD_xii + nuisance (caller's responsibility).
load_refit_inputs <- function(i) {
  BOLD_xii <- ciftiTools::read_cifti(
    sd$BOLD_files[[i]], brainstructures = c("left", "right"), resamp_res = 10000
  )
  nuisance <- as.matrix(utils::read.table(sd$nuisance_files[[i]], header = FALSE))
  list(BOLD_xii = BOLD_xii, nuisance = nuisance, EVs = sd$EVs_list[[i]])
}

# Sub-condition A: subject 1 sequential -- proves single-call equivalence.
cat("  Running subject 1 sequentially (~12s)...\n")
t0 <- Sys.time()
inp1 <- load_refit_inputs(1)
got_refit_s1 <- score_candidates_refit(
  candidate_maps = cands_refit, hrf_grid = expected$hrf_grid,
  BOLD_xii = inp1$BOLD_xii, EVs = inp1$EVs, nuisance = inp1$nuisance,
  TR = 0.72, onsets = TRUE, offsets = TRUE
)
cat("  Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "sec\n")
check("subject 1 refit scores match regularize refit fixture",
      isTRUE(all.equal(got_refit_s1,
                       as.numeric(expected_refit_scores[1, ]),
                       check.attributes = FALSE)))

# Sub-condition B: all 4 subjects in parallel -- proves it survives parLapplyLB.
cat("  Running 4 subjects in parallel via PSOCK (~22s)...\n")
t0 <- Sys.time()
cl <- parallel::makeCluster(2)
parallel::clusterEvalQ(cl, {
  devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)
  library(ciftiTools)
  ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
})
parallel::clusterExport(cl, c("sd", "cands_refit", "expected", "load_refit_inputs"))
got_refit_par <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(n_subjects), function(i) {
  inp <- load_refit_inputs(i)
  score_candidates_refit(
    candidate_maps = cands_refit, hrf_grid = expected$hrf_grid,
    BOLD_xii = inp$BOLD_xii, EVs = inp$EVs, nuisance = inp$nuisance,
    TR = 0.72, onsets = TRUE, offsets = TRUE
  )
}))
parallel::stopCluster(cl)
cat("  Elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "sec\n")
check("all 4 subjects refit scores match regularize refit fixture (parallel)",
      isTRUE(all.equal(got_refit_par, expected_refit_scores, check.attributes = FALSE)))

# Sanity: pick_winning_candidate on refit scores reproduces saved winners.
refit_winners <- apply(got_refit_par, 1, pick_winning_candidate)
check("pick_winning_candidate on refit scores matches saved winners",
      identical(as.integer(refit_winners),
                as.integer(expected$subject_results$winning_candidate_id)))

cat("\n")
if (fail == 0) {
  cat("==== All checks passed ====\n")
} else {
  cat("==== ", fail, " check(s) FAILED ====\n", sep = "")
  quit(status = 1)
}
