# Test: regularize_allHRFs (population-RSS aggregation).
#
# Two layers:
#   1. TOY tests — small synthetic (subject_results, subj_masks, hrf_grid,
#      pop_mask) with hand-calculated expected winners. Fast, deterministic.
#   2. REAL data test — if a 4-subject fixture with RSS retained is present,
#      run the full function end-to-end. Skipped otherwise (fixture wasn't
#      regenerated with save_rss).
#
# Run from repo root:
#   Rscript dev/tests/regularize_pop_rss_test.R

suppressPackageStartupMessages({
  devtools::load_all("~/Documents/Github/hrf-z", quiet = TRUE)
})

fail <- 0
check <- function(name, ok, detail = "") {
  if (isTRUE(ok)) {
    cat("  [PASS]", name, "\n")
  } else {
    cat("  [FAIL]", name, if (nzchar(detail)) paste0(" - ", detail) else "", "\n")
    fail <<- fail + 1
  }
}

# ---- Helpers ----------------------------------------------------------------

# Build a minimal fit_allHRFs subject_results entry with a given RSS matrix.
# RSS is split evenly across L and R hemispheres (rows) so rbind() reconstructs
# the full matrix inside sum_rss_across_subjects.
mk_subj <- function(RSS) {
  n_L <- floor(nrow(RSS) / 2)
  list(glm_result = list(mGLM0s = list(
    cortexL = list(RSS = RSS[seq_len(n_L), , drop = FALSE]),
    cortexR = list(RSS = RSS[(n_L + 1):nrow(RSS), , drop = FALSE])
  )))
}

# Small hrf_grid with 4 candidates on 2 c-values.
toy_grid <- data.frame(
  a1 = c(5, 6, 7, 8),
  b1 = c(1, 1, 1, 1),
  c  = c(0, 0, 1/6, 1/6)
)

# ---- TOY 1: hand-crafted winners ---------------------------------------------
cat("\n--- TOY 1: two subjects, three voxels, four candidates ---\n")

# n_voxels = 4 (so we can split L/R = 2/2 for rbind), n_candidates = 4.
# Subject 1 activated at voxels 1, 2. Subject 2 activated at voxels 2, 3.
# Voxel 4 is not activated by anyone.
# RSS is designed so:
#   voxel 1: subj1 only -> argmin at k=2
#   voxel 2: sum(subj1, subj2) at each k -> k=3 wins
#   voxel 3: subj2 only -> argmin at k=4
#   voxel 4: nobody contributes -> NA before modal unmask

# subj1 RSS: voxels x candidates
subj1_RSS <- rbind(
  c(10, 1, 20, 30),   # voxel 1: k=2 min
  c(5,  4,  3,  2),   # voxel 2: k=4 min alone, but k=3 wins after sum with subj2
  c(9,  9,  9,  9),   # voxel 3: not activated by subj1 (mask says no)
  c(9,  9,  9,  9)    # voxel 4: not activated
)
subj2_RSS <- rbind(
  c(9,  9,  9,  9),   # voxel 1: not activated by subj2
  c(6,  6,  1,  5),   # voxel 2: k=3 min alone; sum with subj1 = (11,10,4,7) -> k=3
  c(20, 30, 40, 2),   # voxel 3: k=4 min
  c(9,  9,  9,  9)    # voxel 4: not activated
)

subject_results <- list(mk_subj(subj1_RSS), mk_subj(subj2_RSS))

# subj_masks: voxels x subjects (TRUE = activated)
subj_masks <- cbind(
  c(TRUE,  TRUE,  FALSE, FALSE),   # subj 1
  c(FALSE, TRUE,  TRUE,  FALSE)    # subj 2
)

agg <- sum_rss_across_subjects(subject_results, subj_masks, n_candidates = 4, verbose = 0)

# Verify sums for voxel 2 (both contributed)
check("voxel 2 sum_rss row is subj1 + subj2",
      all.equal(agg$sum_rss[2, ], c(11, 10, 4, 7)),
      detail = paste("got", paste(agg$sum_rss[2, ], collapse = ",")))
check("voxel 1 sum_rss row is subj1 only",
      all.equal(agg$sum_rss[1, ], subj1_RSS[1, ]))
check("voxel 3 sum_rss row is subj2 only",
      all.equal(agg$sum_rss[3, ], subj2_RSS[3, ]))
check("voxel 4 sum_rss row is zero (nobody activated)",
      all(agg$sum_rss[4, ] == 0))
check("n_contributing: v1=1, v2=2, v3=1, v4=0",
      identical(agg$n_contributing, c(1L, 2L, 1L, 0L)))

# ---- TOY 2: pick_winning_grid_point argmin -----------------------------------
cat("\n--- TOY 2: argmin per voxel ---\n")

# All 4 voxels in pop mask (mask_prop_NA != NA for all).
mask_prop_NA <- c(0.5, 0.9, 0.8, NA)   # voxel 4 out of pop mask (n_contrib=0 anyway)
winning_k <- pick_winning_grid_point(agg$sum_rss, mask_prop_NA, agg$n_contributing, verbose = 0)

check("winning_k[1] = 2 (subj1 argmin)", winning_k[1] == 2L)
check("winning_k[2] = 3 (sum argmin)",   winning_k[2] == 3L)
check("winning_k[3] = 4 (subj2 argmin)", winning_k[3] == 4L)
check("winning_k[4] = NA (no contrib)",  is.na(winning_k[4]))

# ---- TOY 3: modal unmask -----------------------------------------------------
cat("\n--- TOY 3: modal unmask ---\n")

# Build a case where mode is unambiguous: 3 voxels activated with k=(2, 2, 3),
# 2 non-activated. Mode of activated is 2 (appears twice).
winning_k_toy <- c(2L, 2L, 3L, NA_integer_, NA_integer_)
mask_prop_toy <- c(0.5, 0.5, 0.5, NA, NA)
hrf_grid_toy  <- data.frame(a1 = 5:8, b1 = c(1,1,1,1), c = c(0,0,1/6,1/6))

filled <- modal_unmask(winning_k_toy, mask_prop_toy, hrf_grid_toy, verbose = 0)
check("modal fill applied to non-activated (voxel 4)", filled[4] == 2L)
check("modal fill applied to non-activated (voxel 5)", filled[5] == 2L)
check("activated voxels unchanged",
      identical(filled[1:3], c(2L, 2L, 3L)))

# Tie-breaking: which.max returns the FIRST index. So if winning_k = c(1,2), mode = 1.
tie_k     <- c(1L, 2L, NA_integer_)
tie_mask  <- c(0.5, 0.5, NA)
filled_tie <- modal_unmask(tie_k, tie_mask, hrf_grid_toy, verbose = 0)
check("tie-break: first mode (1) wins", filled_tie[3] == 1L)

# ---- TOY 4: warning when RSS is NULL ----------------------------------------
cat("\n--- TOY 4: subject with NULL RSS emits warning ---\n")

subj_null <- list(glm_result = list(mGLM0s = list(cortexL = list(RSS = NULL),
                                                    cortexR = list(RSS = NULL))))
w <- withCallingHandlers(
  sum_rss_across_subjects(list(subj_null, mk_subj(subj1_RSS)), subj_masks, n_candidates = 4, verbose = 0),
  warning = function(w) {
    assign("captured_warning", conditionMessage(w), envir = parent.frame(3))
    invokeRestart("muffleWarning")
  }
)
check("captured a save_rss warning",
      exists("captured_warning") && grepl("save_rss", captured_warning, ignore.case = TRUE),
      detail = if (exists("captured_warning")) captured_warning else "no warning")
# subject 2 in this call is mk_subj(subj1_RSS); its mask (subj_masks[,2]) is
# c(F, T, T, F) — so it contributes to voxels 2 and 3.
check("skipped subject with NULL RSS still let the second subject contribute (voxel 2)",
      all.equal(w$sum_rss[2, ], subj1_RSS[2, ]))
check("skipped subject with NULL RSS still let the second subject contribute (voxel 3)",
      all.equal(w$sum_rss[3, ], subj1_RSS[3, ]))
check("skipped subject contributed nothing at voxel 1 (its mask is FALSE there)",
      all(w$sum_rss[1, ] == 0))

# ---- TOY 5: full pipeline via public entrypoint -----------------------------
cat("\n--- TOY 5: regularize_allHRFs end-to-end on toy data ---\n")

# Build minimal workingHRF / allHRF result structures.
n_voxels <- 4L; n_candidates <- 4L
# activation_masks$masks: voxels x subjects logical matrix
# activation_masks$mask_prop_NA: length-n_voxels vector; NA = out of pop mask
fwh_toy <- list(
  activation_masks = list(
    masks         = subj_masks,
    mask_prop_NA  = c(0.5, 0.9, 0.8, NA)   # voxel 4 out of pop mask
  ),
  subject_results = NULL   # extract_xii_template returns NULL if not present
)

fah_toy <- list(
  subject_results = subject_results,
  hrf_grid = toy_grid
)

# Skip get_hrf_metrics dependency (needs HRF_calc); shim add_grid_metrics via
# a tiny stub for THIS run only.
orig_add <- add_grid_metrics
assignInNamespace("add_grid_metrics",
                  function(grid) { grid$time_to_peak <- grid$a1 - grid$b1; grid$FWHM <- grid$a1; grid },
                  ns = "hrf")

reg <- regularize_allHRFs(fwh_toy, fah_toy, verbose = 0)

assignInNamespace("add_grid_metrics", orig_add, ns = "hrf")

check("result has class regularizeHRFs", inherits(reg, "regularizeHRFs"))
check("pop_avg has 4 rows", nrow(reg$pop_avg) == 4L)
check("pop_avg cols correct (incl. t2p_mean/fwhm_mean for plot compat)",
      identical(sort(names(reg$pop_avg)),
                sort(c("voxel","a1","b1","c","t2p_mean","fwhm_mean"))))
check("voxel 1 winner = a1[k=2]",  reg$pop_avg$a1[1] == toy_grid$a1[2])
check("voxel 2 winner = a1[k=3]",  reg$pop_avg$a1[2] == toy_grid$a1[3])
check("voxel 3 winner = a1[k=4]",  reg$pop_avg$a1[3] == toy_grid$a1[4])
# voxel 4 was outside pop mask AND had no contributors, so modal fills it.
# Mode of (k=2, k=3, k=4) = 2 (first).
check("voxel 4 modal-filled to k=2", reg$pop_avg$a1[4] == toy_grid$a1[2])
check("winning_c is mode of c",
      reg$winning_c == as.numeric(names(sort(table(reg$pop_avg$c), decreasing = TRUE))[1]))

# ---- REAL DATA (fixture) -----------------------------------------------------
cat("\n--- Real fixture regression ---\n")

fixture_dir <- "~/Documents/Github/hrf-z/dev/fixtures"
fwh_path <- file.path(fixture_dir, "fit_workingHRF_result_motorlr_4s.rds")
fah_path <- file.path(fixture_dir, "fit_allHRFs_result_motorlr_4s.rds")

if (file.exists(fwh_path) && file.exists(fah_path)) {
  fah <- readRDS(fah_path)
  has_rss <- !is.null(fah$subject_results[[1]]$glm_result$mGLM0s$cortexL$RSS)
  if (!has_rss) {
    cat("  [SKIP] Fixture fit_allHRFs_result_motorlr_4s.rds was saved with RSS nulled.\n")
    cat("         Regenerate via `Rscript dev/tests/fit_allHRFs_test.R` with save_rss=TRUE to enable.\n")
  } else {
    fwh <- readRDS(fwh_path)
    library(ciftiTools)
    ciftiTools::ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
    reg_real <- regularize_allHRFs(fwh, fah, verbose = 0)
    check("real: class regularizeHRFs",   inherits(reg_real, "regularizeHRFs"))
    check("real: pop_avg is nonempty",    nrow(reg_real$pop_avg) > 0)
    check("real: pop_avg has snapped cols",
          all(c("a1","b1","c") %in% names(reg_real$pop_avg)))
    check("real: all winning a1 are on grid",
          all(reg_real$pop_avg$a1 %in% fah$hrf_grid$a1))
    check("real: all winning b1 are on grid",
          all(reg_real$pop_avg$b1 %in% fah$hrf_grid$b1))
    check("real: mask_prop_NA passed through",
          identical(reg_real$mask_prop_NA, fwh$activation_masks$mask_prop_NA))
    check("real: winning_k length = n_voxels",
          length(reg_real$winning_k) == length(reg_real$mask_prop_NA))
  }
} else {
  cat("  [SKIP] Fixture files not found at", fixture_dir, "\n")
}

# ---- Summary ---------------------------------------------------------------
cat("\n")
if (fail == 0) {
  cat("All checks passed.\n")
  quit(status = 0)
} else {
  cat("FAILED:", fail, "check(s).\n")
  quit(status = 1)
}
