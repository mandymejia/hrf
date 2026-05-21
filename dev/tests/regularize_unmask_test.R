# Test: regularize_allHRFs unmask step.
# Runs regularize end-to-end (seffects=FALSE for speed) on the 4-subject
# MOTOR_LR fixture and verifies the new unmask step in Step 3.5:
#   - pop_avg expands to cover all valid cortex voxels
#   - pop_avg still carries the snapped a1/b1/c cols after re-snap
#   - candidate_maps reflect the expanded voxel set
# Also writes before/after PNGs for visual eyeballing.
#
# Run from repo root:
#   Rscript dev/tests/regularize_unmask_test.R

library(ciftiTools)
ciftiTools.setOption('wb_path', '/Applications/workbench/bin_macosxub/wb_command')
devtools::load_all('~/Documents/Github/hrf-z', quiet = TRUE)

OUT <- 'dev/tests/out/plots/regularize_unmask_test'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ** Inputs **
fwh <- readRDS('dev/fixtures/fit_workingHRF_result_motorlr_4s.rds')
fah <- readRDS('dev/fixtures/fit_allHRFs_result_motorlr_4s_norss.rds')
sd  <- readRDS('dev/fixtures/session_data_4s/session_data_motor_lr_4s.rds')

# ** Pre-unmask baseline (saved before unmask was introduced) **
baseline <- readRDS('dev/fixtures/regularize_allHRFs_result_motorlr_4s_norss.rds')

# ** Run regularize with the new unmask step **
cat('\n=== Running regularize_allHRFs (seffects=FALSE) ===\n')
got <- regularize_allHRFs(
  workingHRF_results = fwh,
  allHRF_results     = fah,
  BOLD     = sd$BOLD_files,
  EVs      = sd$EVs_list,
  nuisance = sd$nuisance_files,
  TR       = 0.72,
  onsets   = TRUE,
  offsets  = TRUE,
  seffects = FALSE,
  verbose  = 1
)

# ** Asserts **
fail <- 0
check <- function(name, ok, detail = '') {
  if (ok) cat('  [PASS]', name, '\n')
  else { cat('  [FAIL]', name, if (nzchar(detail)) paste0(' — ', detail) else '', '\n'); fail <<- fail + 1 }
}

cat('\n=== Unit checks ===\n')
expected_cols <- c('voxel', 't2p_mean', 'fwhm_mean', 'a1_snapped', 'b1_snapped', 'c_snapped')
check('pop_avg has all expected cols',
      all(expected_cols %in% colnames(got$pop_avg)),
      paste('cols:', paste(colnames(got$pop_avg), collapse = ',')))

check('pop_avg expanded vs baseline',
      nrow(got$pop_avg) > nrow(baseline$pop_avg),
      sprintf('got=%d baseline=%d', nrow(got$pop_avg), nrow(baseline$pop_avg)))

check('a1_snapped on grid',
      all(got$pop_avg$a1_snapped %in% unique(got$hrf_grid$a1)))
check('b1_snapped on grid',
      all(got$pop_avg$b1_snapped %in% unique(got$hrf_grid$b1)))
check('c_snapped is single value (winning_c)',
      length(unique(got$pop_avg$c_snapped)) == 1)

check('candidate_maps voxel count matches pop_avg',
      all(sapply(got$candidate_maps, nrow) == nrow(got$pop_avg)))

check('winning_c unchanged vs baseline',
      identical(got$winning_c, baseline$winning_c))

cat('\n=== Pop_avg row counts ===\n')
cat('  baseline (pre-unmask): ', nrow(baseline$pop_avg), '\n', sep = '')
cat('  got      (post-unmask):', nrow(got$pop_avg), '\n', sep = '')

# ** Plots **
xii_template <- extract_xii_template(fwh)
N <- length(got$mask_prop_NA)
build <- function(df, col) {
  v <- rep(NA, N); v[df$voxel] <- df[[col]]
  newdata_xifti(xii_template, v)
}

ZLIM  <- list(a1_snapped = c(4, 8), b1_snapped = c(0.5, 1.5), c_snapped = c(0, 0.17))
CMODE <- c(a1_snapped = 'diverging', b1_snapped = 'diverging', c_snapped = 'sequential')
mat   <- list(lit = TRUE, smooth = FALSE)

cat('\n=== Plots ===\n')
for (col in names(ZLIM)) {
  plot(build(baseline$pop_avg, col),
       zlim = ZLIM[[col]], color_mode = CMODE[[col]], material = mat,
       shadows = 1, NA_color = '#505560',
       title = paste0('motor_lr_4s | ', col, ' | BASELINE (no unmask)'),
       fname = file.path(OUT, paste0('motor_lr_4s_', col, '_baseline.png')))
  plot(build(got$pop_avg, col),
       zlim = ZLIM[[col]], color_mode = CMODE[[col]], material = mat,
       shadows = 1, NA_color = '#505560',
       title = paste0('motor_lr_4s | ', col, ' | NEW (median fill + smooth 5mm + snap)'),
       fname = file.path(OUT, paste0('motor_lr_4s_', col, '_new.png')))
}

cat('\n', if (fail == 0) 'ALL CHECKS PASSED' else paste(fail, 'CHECK(S) FAILED'), '\n', sep = '')
cat('6 PNGs in', OUT, '\n')
if (fail > 0) quit(status = 1)
