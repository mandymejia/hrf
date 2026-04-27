# `fit_workingHRF()` — arguments

Fits a per-subject GLM with a single canonical HRF and identifies activated voxels.

## Signature

```r
fit_workingHRF(
  BOLD, EVs, TR,
  brainstructures = c("left", "right"),
  resamp_res = NULL,
  hpf = 0.01,
  nuisance = NULL,
  hrf_params = list(a1 = 6, b1 = 1, c = 1/6, a2 = 16, b2 = 1),
  derivatives = TRUE,
  onsets = TRUE, offsets = TRUE,
  scrub = NULL,
  alpha = 0.01,
  min_active_subjects = 20,
  smoothing = TRUE, surf_FWHM = 5,
  verbose = 1, n_cores = 1,
  log_dir = "logs"
)
```

## Arguments

| Arg | What it is |
|---|---|
| `BOLD` | character vector of CIFTI file paths, one per subject |
| `EVs` | list of event tables (one per subject) |
| `TR` | repetition time (seconds) |
| `brainstructures` | which surfaces to model: `"left"`, `"right"`, both |
| `resamp_res` | CIFTI resampling resolution (NULL = no resample) |
| `hpf` | high-pass filter cutoff in Hz (DCT bases) |
| `nuisance` | optional vector of nuisance regressor file paths |
| `hrf_params` | canonical HRF — defaults match SPM |
| `derivatives` | include time + dispersion derivatives in design |
| `onsets`, `offsets` | model event onsets/offsets as separate regressors |
| `scrub` | optional list of timepoint indices to scrub |
| `alpha` | significance threshold for activation mask |
| `min_active_subjects` | population mask requires at least this many subjects active per voxel |
| `smoothing`, `surf_FWHM` | optional surface smoothing |
| `n_cores` | parallel workers |
| `log_dir` | where parallel cluster logs go |

## Values used in `1a_analysis.R`

```r
fit_workingHRF(
  BOLD = session_data$BOLD_files[seq_len(n_use)],
  EVs  = session_data$EVs_list[seq_len(n_use)],
  TR = 0.72,
  brainstructures = c('left', 'right'),
  resamp_res = 10000,
  hpf = 0.01,
  nuisance = session_data$nuisance_files[seq_len(n_use)],
  hrf_params = list(a1 = 6, b1 = 1, c = 1/6, a2 = 16, b2 = 1),
  derivatives = TRUE,
  scrub = NULL,
  onsets  = if(t == 3) FALSE else TRUE,   # gambling has no onsets
  offsets = if(t == 3) FALSE else TRUE,   # gambling has no offsets
  alpha = 0.05,
  min_active_subjects = 20,
  verbose = 2,
  n_cores = 32,
  log_dir = safe_path(log_root, "fit_workingHRF", paste0(sess, "_", n_use, "s"))
)
```

Notable: `alpha = 0.05` (loose mask, lots of voxels pass), `n_cores = 32`, `resamp_res = 10000` (matches HCP defaults).
