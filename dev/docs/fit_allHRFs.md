# `fit_allHRFs()` — arguments

Fits a per-subject GLM across a grid of HRFs and finds the best fit per voxel.

## Signature

```r
fit_allHRFs(
  BOLD, EVs, TR,
  hrf_grid = generate_hrf_grid,
  brainstructures = c("left", "right"),
  resamp_res = NULL,
  hpf = 0.01,
  nuisance = NULL,
  onsets = TRUE, offsets = TRUE,
  scrub = NULL,
  smoothing = TRUE, surf_FWHM = 5,
  verbose = 1, n_cores = 1,
  save_rss = FALSE,
  log_dir = "logs", work_dir = NULL
)
```

## Arguments

| Arg | What it is |
|---|---|
| `BOLD`, `EVs`, `TR` | same as `fit_workingHRF` |
| `hrf_grid` | data.frame of HRF parameter combos (`a1`, `b1`, `c`, `a2`, `b2`); default = `generate_hrf_grid()` (~100+ rows) |
| `brainstructures`, `resamp_res`, `hpf`, `nuisance` | same as `fit_workingHRF` |
| `onsets`, `offsets`, `scrub` | same as `fit_workingHRF` |
| `smoothing`, `surf_FWHM` | optional surface smoothing |
| `n_cores` | parallel workers |
| `save_rss` | **important.** If `TRUE`, save per-voxel RSS for every grid HRF — lets `regularize_allHRFs` skip refits later |
| `log_dir` | parallel cluster logs |
| `work_dir` | per-subject `.qs` files dumped here (each holds GLM results, design matrices, optional RSS) |

## Values used in `1a_analysis.R`

```r
fit_allHRFs(
  BOLD = session_data$BOLD_files[seq_len(n_use)],
  EVs  = session_data$EVs_list[seq_len(n_use)],
  TR = 0.72,
  brainstructures = c('left', 'right'),
  resamp_res = 10000,
  hpf = 0.01,
  nuisance = session_data$nuisance_files[seq_len(n_use)],
  onsets  = if(t == 3) FALSE else TRUE,
  offsets = if(t == 3) FALSE else TRUE,
  scrub = NULL,
  save_rss = TRUE,
  verbose = 2,
  n_cores = 32,
  log_dir  = safe_path(log_root,  "fit_allHRFs", paste0(sess, "_", n_use, "s")),
  work_dir = safe_path(output_root, "fit_allHRFs", "work", paste0(sess, "_", n_use, "s"))
)
```

Notable: `save_rss = TRUE` is what unlocks the fast lookup mode in `regularize_allHRFs`. `work_dir` ends up holding `~n_subjects` `.qs` files.
