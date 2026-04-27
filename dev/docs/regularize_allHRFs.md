# `regularize_allHRFs()` — arguments

Builds the population-averaged HRF map and (optionally) the per-subject **adapted HRF** maps.

## Signature

```r
regularize_allHRFs(
  workingHRF_results,
  allHRF_results,
  BOLD = NULL, EVs = NULL, nuisance = NULL, TR = NULL,
  a1_offsets = c(-2, -1, 0, 1, 2),
  b1_offsets = c(-0.5, -0.25, 0, 0.25, 0.5),
  seffects = TRUE,
  onsets = FALSE, offsets = FALSE,
  verbose = 1
)
```

## Arguments

| Arg | What it is |
|---|---|
| `workingHRF_results` | full result object returned by `fit_workingHRF()` |
| `allHRF_results` | full result object returned by `fit_allHRFs()` |
| `BOLD`, `EVs`, `nuisance`, `TR` | only required when refit mode is used (i.e. `allHRF_results` was made with `save_rss = FALSE`). In RSS-lookup mode these are ignored |
| `a1_offsets` | shifts in `a1` (time-to-peak) used to build candidate maps |
| `b1_offsets` | shifts in `b1` (width) used to build candidate maps |
| `seffects` | if `TRUE`, fit all 25 candidate maps to each subject and pick the winner per subject (= the **adapted HRF**). If `FALSE`, only the population avg is returned |
| `onsets`, `offsets` | only used in refit mode (must match the values used during `fit_allHRFs`) |
| `verbose` | 0 = silent, 1 = step messages, 2 = detailed |

## Mode auto-detection

`regularize_allHRFs` reads `allHRF_results$call_info$save_rss`:

- `save_rss == TRUE` → **lookup mode** (fast, reads pre-computed RSS from `.qs`)
- `save_rss == FALSE` → **refit mode** (slow, re-runs multiGLM per candidate; needs `BOLD`/`EVs`/`TR`)

## Values used in `1b_analysis.R`

```r
workingHRF_results <- readRDS(file.path(output_root, "fit_workingHRF",
  paste0("fit_workingHRF_result_", sess, "_", n_use, "s.rds")))
allHRF_results <- readRDS(file.path(output_root, "fit_allHRFs",
  paste0("fit_allHRFs_result_", sess, "_", n_use, "s.rds")))

regularize_allHRFs(
  workingHRF_results,
  allHRF_results,
  seffects = TRUE
)
```

That's it. Defaults handle everything else because we always run in `save_rss = TRUE` lookup mode.
