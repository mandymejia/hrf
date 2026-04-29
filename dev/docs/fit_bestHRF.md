# `fit_bestHRF()` — arguments

Final per-subject GLM using the adapted (or population) HRF at each voxel. Returns betas, contrasts, t-stats, p-values.

## Signature

```r
fit_bestHRF(
  regularize_result,
  BOLD_file,
  EVs,
  nuisance_file = NULL,
  TR,
  use = c("subject", "population"),
  subject_idx = NULL,
  contrasts = NULL,
  working_hrf = list(a1 = 6, b1 = 1, c = 1/6),
  brainstructures = c("left", "right"),
  resamp_res = 10000,
  hpf = 0.01,
  onsets = FALSE, offsets = FALSE,
  p_adjust_method = "BH",
  verbose = 1
)
```

## Arguments

| Arg | What it is |
|---|---|
| `regularize_result` | full result from `regularize_allHRFs()` |
| `BOLD_file` | **single** CIFTI file path for the subject you want to fit |
| `EVs` | event table for that subject |
| `nuisance_file` | optional path to nuisance regressors for that subject |
| `TR` | repetition time (must match upstream stages) |
| `use` | `"subject"` (default) → use that subject's adapted HRF map; `"population"` → use the population-average map. `"subject"` errors out if regularize was run with `seffects = FALSE` |
| `subject_idx` | required when `use = "subject"`. Picks which subject's adapted HRF map to use |
| `contrasts` | optional contrast matrix `A` (rows = contrasts, cols = task regressors). `NULL` = identity (each beta tested individually) |
| `working_hrf` | the canonical HRF used for the `working` baseline section |
| `brainstructures`, `resamp_res`, `hpf` | must match upstream stages |
| `onsets`, `offsets` | must match the values used in `fit_allHRFs()` and `regularize_allHRFs()` |
| `p_adjust_method` | passed to `p.adjust()` — `"BH"`, `"bonferroni"`, etc. |

## What you need to provide

To run for a single subject:

1. **Pipeline output:** `regularize_result <- readRDS(".../regularize_allHRFs_result_<sess>_<n>s.rds")`
2. **The subject's BOLD file:** path to `.dtseries.nii`
3. **The subject's EVs:** event table for that scan
4. **Their nuisance file** (optional but matches upstream)
5. **Subject index:** `subject_idx` matching the position in the original session_data subject list
6. **Contrasts** (optional): a numeric matrix where each row is one contrast across task regressors. Default tests each regressor individually.

## Example

```r
regularize_result <- readRDS("/N/project/hrf_adaptation/validation_final/regularize_allHRFs/regularize_allHRFs_result_tfMRI_MOTOR_LR_1080s.rds")

# Custom contrast: regressor 1 minus regressor 2
A <- matrix(c(1, -1, 0, 0, 0, 0), nrow = 1)

result <- fit_bestHRF(
  regularize_result = regularize_result,
  BOLD_file = "/path/to/subject_001/MOTOR_LR.dtseries.nii",
  EVs = session_data$EVs_list[[1]],
  nuisance_file = session_data$nuisance_files[1],
  TR = 0.72,
  use = "subject",
  subject_idx = 1,
  contrasts = A,
  brainstructures = c("left", "right"),
  resamp_res = 10000,
  hpf = 0.01,
  onsets = TRUE, offsets = TRUE
)
```

## What you get back

A `bestHRF` object — see [pipeline.md](pipeline.md#structure-of-the-besthrf-return-object) for the full tree.
