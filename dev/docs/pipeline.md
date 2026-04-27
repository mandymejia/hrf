# HRF Adaptation Pipeline

Four stages, run in order. Outputs feed forward.

```
BOLD + EVs ──┬──→ fit_workingHRF ──→ activation masks ────┐
             │                                             ├──→ regularize_allHRFs ──→ HRF map ──→ fit_bestHRF ──→ betas, t-stats, p-values
             └──→ fit_allHRFs    ──→ best HRF per voxel ──┘
```

## `fit_workingHRF`
Fits a GLM per subject using a **single canonical HRF** (`a1=6, b1=1, c=1/6`).
Identifies voxels with significant task activation via F-stats.

- **In:** BOLD CIFTI, event timings, nuisance, TR
- **Out:** per-subject GLM results, activation masks

## `fit_allHRFs`
Fits a GLM per subject across a **grid of HRF parameterizations** (~100+).
For each voxel × subject, finds the best-fitting HRF in the grid.

- **In:** same BOLD/events, plus an HRF grid
- **Out:** best `(a1, b1, c)` per voxel × subject (or RSS table for fast lookup)

## `regularize_allHRFs`
Combines both upstream outputs to produce a smoothed, population-aware HRF map.

1. Uses `fit_workingHRF` masks to keep only active voxels per subject
2. Uses `fit_allHRFs` best params to compute a population-averaged HRF (per voxel)
3. Optionally generates 25 candidate maps (offsets around population avg) and picks the best per subject

- **In:** `fit_workingHRF` result + `fit_allHRFs` result
- **Out:** `pop_avg` HRF map, per-subject candidate maps (if `seffects=TRUE`)

## `fit_bestHRF`
Final per-subject GLM using each voxel's personalized HRF. Computes contrasts, t-stats, and adjusted p-values.

- **In:** `regularize_allHRFs` result + BOLD/events/contrasts
- **Out:** `bestHRF` object with `working` and `population`/`adapted` sections, each holding betas, contrasts (est, SE, t, p, p-adj), and HRF assignments
