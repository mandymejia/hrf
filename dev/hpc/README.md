# Running `hrf` on HPC

Samples in this folder:
- `fit_bestHRF_sample.R` — basic run with `use = "population"`, saves working + population beta plots
- `fit_bestHRF_sample_contrast.R` — same but with a custom `hands - feet` contrast, saves t-stat plots

## `use` argument (fit_bestHRF)

- `use = "subject"` — uses that subject's adapted HRF map (the winning candidate from `regularize_allHRFs`). Requires `subject_idx` and the regularize result must have been built with `seffects = TRUE`.
- `use = "population"` — uses the population-average HRF map for everyone. Doesn't need `subject_idx`. Works on any regularize result.

Both samples here use `"population"` because it's the simplest sanity check and runs against any regularize output. Switch to `"subject"` + `subject_idx = N` once you want per-subject adapted fits. See [`../docs/fit_bestHRF.md`](../docs/fit_bestHRF.md) for the full arg breakdown.

## 1. Launch RStudio

From a ThinkLinc terminal — **not** the Applications menu (that loads an old R version):

```bash
module purge
module load StdEnv
module load gnu/12.2.0
module load r/4.5.1
module load rstudio/2025.04
module load java/15.0.2
module load mvapich/2.3.5
module load zlib/1.2.13
module load udunits/2.2.28
module load hdf5/1.12
module load netcdf/4.7
module load proj/8.0.1
module load gdal/3.6.4
module load gsl/2.6
module load sqlite/3.35.5
module load geos/3.8.2
rstudio
```

## 2. Install the package

In the RStudio console (one time):

```r
devtools::install("~/Documents/Github/hrf-z")
```

This pulls in every CRAN dep automatically. **Use `install`, not `load_all`** — parallel workers can only see packages that are properly installed.

## 3. Run

```r
source("~/Documents/Github/hrf-z/dev/hpc/fit_bestHRF_sample.R")
# or:
source("~/Documents/Github/hrf-z/dev/hpc/fit_bestHRF_sample_contrast.R")
```
