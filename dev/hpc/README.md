# Running `hrf` on HPC

Sample: `fit_bestHRF_sample.R` runs `fit_bestHRF()` on one subject and saves working + population beta plots.

## 1. Launch RStudio

From a ThinkLinc terminal — **not** the Applications menu (that loads an old R version):

```bash
module load rstudio/2025.04
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
```
