
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hrf

<!-- # hrf <img src="man/figures/logo.png" align="right" height="139" alt="hrf sticker" /> -->
<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/hrf)](https://cran.r-project.org/package=hrf)
[![R-CMD-check](https://github.com/mandymejia/hrf/workflows/R-CMD-check/badge.svg)](https://github.com/mandymejia/hrf/actions)
<!-- [![Codecov test coverage](https://codecov.io/gh/mandymejia/hrf/branch/master/graph/badge.svg)](https://app.codecov.io/gh/mandymejia/hrf?branch=master) -->
<!-- badges: end -->

The `hrf` R package includes these main functions:

- `HRF_main` for HRF computation
- `multiGLM_fun` for testing different design matrices
- `multiGLM` for running `multiGLM_fun` on CIFTI data

## Important Note on Dependencies:

`hrf::multiGLM` requires the `ciftiTools` package, which requires an
installation of Connectome Workbench. It can be installed from the [HCP
website](https://www.humanconnectome.org/software/get-connectome-workbench).
