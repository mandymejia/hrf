# TODO - Package Preparation
- Remove toctoc from DESCRIPTION suggests
- Remove any references to tictoc:: (fit_ and regularize_)

# Notes
- If you ever need to modify `default_hrf_grid` you must run
```R
default_hrf_grid <- generate_default_hrf_grid()
usethis::use_data(default_hrf_grid, overwrite = TRUE)
```