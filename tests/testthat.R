# Build --> Install and Restart
# [Edit these] paths to Workbench and result output.
my_wb <- "../workbench"
save_path <- "tests/output" # in gitignore

library(testthat)
library(hrf)

tests_dir <- "testthat"
if (!endsWith(getwd(), "tests")) { tests_dir <- file.path("tests", tests_dir) }

source(file.path(tests_dir, "test-auto.R"))
