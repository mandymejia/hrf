#' @name hrf-package-globals
#' @title Global Imports and Known Variables
#' @description
#' This file centralizes package-wide imports and declares known variable names
#' used inside non-standard evaluation (NSE) contexts, such as `subset()`, `transform()`,
#' or dynamically generated columns. Declaring them here informs R CMD check
#' that these variables are intentionally used within the package.
#'
#' @importFrom stats runif
#' @importFrom utils read.table
#' @importFrom utils object.size
NULL

# ======================================================================
# Declare known variables and functions for R CMD check
# ======================================================================

# 1. Internal ciftiTools helpers
# These functions are intentionally used within the package but are not exported
# from ciftiTools, so we explicitly declare them here for clarity.
utils::globalVariables(c(
  "ciftiTools.getOption",
  "ciftiTools.setOption"
))

# 2. Data frame column names created dynamically in generate_default_hrf_grid()
# These names are generated programmatically, so we explicitly register them
# to avoid false-positive NOTES from R CMD check.
utils::globalVariables(c(
  "time_to_peak",
  "peak2_time"
))