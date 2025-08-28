# ======================================================================
# hrf-package-globals.R
# Centralized package-wide imports, dynamic variable declarations,
# and Rcpp integration settings for the 'hrf' package.
# ======================================================================

#' @name hrf-package-globals
#' @title Global Imports and Known Variables
#' @description
#' This file centralizes package-wide imports, Rcpp integration,
#' and declares known variable names used inside non-standard evaluation (NSE)
#' contexts. NSE occurs in functions like `subset()`, `transform()`, and
#' dynamically generated column names in data frames. Declaring these names here
#' informs R CMD check that they are intentional and suppresses false-positive NOTES.
#'
#' @importFrom stats runif
#' @importFrom utils read.table
#' @importFrom utils object.size
#' @importFrom dplyr group_by summarize left_join filter
#' @importFrom magrittr %>%
NULL

# ======================================================================
# Rcpp Integration for Compiled Code
# ======================================================================

#' @useDynLib hrf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
# Explanation:
# - @useDynLib ensures that the compiled shared library for the package
#   (generated from src/*.cpp) is loaded automatically.
# - .registration = TRUE enables symbol registration for better performance
#   and compliance with CRAN best practices.
# - @importFrom Rcpp sourceCpp ensures Rcpp’s generated wrappers work correctly.

# ======================================================================
# Declare Known Variables and Functions for R CMD check
# ======================================================================

# 1. Internal ciftiTools helpers
utils::globalVariables(c(
  "ciftiTools.getOption",
  "ciftiTools.setOption"
))

# 2. Dynamically Generated Column Names
utils::globalVariables(c(
  "time_to_peak",
  "peak2_time",
  "voxel",
  "mask",
  "param0",
  "param_mean0"
))