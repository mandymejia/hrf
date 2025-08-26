#' Declare global imports and known global variables
#'
#' This file is used to manage package-level imports and to suppress 
#' R CMD check notes about symbols that are known to exist but are 
#' not visible to static analysis (e.g., internal functions from other packages).
#'
#' @importFrom stats runif
#' @importFrom utils read.table
NULL

# Suppress R CMD check notes for known non-exported functions from ciftiTools
utils::globalVariables(c(
  "ciftiTools.getOption",
  "ciftiTools.setOption"
))