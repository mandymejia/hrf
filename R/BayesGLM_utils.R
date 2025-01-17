#' Scale the design matrix
#'
#' @param design_mat The original (unscaled) design matrix that is T x K, where
#'     T is the number of time points, and k is the number of field covariates
#'
#' @return A scaled design matrix
#'
#' @keywords internal
#'
scale_design_mat <- function(design_mat) {
  stopifnot(is.matrix(design_mat))
  apply(design_mat,2,function(field) {
    #(field - mean(field)) / max(abs(field))
    field / max((field))
  })
}

#' Bayes GLM arg checks
#'
#' Checks arguments for \code{BayesGLM} and \code{fit_bayesglm}
#'
#' Avoids duplicated code between \code{BayesGLM} and \code{fit_bayesglm}
#'
#' @param scale_BOLD See \code{BayesGLM}.
#' @param Bayes,EM See \code{BayesGLM}.
#' @param ar_order,ar_smooth,aic See \code{BayesGLM}.
#' @param n_threads See \code{BayesGLM}.
#' @param return_INLA See \code{BayesGLM}.
#' @param verbose See \code{BayesGLM}.
# @param combine_sessions See \code{BayesGLM}.
#' @param meanTol,varTol,emTol See \code{BayesGLM}.
#'
#' @return The arguments that may have changed, in a list: \code{scale_BOLD},
#'  \code{do_Bayesian}, \code{do_EM}, and \code{do_pw}.
#'
#' @importFrom fMRItools is_1
#' @keywords internal
BayesGLM_argChecks <- function(
    scale_BOLD=c("mean", "sd", "none"),
    Bayes=TRUE,
    EM=FALSE,
    ar_order=6,
    ar_smooth=5,
    aic=FALSE,
    n_threads=4,
    return_INLA=c("trimmed", "full", "minimal"),
    verbose=1,
    meanTol=1e-6,
    varTol=1e-6,
    emTol=1e-3
){

  if (isTRUE(scale_BOLD)) {
    message("Setting `scale_BOLD` to 'mean'"); scale_BOLD <- "mean"
  }
  if (isFALSE(scale_BOLD)) {
    message("Setting `scale_BOLD` to 'none'"); scale_BOLD <- "none"
  }
  scale_BOLD <- match.arg(scale_BOLD, c("mean", "sd", "none"))

  stopifnot(fMRItools::is_1(Bayes, "logical"))
  stopifnot(fMRItools::is_1(EM, "logical"))
  if (EM && !Bayes) {
    warning("EM is a Bayesian method: setting `Bayes` to `TRUE`.")
    Bayes <- TRUE
  }
  if (Bayes) {
    stop()
  }

  if(EM) stop("EM not available.") #not currently available

  if (is.null(ar_order)) ar_order <- 0
  stopifnot(fMRItools::is_1(ar_order, "numeric"))
  do_pw <- ar_order > 0
  if (is.null(ar_smooth)) ar_smooth <- 0
  stopifnot(fMRItools::is_1(ar_smooth, "numeric"))
  stopifnot(fMRItools::is_1(aic, "logical"))

  stopifnot(fMRItools::is_1(n_threads, "numeric"))
  stopifnot(n_threads <= parallel::detectCores())

  if (isTRUE(return_INLA)) {
    message("Setting `return_INLA` to 'trimmed'"); return_INLA <- "trimmed"
  }
  if (isFALSE(return_INLA)) {
    message("Setting `return_INLA` to 'minimal'"); return_INLA <- "minimal"
  }
  return_INLA <- match.arg(return_INLA, c("trimmed", "full", "minimal"))

  if (isTRUE(verbose)) { verbose <- 2 }
  if (isFALSE(verbose)) { verbose <- 0 }
  stopifnot(fMRItools::is_posNum(verbose, zero_ok=TRUE))

  stopifnot(fMRItools::is_posNum(meanTol))
  stopifnot(fMRItools::is_posNum(varTol))
  stopifnot(fMRItools::is_posNum(emTol))

  # Return new parameters, and parameters that may have changed.
  list(
    scale_BOLD=scale_BOLD,
    Bayes=Bayes,
    EM = EM,
    do_pw = do_pw,
    return_INLA=return_INLA
  )
}
