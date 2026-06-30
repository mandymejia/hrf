#' Plot method for combined hrfs objects
#'
#' Routes \code{type} to the matching sub-object plot method:
#' \itemize{
#'   \item working-HRF: \code{"proportion"}, \code{"binary"}, \code{"mask"}
#'   \item allHRF:      \code{"hrfs"}, \code{"param_grid"}, \code{"single_hrf"},
#'                       \code{"multiple_hrf"}
#'   \item regularize:  \code{"pop_avg"}, \code{"pop_avg_continuous"},
#'                       \code{"mean"}, \code{"mean_all"}, \code{"param_heatmap"}
#' }
#' \code{type = "design"} is defined by both workingHRF AND allHRFs; pass
#' \code{which = "working"} (default) or \code{which = "all"} to disambiguate.
#'
#' @param x A \code{"hrfs"} object from \code{\link{fit_allHRFs}}.
#' @param type Character plot type (see above).
#' @param which Used only when \code{type = "design"}: \code{"working"} or \code{"all"}.
#' @param ... Forwarded to the sub-object plot method.
#'
#' @return Invisibly returns the underlying plot result.
#' @export
plot.hrfs <- function(x, type, which = c("working", "all"), ...) {
  if (identical(type, "design")) {
    which <- match.arg(which)
    sub <- if (which == "working") "fit_workingHRF" else "fit_allHRFs"
    return(plot(x[[sub]], type = "design", ...))
  }

  routing <- c(
    proportion         = "fit_workingHRF",
    binary             = "fit_workingHRF",
    mask               = "fit_workingHRF",
    hrfs               = "fit_allHRFs",
    param_grid         = "fit_allHRFs",
    single_hrf         = "fit_allHRFs",
    multiple_hrf       = "fit_allHRFs",
    pop_avg            = "regularize_allHRFs",
    pop_avg_continuous = "regularize_allHRFs",
    mean_all           = "regularize_allHRFs",
    mean               = "regularize_allHRFs",
    param_heatmap      = "regularize_allHRFs"
  )

  sub <- routing[[type]]
  if (is.null(sub)) {
    stop("plot.hrfs: unknown type '", type,
         "'. Valid: design, ", paste(names(routing), collapse = ", "))
  }
  plot(x[[sub]], type = type, ...)
}
