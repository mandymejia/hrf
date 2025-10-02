#' Plot method for allHRFs objects
#'
#' Creates diagnostic plots for grid-based HRF fitting results. Supports
#' design matrix visualization, and HRF grid plots (raw, tapered, parameter grid).
#'
#' @param x An object of class \code{"allHRFs"} from \code{\link{fit_allHRFs}}.
#' @param type Character. Type of plot to generate. Options are:
#'   \itemize{
#'     \item \code{"design"} – plot subject design matrix.
#'     \item \code{"hrfs"} – plot all raw HRFs from the parameter grid.
#'     \item \code{"hrfs_tapered"} – plot tapered HRFs.
#'     \item \code{"param_grid"} – plot HRF parameter grid heatmaps.
#'   }
#' @param subject Integer. Subject index for design matrix plots.
#' @param hrf_idx Integer. HRF index for design matrix plots.
#' @param ... Additional arguments passed on.
#'
#' @return Invisibly returns the result of the specific plotting function.
#' @export
plot.allHRFs <- function(x, type = c("design", "hrfs", "hrfs_tapered", "param_grid"),
                         subject = 1, hrf_idx = 1, ...) {
  type <- match.arg(type)

  switch(type,
         design       = plot_design_fit_all(x, hrf_idx = hrf_idx, subject = subject, ...),
         hrfs         = plot.hrf_grid(x$hrf_grid, type = "hrfs", ...),
         hrfs_tapered = plot.hrf_grid(x$hrf_grid, type = "hrfs_tapered", ...),
         param_grid   = plot.hrf_grid(x$hrf_grid, type = "param_grid", ...)
  )
}

#' Plot design matrix for specific HRF parameter combination
#'
#' Internal function to visualize the design matrix for a specific subject and
#' HRF parameter combination from allHRFs results. Loads the saved design data
#' and renders it using the core plotting infrastructure.
#'
#' @param x An object of class \code{"allHRFs"} from \code{\link{fit_allHRFs}}.
#' @param hrf_idx Integer. Index of the HRF parameter combination to plot.
#' @param subject Integer. Index of the subject to plot.
#' @param ... Additional arguments passed to \code{\link{plot_design_fit_core}}.
#'
#' @return Result from \code{\link{plot_design_fit_core}}.
#'
#' @keywords internal
plot_design_fit_all <- function (x, hrf_idx = 1, subject = 1, ...) {
  # Get file paths from the allHRFs result
  file_paths <- attr(x, "result_paths")

  # Load the full subject data (including design_3D)
  full_subject_data <- qs::qread(file_paths[[subject]])

  # Extract the specific design matrix for the requested HRF index
  design_dHRFs <- full_subject_data[["design_3D"]][["list"]][[hrf_idx]]

  plot_design_fit_core(design_dHRFs, ...)
}




#' Preview HRF Design Matrix from Parameter Grid
#'
#' This function previews the design matrix for a single HRF configuration from a grid,
#' including optional onset/offset modeling. It is primarily used for inspecting how the
#' modeled regressors vary with different HRF parameters.
#'
#' @param hrf_grid A data frame or list with HRF parameter combinations. Must contain
#'   columns or list elements named \code{a1}, \code{a2}, \code{b1}, \code{b2}, and \code{c}.
#' @param hrf_idx Integer index into the HRF grid specifying which HRF to preview.
#' @inheritParams EVs_Param
#' @inheritParams nT_Param
#' @inheritParams onsets_Param
#' @inheritParams offsets_Param
#' @param TR Repetition time (TR) in seconds.
#' @param ... Passed to \code{\link[ggplot2]{theme}}, \code{\link[ggplot2]{facet_grid}},
#'   or other internal ggplot functions for additional customization.
#'
#' @return A named list of \code{ggplot2} objects: \code{main}, \code{dhrf}, and \code{ddhrf}
#'   (if applicable). Also prints the plots.
#'
#' @import ggplot2
#' @export
plot_hrf_preview <- function(hrf_grid, EVs, nT, TR, hrf_idx = 1, onsets = TRUE, offsets = TRUE, ...) {
  design_dHRFs <- make_design(
    EVs = EVs, nTime = nT, TR = TR,
    onset = onsets, offset = offsets,
    a1 = hrf_grid$a1[hrf_idx], a2 = hrf_grid$a2[hrf_idx],
    b1 = hrf_grid$b1[hrf_idx], b2 = hrf_grid$b2[hrf_idx],
    c = hrf_grid$c[hrf_idx],
  )$design

  plot_design_fit_core(design_dHRFs, ...)
}
