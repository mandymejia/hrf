
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
plot_hrf_preview <- function(hrf_grid, EVs, nT, TR, hrf_idx = 1, onsets = FALSE, offsets = FALSE, ...) {
  design_dHRFs <- make_design(
    EVs = EVs, nTime = nT, TR = TR,
    onset = onsets, offset = offsets,
    a1 = hrf_grid$a1[hrf_idx], a2 = hrf_grid$a2[hrf_idx],
    b1 = hrf_grid$b1[hrf_idx], b2 = hrf_grid$b2[hrf_idx],
    c = hrf_grid$c[hrf_idx],
  )$design
  
  plot_design_fit_core(design_dHRFs, ...)
}