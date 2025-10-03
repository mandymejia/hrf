#' Plot method for regularizeHRFs objects
#'
#' Creates diagnostic plots for HRF regularization results.
#' Currently supports visualization of the population activation mask
#' (\code{mask_prop_NA}).
#'
#' @param x An object of class \code{"regularizeHRFs"} from
#'   \code{\link{regularize_allHRFs}}.
#' @param type Character. Type of plot to generate. Options:
#'   \itemize{
#'     \item \code{"mask"} – plot the group-level mask (mask_prop_NA).
#'   }
#' @param ... Additional arguments passed to internal plotting functions.
#'
#' @return Invisibly returns the result of the plotting call.
#' @export
plot.regularizeHRFs <- function(x, type = c("mask"), ...) {
  type <- match.arg(type)

  switch(type,
         mask = plot_mask_prop(x, ...)
  )
}

#' Plot mask_prop_NA (population activation mask)
#'
#' Visualize the population-level activation mask across subjects,
#' using a provided xifti template for surface geometry.
#'
#' @param x A \code{regularizeHRFs} object containing \code{mask_prop_NA}.
#' @param xii A template \code{xifti} object (usually stored as an attribute of \code{x}).
#' @param fname Optional output filename. If provided, the plot is saved.
#' @param title Plot title. Default is 'Subjects with Significant "Working HRF" Activation > 10%'.
#' @param legend_fname Optional filename to save the color legend separately.
#' @param shadows Numeric, shadow depth for rendering. Default is 1.
#' @param material List of rendering options. Default is \code{list(lit = TRUE, smooth = FALSE)}.
#' @param NA_color Color for NA values and medial wall. Default is "#505560".
#' @param ... Additional arguments passed to \code{ciftiTools::view_xifti_surface}.
#'
#' @return Invisibly returns the plot object.
#' @importFrom ciftiTools newdata_xifti
#' @importFrom ciftiTools view_xifti_surface
#' @keywords internal
plot_mask_prop <- function(x,
                           xii,
                           fname = NULL,
                           title = NULL,
                           legend_fname = NULL,
                           shadows = 1,
                           material = list(lit = TRUE, smooth = FALSE),
                           NA_color="#505560",
                           ...  # shadows, material, etc will pass through
) {
  stopifnot("mask_prop_NA not found in object" = "mask_prop_NA" %in% names(x))

  xii <- attr(x, "xii")

  thr <- if ("threshold" %in% names(x)) x$threshold else 0.10

  if (is.null(title)) {
    title <- paste0("Subjects with Significant \"Working HRF\" Activation > ", thr * 100, "%")
  }

  plot_obj <- plot(
    newdata_xifti(xii, x$mask_prop_NA),
    title = title,
    color_mode = "qualitative",
    fname = fname,
    legend_fname = legend_fname,
    shadows = 1,
    material = list(lit = TRUE, smooth = FALSE),
    NA_color="#505560",
    ...
  )
  invisible(plot_obj)
}
