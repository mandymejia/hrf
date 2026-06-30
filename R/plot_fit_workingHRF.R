#' Plot method for `workingHRF` objects
#'
#' Dispatches to one of four plot types for visualizing a \code{workingHRF}
#' sub-object (e.g. \code{combo$fit_workingHRF} from \code{fit_allHRFs}):
#'
#' \describe{
#'   \item{\code{type = "design"}}{Plots the design matrix for a specific subject.}
#'   \item{\code{type = "proportion"}}{Displays the cross-subject activation proportion map.}
#'   \item{\code{type = "binary"}}{Displays the binary activation mask after thresholding.}
#'   \item{\code{type = "mask"}}{Displays the population-level activation mask (mask_prop_NA).}
#' }
#'
#' @param x A \code{workingHRF} object.
#' @param type Type of plot to display. One of \code{"design"}, \code{"proportion"}, \code{"binary"}, or \code{"mask"}.
#' @param subject Integer. Subject index for the design matrix plot (only used if \code{type = "design"}).
#' @param ... Additional arguments passed to the specific plot function. If you wish to use you must set all positional arguments.
#'
#' @return A plot is rendered; invisibly returns the underlying object used for plotting.
#'
#' @examples
#' \dontrun{
#' combo <- fit_allHRFs(...)
#' plot(combo$fit_workingHRF, type = "design", subject = 1)
#' plot(result, type = "proportion")
#' plot(result, type = "binary", threshold = 0.2)
#' plot(result, type = "mask")
#' }
#'
#' @export
plot.workingHRF <- function(x, type = c("design", "proportion", "binary", "mask"), subject = 1, ...) {
  type <- match.arg(type)
  cat("Plotting type:", type, "\n")

  switch(type,
    design = plot_design_fit(x, subject = subject, ...),
    proportion = plot_activation_proportion(x, ...),
    binary = plot_binary_map(x, ...),
    mask = plot_mask_prop(x, ...)
  )
}

#' Plot HRF design matrix for a single subject
#'
#' Visualizes the design matrix from a working-HRF fit for one
#' subject. Dispatches to one of two layouts:
#' \describe{
#'   \item{\code{style = "default"}}{Faceted panels via \code{\link{plot_design_unlapped}}:
#'     one ggplot per HRF type (main / dHRF / ddHRF), each faceted by task.}
#'   \item{\code{style = "overlapping"}}{Single-panel overlay via
#'     \code{\link{plot_design}} (\code{method = "lineplot"}, \code{style = "overlapping"}):
#'     all regressors plotted together, with color/linetype per field.}
#' }
#'
#' @param x A \code{workingHRF} object (e.g. \code{combo$fit_workingHRF}).
#' @param subject Integer. Index of the subject whose design matrix should be plotted.
#' @param style \code{"default"} (faceted) or \code{"overlapping"} (single overlay).
#' @param ... Additional arguments forwarded to the underlying plot function
#'   (\code{\link{plot_design_unlapped}} or \code{\link{plot_design}}). Typical
#'   args for \code{"overlapping"}: \code{drop_derivatives}, \code{drop_onset_offset},
#'   \code{colors}, \code{linetype}.
#'
#' @return Whatever the dispatched plot function returns (invisibly).
#'
#' @import ggplot2
#' @keywords internal
plot_design_fit <- function(x, subject = 1, style = c("default", "overlapping"), ...) {
  style <- match.arg(style)
  design_dHRFs <- x[["subject_results"]][[subject]][["design_matrix"]]
  if (style == "overlapping") {
    plot_design(design_dHRFs, method = "lineplot", style = "overlapping", ...)
  } else {
    plot_design_unlapped(design_dHRFs, ...)
  }
}

#' Plot cross-subject activation proportion map
#'
#' Displays a surface plot of the proportion of subjects showing significant activation
#' at each location, based on F-test p-values. The surface is generated using the
#' first successful subject's geometry as a template.
#'
#' @param x A \code{workingHRF} object.
#' @param alpha Optional numeric. Significance threshold for p-values (e.g., 0.001). 
#'   If not provided, uses the \code{alpha} stored in \code{x$activation_masks}.
#' @param colors Color map for the surface (default: \code{"viridis"}).
#' @param title Optional plot title. Auto-generated if not provided.
#' @param fname Optional filename to save the plot (e.g., \code{"plot.png"}).
#' @param width,height Pixel dimensions for saved plot (default: 1200x800).
#' @param shadows Shadow depth for surface rendering (default: 0.3).
#' @param material List of lighting/material properties (e.g., \code{list(lit = TRUE)}).
#' @param ... Additional arguments passed to \code{\link[ciftiTools]{view_xifti_surface}}.
#'
#' @importFrom ciftiTools view_xifti_surface
#'
#' @return Invisibly returns a \code{xifti} object containing the activation proportions.
#'         The surface plot is rendered or saved as specified.
#'
#' @keywords internal
plot_activation_proportion <- function(x, alpha = NULL, colors = 'viridis',
                                                  title = NULL, fname = NULL,
                                                  width = 1200, height = 800,
                                                  shadows = 1,
                                                  material = list(lit = TRUE, smooth = FALSE),
                                                  ...) {

  # Use alpha from object if not provided
  if(is.null(alpha)) {
    alpha <- x$activation_masks$alpha
    proportions <- x$activation_masks$prop
  } else {
    # Recalculate proportions with new alpha
    proportions <- recalculate_proportions(x, alpha)
  }

  # Get template from first successful subject
  successful_subjects <- which(sapply(x$subject_results, function(s) s$status == "success"))
  if(length(successful_subjects) == 0) {
    stop("No successful subjects found")
  }

  xii_template <- x$subject_results[[successful_subjects[1]]]$glm_results$pvalF_xii

  # Create title if not provided
  if(is.null(title)) {
    title <- paste0('% Subjects Active (F-test p < ', alpha, ')')
  }

  # Create the xifti object with new data
  result_xifti <- ciftiTools::newdata_xifti(xii_template, proportions)

  # Use view_xifti_surface with built-in saving
  plot_result <- view_xifti_surface(
    xifti = result_xifti,
    zlim = c(0, 1),
    colors = colors,
    title = title,
    fname = fname,
    width = width,
    height = height,
    shadows = shadows,      # Control shadows like in the example
    material = material,    # Control material properties
    bg = "white",           # Background color
    NA_color="#505560",  
    ...
  )
  return(invisible(result_xifti))
}

#' Plot binary activation mask across subjects
#'
#' Displays a surface plot of the binary activation mask based on a user-specified
#' proportion threshold and significance level. Each brain location is marked as active (1)
#' if the proportion of subjects with F-test p-values below `alpha` exceeds `threshold`.
#'
#' The surface is rendered using the geometry from the first successful subject.
#'
#' @param x A \code{workingHRF} object.
#' @param threshold Numeric. Proportion threshold for declaring activation (default: \code{0.1}).
#' @param alpha Optional numeric. Significance threshold for F-test p-values (e.g., \code{0.001}). 
#'   If not provided, uses the \code{alpha} stored in \code{x$activation_masks}.
#' @param colors Color map for the surface (default: \code{"viridis"}).
#' @param title Optional plot title. Auto-generated if not provided.
#' @param fname Optional filename to save the plot (e.g., \code{"plot.png"}).
#' @param width,height Pixel dimensions for saved plot (default: \code{1200x800}).
#' @param shadows Shadow depth for surface rendering (default: \code{1}).
#' @param material List of lighting/material properties for rendering
#'   (e.g., \code{list(lit = TRUE, smooth = FALSE)}).
#' @param ... Additional arguments passed to \code{\link[ciftiTools]{view_xifti_surface}}.
#'
#' @importFrom ciftiTools view_xifti_surface
#'
#' @return Invisibly returns a \code{xifti} object containing the binary activation mask.
#'         The surface plot is rendered or saved as specified.
#'
#' @keywords internal
plot_binary_map <-  function(x, threshold = NULL, alpha = NULL,
                                                 colors = 'viridis', title = NULL,
                                                 fname = NULL,
                                                 width = 1200, height = 800,
                                                 shadows = 1,
                                                 material = list(lit = TRUE, smooth = FALSE),
                                                 ...) {

  if(is.null(threshold)) {
    min_active <- x$call_info$min_active_subjects
    if(!is.null(min_active)) {
      total_subjects <- sum(vapply(x[["subject_results"]], \(s) identical(s[["status"]], "success"), logical(1)))
      threshold <- round(min_active / total_subjects, 2)
    } else {
      # Default for old results - use 10% or minimum 2 subjects, whichever is larger
      total_subjects <- sum(vapply(x[["subject_results"]], \(s) identical(s[["status"]], "success"), logical(1)))
      min_active <- max(2, ceiling(total_subjects * 0.1))
      threshold <- round(min_active / total_subjects, 2)
      cat("Using threshold:", threshold, "\n")
    }
  }
  # Use alpha from object if not provided
  if(is.null(alpha)) {
    alpha <- x$activation_masks$alpha
    proportions <- x$activation_masks$prop
  } else {
    # Recalculate proportions with new alpha
    proportions <- recalculate_proportions(x, alpha)
  }

  # Get template from first successful subject
  successful_subjects <- which(sapply(x$subject_results, function(s) s$status == "success"))
  if(length(successful_subjects) == 0) {
    stop("No successful subjects found")
  }

  xii_template <- x$subject_results[[successful_subjects[1]]]$glm_results$pvalF_xii

  # Create binary threshold map
  binary_map <- as.numeric(proportions > threshold)

  # Create title if not provided
  if(is.null(title)) {
    threshold_pct <- threshold * 100
    title <- paste0('>', threshold_pct, '% Subjects Active (F-test p < ', alpha, ')')
  }

  # Create the xifti object with binary data
  result_xifti <- ciftiTools::newdata_xifti(xii_template, binary_map)

  # Use view_xifti_surface with built-in saving (same as activation_proportion)
  plot_result <- view_xifti_surface(
    xifti = result_xifti,
    zlim = c(0, 1),
    colors = colors,
    title = title,
    fname = fname,
    width = width,
    height = height,
    shadows = shadows,      # Control shadows
    material = material,    # Control material properties
    bg = "white",           # Background color
    NA_color="#505560",  
    ...
  )

  return(invisible(result_xifti))
}


#' Recompute activation proportions across subjects
#'
#' Calculates the proportion of subjects showing significant activation at each brain
#' location, using a user-specified alpha threshold. This is useful when overriding
#' the default threshold stored in a \code{workingHRF} object.
#'
#' @param x A \code{workingHRF} object.
#' @param alpha Numeric significance threshold (e.g., 0.001).
#'
#' @return A numeric vector of length equal to the number of brain locations,
#'         representing the proportion of subjects with p < alpha at each location.
#'
#' @keywords internal
recalculate_proportions <- function(x, alpha) {
  # Extract p-values from all subjects
  n_subjects <- length(x$subject_results)
  n_locations <- x$subject_results[[1]]$n_locations

  pvals_matrix <- matrix(NA, nrow = n_locations, ncol = n_subjects)

  for(i in 1:n_subjects) {
    if(x$subject_results[[i]]$status == "success") {
      pvals <- as.vector(as.matrix(x$subject_results[[i]]$glm_results$pvalF_xii))
      pvals_matrix[, i] <- pvals
    }
  }

  # Create new masks and proportions
  masks <- (pvals_matrix < alpha)
  prop <- rowSums(masks, na.rm = TRUE) / rowSums(!is.na(masks))

  return(prop)
}

#' Plot mask_prop_NA (population activation mask)
#'
#' Visualize the population-level activation mask across subjects,
#' using a provided xifti template for surface geometry.
#'
#' @param x A \code{regularizeHRFs} object containing \code{mask_prop_NA}.
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
                           fname = NULL,
                           title = NULL,
                           legend_fname = NULL,
                           shadows = 1,
                           material = list(lit = TRUE, smooth = FALSE),
                           NA_color="#505560",
                           ...  # shadows, material, etc will pass through
) {
  stopifnot("mask_prop_NA not found in object" = "mask_prop_NA" %in% names(x[["activation_masks"]]))

  # Extract xii from workingHRF_results (which is x)
  xii <- x[["subject_results"]][[1]][["glm_results"]][["bestmodel_xii"]]
  xii <- ciftiTools::convert_xifti(xii, "dscalar")
  
  # Calculate threshold from x (which is workingHRF_results)
  min_active_subjects <- x[["call_info"]][["min_active_subjects"]]
  n_successful <- sum(vapply(x[["subject_results"]], 
                             \(res) identical(res[["status"]], "success"), 
                             logical(1)))
  thr <- round(min_active_subjects / n_successful, 2)
  message("Using group mask threshold: ", thr, " (", min_active_subjects, " out of ", n_successful, " successful subjects)")
  if (is.null(title)) {
    title <- paste0("Subjects with Significant \"Working HRF\" Activation > ", min_active_subjects, " subjects")
  }

  plot_obj <- plot(
    newdata_xifti(xii, x[["activation_masks"]][["mask_prop_NA"]]),
    title = title,
    color_mode = "qualitative",
    fname = fname,
    legend_fname = legend_fname,
    shadows = shadows,
    material = material,
    NA_color=NA_color,
    ...
  )
  invisible(plot_obj)
}
