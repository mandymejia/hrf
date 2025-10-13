#' Plot method for regularizeHRFs objects
#'
#' Creates diagnostic plots for HRF regularization results.
#' Supports visualization of:
#' \itemize{
#'   \item \code{"mask"} – group-level population mask (\code{mask_prop_NA})
#'   \item \code{"variance"} – residual variance maps for model fits
#'   \item \code{"mean"} – population mean parameter maps (\code{param_mean0})
#' }
#'
#' @param x An object of class \code{"regularizeHRFs"} from
#'   \code{\link{regularize_allHRFs}}.
#' @param type Character. Type of plot to generate. Options:
#'   \itemize{
#'     \item \code{"mask"}
#'     \item \code{"variance"}
#'     \item \code{"mean"} – population mean over all subjects
#'     \item \code{"mean_filtered"} – population mean over subjects with activation only
#'   }
#' @param ... Additional arguments passed to internal plotting functions
#'   (e.g., \code{param}, \code{model}, \code{method}, \code{fname}).
#'
#' @return Invisibly returns the result of the plotting call.
#' @export
plot.regularizeHRFs <- function(x, type = c("mask", "variance", "mean", "mean_filtered"), ...) {
  type <- match.arg(type)

  switch(
    type,
    mask          = plot_mask_prop(x, ...),
    variance      = plot_variance(x, ...),
    mean          = plot_mean_param(x, ...),
    mean_filtered = plot_mean_param_filtered(x, ...)
  )
}

#' Plot Population Mean Parameter Maps
#'
#' Creates brain surface visualizations of population-averaged HRF parameter
#' estimates, averaging across all subjects and applying the population-level mask.
#'
#' @param x A \code{regularizeHRFs} object from \code{\link{regularize_allHRFs}}.
#' @param param Character string specifying which parameter to plot ("a1", "b1", or "c").
#' @param fname Optional file path to save the plot. If \code{NULL} (default),
#'   plot is displayed but not saved.
#' @param title Optional custom title for the plot. If \code{NULL}, a default
#'   title is generated.
#' @param shadows Numeric value controlling shadow intensity (default = 1).
#' @param material List with plotting material properties (default uses lit surfaces).
#' @param NA_color Color to use for NA values (default = "#505560").
#' @param ... Additional arguments passed to the plot function.
#'
#' @return Invisibly returns the plot result.
#' @keywords internal
#' @importFrom ciftiTools newdata_xifti
plot_mean_param <- function(x,
                            param,
                            fname = NULL,
                            title = NULL,
                            shadows = 1,
                            material = list(lit = TRUE, smooth = FALSE),
                            NA_color = "#505560",
                            ...) {
  
  # Get the best_params_df from your results
  best_params_df <- x$best_params_df
  mask_prop_NA <- x$mask_prop_NA
  
  # Recreate the "ALL subjects" average (no subject activation filter)
  best_params_df_avg <- best_params_df %>%
    filter(!is.na(mask)) %>%  # Only exclude completely missing data
    group_by(voxel) %>%
    summarize(param_mean = mean(.data[[param]], na.rm = TRUE))
  
  # Apply the population mask
  full_vector <- rep(NA, length(mask_prop_NA))
  full_vector[best_params_df_avg$voxel] <- best_params_df_avg$param_mean * mask_prop_NA[best_params_df_avg$voxel]
  
  # Geometry template
  xii <- attr(x, "xii")
  
  # Visualization parameters
  param_range <- list(
    a1 = c(4, 8),
    b1 = c(0.5, 1.5),
    c  = c(0, 0.17)
  )
  
  zlim_values <- param_range[[param]]
  color_mode <- if (param == "c") "sequential" else "diverging"
  
  if (is.null(title)) {
    title <- paste0("Mean of ", param, " (over all subjects)")
  }
  
  # Create xifti with data
  new_xifti <- ciftiTools::newdata_xifti(xii, full_vector)
  
  plot_result <- plot(
    new_xifti,
    zlim = zlim_values,
    color_mode = color_mode,
    fname = fname,
    title = title,
    shadows = shadows,
    material = material,
    NA_color = NA_color,
    ...
  )
  
  invisible(plot_result)
}

#' Plot Filtered Population Mean Parameter Maps
#'
#' Creates brain surface visualizations of population-averaged HRF parameter
#' estimates, using only subjects with activation at each voxel.
#'
#' @param x A \code{regularizeHRFs} object from \code{\link{regularize_allHRFs}}.
#' @param param Character string specifying which parameter to plot ("a1", "b1", or "c").
#' @param fname Optional file path to save the plot. If \code{NULL} (default),
#'   plot is displayed but not saved.
#' @param title Optional custom title for the plot. If \code{NULL}, a default
#'   title is generated.
#' @param shadows Numeric value controlling shadow intensity (default = 1).
#' @param material List with plotting material properties (default uses lit surfaces).
#' @param NA_color Color to use for NA values (default = "#505560").
#' @param ... Additional arguments passed to the plot function.
#'
#' @return Invisibly returns the plot result.
#' @keywords internal
#' @importFrom ciftiTools newdata_xifti
plot_mean_param_filtered <- function(x,
                                     param,
                                     fname = NULL,
                                     title = NULL,
                                     shadows = 1,
                                     material = list(lit = TRUE, smooth = FALSE),
                                     NA_color = "#505560",
                                     ...) {
  
  # Get the best_params_df from your results
  best_params_df <- x$best_params_df
  mask_prop_NA <- x$mask_prop_NA
  
  # Recreate the "subjects with activation" average
  best_params_df_avg <- best_params_df %>%
    filter(mask == TRUE) %>%  # Only subjects with activation
    group_by(voxel) %>%
    summarize(param_mean = mean(.data[[param]], na.rm = TRUE))
  
  # Apply the population mask
  full_vector <- rep(NA, length(mask_prop_NA))
  full_vector[best_params_df_avg$voxel] <- best_params_df_avg$param_mean * mask_prop_NA[best_params_df_avg$voxel]
  
  # Geometry template
  xii <- attr(x, "xii")
  
  # Visualization parameters
  param_range <- list(
    a1 = c(4, 8),
    b1 = c(0.5, 1.5),
    c  = c(0, 0.17)
  )
  
  zlim_values <- param_range[[param]]
  color_mode <- if (param == "c") "sequential" else "diverging"
  
  if (is.null(title)) {
    title <- paste0("Mean of ", param, " (over subjects with activation)")
  }
  
  # Create xifti with filtered data
  new_xifti <- ciftiTools::newdata_xifti(xii, full_vector)
  
  plot_result <- plot(
    new_xifti,
    zlim = zlim_values,
    color_mode = color_mode,
    fname = fname,
    title = title,
    shadows = shadows,
    material = material,
    NA_color = NA_color,
    ...
  )
  
  invisible(plot_result)
}

#' Plot Residual Variance Maps
#'
#' Creates brain surface plots showing the residual variance from OLS or WLS
#' regression models used in HRF regularization.
#'
#' @param x A \code{regularizeHRFs} object from \code{\link{regularize_allHRFs}}.
#' @param param Character string specifying which parameter to plot ("a1", "b1", or "c").
#' @param model Character string specifying which model to plot ("A" or "B").
#' @param method Character string specifying estimation method ("OLS" or "WLS").
#' @param fname Optional file path to save the plot.
#' @param title Optional custom title for the plot. If \code{NULL}, a default
#'   title is generated.
#' @param legend_fname Optional file path to save the color legend separately.
#' @param shadows Numeric value controlling shadow intensity (default = 1).
#' @param material List with plotting material properties (default uses lit surfaces).
#' @param NA_color Color to use for NA values (default = "#505560").
#' @param ... Additional arguments passed to the plot function.
#'
#' @return Invisibly returns the plot result.
#' @keywords internal
plot_variance <- function(x,
                          param,
                          model = c("A", "B"),
                          method = c("OLS", "WLS"),
                          fname,
                          title = NULL,
                          legend_fname = NULL,
                          shadows = 1,
                          material = list(lit = TRUE, smooth = FALSE),
                          NA_color="#505560", ...) {
  model  <- match.arg(model)
  method <- match.arg(method)

  # Pull the correct sublist
  if (!"regularized_params" %in% names(x)) {
    stop("Object does not contain $regularized_params")
  }
  params_list <- x$regularized_params

  if (!param %in% names(params_list)) {
    stop("No parameter named ", param, " in x$regularized_params")
  }

  # Build key for residual variance field
  field_name <- paste0("residual_variance_xii_", model, "_", method)

  # Extract the right xifti object
  xii_obj <- params_list[[param]][[field_name]]

  if (is.null(xii_obj)) {
    stop("No variance data found for param=", param, ", model=", model, ", method=", method)
  }

  # Get variance ranges for both models to set consistent zlim
  other_model <- if (model == "A") "B" else "A"
  field_name_other <- paste0("residual_variance_xii_", other_model, "_", method)

  xii_other <- params_list[[param]][[field_name_other]]
  var_range_self  <- range(unlist(xii_obj$data), na.rm = TRUE)
  var_range_other <- range(unlist(xii_other$data), na.rm = TRUE)
  zlim_values <- range(c(var_range_self, var_range_other))

  if (is.null(title)) {
    title <- paste0("Variance of Model ", model,
                    " Predictions for param ", param,
                    " using ", method)
  }

  # Plot
  plot_obj <- plot(
    xii_obj,
    title = title,
    color_mode = "sequential",
    zlim = zlim_values,
    fname = fname,
    legend_fname = legend_fname,
    shadows = 1,
    material = list(lit = TRUE, smooth = FALSE),
    NA_color="#505560",
    ...
  )

  invisible(plot_obj)
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
