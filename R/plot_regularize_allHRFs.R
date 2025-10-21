#' Plot method for regularizeHRFs objects
#'
#' Creates diagnostic plots for HRF regularization results.
#' Supports visualization of:
#' \itemize{
#'   \item \code{"variance"} – residual variance maps for model fits
#'   \item \code{"mean"} – population mean parameter maps (\code{param_mean0})
#'   \item \code{"mean_filtered"} – population mean over subjects with activation only
#' }
#'
#' @param x An object of class \code{"regularizeHRFs"} from
#'   \code{\link{regularize_allHRFs}}.
#' @param type Character. Type of plot to generate. Options:
#'   \itemize{
#'     \item \code{"variance"}
#'     \item \code{"mean_all"} – population mean over all subjects
#'     \item \code{"mean"} – population mean over subjects with activation only
#'   }
#' @param ... Additional arguments passed to internal plotting functions
#'   (e.g., \code{param}, \code{model}, \code{method}, \code{fname}).
#'
#' @return Invisibly returns the result of the plotting call.
#' @export
plot.regularizeHRFs <- function(x, type = c("variance", "mean_all", "mean", "precision"), ...) {
  type <- match.arg(type)

  switch(
    type,
    variance      = plot_variance(x, ...),
    precision     = plot_precision(x, ...),
    mean_all      = plot_mean_param(x, ...),
    mean          = plot_mean_param_filtered(x, ...)
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

#' Plot Precision (Inverse Variance) Maps
#'
#' Creates brain surface plots showing the *precision* (i.e., inverse residual variance)
#' derived from OLS or WLS regression models used in HRF regularization.
#' This visualization highlights areas of higher certainty in the parameter estimates
#' across the cortical surface.
#'
#' @param x A \code{regularizeHRFs} object from \code{\link{regularize_allHRFs}}.
#' @param param Character string specifying which HRF parameter to plot
#'   ("a1", "b1", or "c").
#' @param model Character string specifying which model to plot ("A" or "B").
#' @param method Character string specifying estimation method ("OLS" or "WLS").
#' @param fname Optional file path to save the plot. If \code{NULL}, the plot
#'   is displayed interactively but not saved.
#' @param title Optional custom plot title. If \code{NULL}, a default title is
#'   generated based on model, parameter, and method.
#' @param legend_fname Optional file path to save the color legend separately.
#' @param shadows Numeric value controlling shadow intensity (default = 1).
#' @param material List specifying surface rendering properties. Defaults to
#'   \code{list(lit = TRUE, smooth = FALSE)} for a lightly shaded surface.
#' @param NA_color Color to use for NA or masked values (default = "#505560").
#' @param ... Additional arguments passed to the underlying
#'   \code{\link[ciftiTools]{plot}} method.
#'
#' @details
#' Precision maps are generated by inverting the residual variance values from
#' the HRF regularization results:
#' \deqn{Precision = 1 / Variance}
#'
#' This transformation provides a spatial view of the confidence or reliability
#' of the parameter estimates — higher precision values indicate more stable,
#' less variable fits across subjects.
#'
#' The function automatically ensures consistent color scaling between models
#' ("A" and "B") to facilitate visual comparison.
#'
#' @return Invisibly returns the plot result.
#'
#' @keywords internal
plot_precision <- function(x,
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

  params_list <- x$regularized_params
  if (!param %in% names(params_list)) stop("No parameter named ", param)

  field_name <- paste0("residual_variance_xii_", model, "_", method)
  xii_obj <- params_list[[param]][[field_name]]
  if (is.null(xii_obj)) stop("No variance data found.")

  # Invert the variance values
  dat <- as.numeric(as.matrix(xii_obj))
  dat_inv <- 1 / dat
  dat_inv[!is.finite(dat_inv)] <- NA
  xii_inv <- ciftiTools::newdata_xifti(xii_obj, dat_inv)

  # Match color scale to other model
  other_model <- if (model == "A") "B" else "A"
  xii_other <- params_list[[param]][[paste0("residual_variance_xii_", other_model, "_", method)]]
  zlim_values <- range(c(dat_inv, 1 / as.numeric(as.matrix(xii_other))), na.rm = TRUE)

  if (is.null(title)) {
    title <- paste0("Precision (1 / Variance) of Model ", model,
                    " Predictions for param ", param,
                    " using ", method)
  }

  plot_obj <- plot(
    xii_inv,
    title = title,
    color_mode = "sequential",
    zlim = zlim_values,
    fname = fname,
    legend_fname = legend_fname,
    shadows = shadows,
    material = material,
    NA_color = NA_color,
    ...
  )

  invisible(plot_obj)
}