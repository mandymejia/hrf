#' Plot method for regularizeHRFs objects
#'
#' Creates diagnostic plots for HRF regularization results.
#' Supports visualization of:
#' \itemize{
#'   \item \code{"variance"} - residual variance maps for model fits
#'   \item \code{"precision"} - inverse variance (precision) maps
#'   \item \code{"mean_all"} - population mean parameter maps
#'   \item \code{"mean"} - population mean over subjects with activation only
#'   \item \code{"param_heatmap"} - heatmap of parameter frequency distributions
#' }
#'
#' @param x An object of class \code{"regularizeHRFs"} from
#'   \code{\link{regularize_allHRFs}}.
#' @param type Character. Type of plot to generate. Options:
#'   \itemize{
#'     \item \code{"variance"}
#'     \item \code{"precision"}
#'     \item \code{"mean_all"} - population mean over all subjects
#'     \item \code{"mean"} - population mean over subjects with activation only
#'     \item \code{"param_heatmap"} - frequency heatmap of parameter combinations
#'   }
#' @param ... Additional arguments passed to internal plotting functions
#'   (e.g., \code{param}, \code{model}, \code{method}, \code{fname}).
#'
#' @return Invisibly returns the result of the plotting call.
#' @export
plot.regularizeHRFs <- function(x, type = c("variance", "mean_all", "mean", "precision", "param_heatmap"), ...) {
  type <- match.arg(type)

  switch(
    type,
    variance      = plot_variance(x, ...),
    precision     = plot_precision(x, ...),
    mean_all      = plot_mean_param(x, ...),
    mean          = plot_mean_param_filtered(x, ...),
    param_heatmap = plot_param_heatmap(x, ...)
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
#' Creates brain surface visualizations of the residual variance from OLS or WLS
#' regression models used during HRF parameter regularization. Each map shows
#' the voxel-wise variance of parameter fits across subjects for either the
#' intercept-slope (\code{"IS"}) or intercept-only (\code{"IO"}) regularization fit.
#'
#' @param x A \code{regularizeHRFs} object from \code{\link{regularize_allHRFs}}.
#' @param param Character string specifying which HRF parameter to plot
#'   (\code{"a1"}, \code{"b1"}, or \code{"c"}).
#' @param model Character string specifying which regularization fit to visualize:
#'   \code{"IS"} for the intercept-slope fit, or \code{"IO"} for the intercept-only fit.
#' @param method Character string specifying the estimation method used
#'   (\code{"OLS"} or \code{"WLS"}).
#' @param fname Optional file path to save the plot output. If \code{NULL}
#'   (default), the plot is displayed interactively but not saved.
#' @param title Optional custom plot title. If \code{NULL}, a descriptive title
#'   is generated automatically.
#' @param legend_fname Optional file path to save the color legend separately.
#' @param shadows Numeric value controlling surface shadow intensity (default = 1).
#' @param material List of surface rendering options (default uses
#'   \code{list(lit = TRUE, smooth = FALSE)}).
#' @param NA_color Color to use for NA or masked values (default = "#505560").
#' @param ... Additional arguments passed to the underlying plotting function.
#'
#' @details
#' The function compares voxel-level variance across both regularization fits
#' (\code{"IS"} and \code{"IO"}) and harmonizes color scaling to ensure consistent
#' visual interpretation between the two.
#'
#' Lower variance values indicate more stable parameter estimates across subjects,
#' whereas higher values reflect increased variability or lower reliability in
#' the estimated HRF parameters.
#'
#' @return Invisibly returns the plot object.
#'
#' @keywords internal
plot_variance <- function(x,
                          param,
                          model = c("IS", "IO"),
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
  other_model <- if (model == "IS") "IO" else "IS"
  field_name_other <- paste0("residual_variance_xii_", other_model, "_", method)

  xii_other <- params_list[[param]][[field_name_other]]
  var_range_self  <- range(unlist(xii_obj$data), na.rm = TRUE)
  var_range_other <- range(unlist(xii_other$data), na.rm = TRUE)
  zlim_values <- range(c(var_range_self, var_range_other))
  pretty_model <- if (model == "IS") "intercept-slope" else "intercept-only"

  if (is.null(title)) {
    title <- paste0("Variance of Model ", pretty_model,
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
#' Creates brain surface plots showing the *precision* (inverse residual variance)
#' derived from OLS or WLS regression models used in HRF regularization.
#' These maps highlight areas of higher certainty or reliability in the estimated
#' HRF parameters across the cortical surface.
#'
#' @param x A \code{regularizeHRFs} object from \code{\link{regularize_allHRFs}}.
#' @param param Character string specifying which HRF parameter to plot
#'   (\code{"a1"}, \code{"b1"}, or \code{"c"}).
#' @param model Character string specifying which regularization fit to visualize:
#'   \code{"IS"} for the intercept-slope fit, or \code{"IO"} for the intercept-only fit.
#' @param method Character string specifying estimation method (\code{"OLS"} or \code{"WLS"}).
#' @param fname Optional file path to save the plot. If \code{NULL}, the plot
#'   is displayed interactively but not saved.
#' @param title Optional custom plot title. If \code{NULL}, a descriptive title
#'   is generated automatically.
#' @param legend_fname Optional file path to save the color legend separately.
#' @param shadows Numeric value controlling shadow intensity (default = 1).
#' @param material List specifying surface rendering properties (default:
#'   \code{list(lit = TRUE, smooth = FALSE)}).
#' @param NA_color Color to use for NA or masked values (default = "#505560").
#' @param ... Additional arguments passed to the underlying plotting function.
#'
#' @details
#' Precision maps are computed as:
#' \deqn{Precision = 1 / Variance}
#'
#' Higher precision values indicate more stable and reliable parameter estimates
#' across subjects. The function automatically aligns color scales across the two
#' regularization fits (\code{"IS"} and \code{"IO"}) to allow direct visual comparison.
#'
#' @return Invisibly returns the plot object.
#'
#' @keywords internal
plot_precision <- function(x,
                           param,
                           model = c("IS", "IO"),
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
  other_model <- if (model == "IS") "IO" else "IS"
  xii_other <- params_list[[param]][[paste0("residual_variance_xii_", other_model, "_", method)]]
  zlim_values <- range(c(dat_inv, 1 / as.numeric(as.matrix(xii_other))), na.rm = TRUE)
  pretty_model <- if (model == "IS") "intercept-slope" else "intercept-only"
  if (is.null(title)) {
    title <- paste0("Precision (1 / Variance) of Model ", pretty_model,
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


#' Plot parameter frequency heatmap
#'
#' Generates a heatmap showing the frequency of HRF parameter combinations
#' across the fitted grid for a given model type. Optionally adds a custom
#' plot title; if not provided, an informative default title is generated
#' based on the selected model.
#'
#' @param x A \code{regularizeHRFs} object from
#'   \code{\link{regularize_allHRFs}}.
#' @param model Character. Which model to use for plotting parameter frequencies:
#'   \code{"intercept_only"}, \code{"intercept_slope"}, or \code{"best_params"}.
#' @param title Character or \code{NULL}. Custom plot title. If \code{NULL}
#'   (default), a descriptive title is generated automatically based on
#'   \code{model}.
#' @param ... Any additional arguments to pass into ggplot formation via list(...).
#'
#' @return A \code{ggplot} object showing the parameter frequency heatmap.
#' @details
#' This plot visualizes how often each combination of HRF parameters
#' (\code{a1}, \code{b1}, \code{c}) occurs in the regularized results, with
#' relative frequency encoded by fill intensity. Facets indicate whether
#' the model includes an undershoot component (\code{c = 1/6}) or not
#' (\code{c = 0}).
#'
#' @examples
#' \dontrun{
#' plot_param_heatmap(result, model = "intercept_slope")
#' plot_param_heatmap(result, model = "best_params", title = "Custom Heatmap Title")
#' }
#'
#' @importFrom dplyr mutate distinct summarise group_by select left_join
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_tile geom_text
#'   geom_rect scale_fill_viridis_c scale_x_continuous scale_y_continuous
#'   facet_grid ggtitle
#' @importFrom scales percent_format
#' @importFrom ggthemes theme_few
#'
#' @keywords internal
#' @export
plot_param_heatmap <- function(x,
                               model = c("intercept_only", "intercept_slope", "best_params"),
                               title = NULL,
                               ...) {
  model <- match.arg(model)

  if (!inherits(x, "regularizeHRFs")) {
    stop("Input object must be of class 'regularizeHRFs'.")
  }

  # Define rounded_params and hrf_grid based on model type
  if (model == "intercept_only") {
    rounded_params <- x[["rounded_params"]][["rounded_intercept_only"]]
  } else if (model == "intercept_slope") {
    rounded_params <- x[["rounded_params"]][["rounded_intercept_slope"]]
  } else if (model == "best_params") {
    rounded_params <- x[["best_params_df"]]
  }

  if (is.null(title)) {
    model_label <- switch(model,
      intercept_only = "Intercept-Only Model",
      intercept_slope = "Intercept + Slope Model",
      best_params = "Best Parameter Estimates"
    )
    title <- paste("Parameter Frequency Heatmap -", model_label)
  }

  hrf_grid <- x[["hrf_grid"]]

  rounded_params <- rounded_params %>%
    mutate(c = round(.data$c, 7))

  # Get the complete grid of a1, b1, c combinations from hrf_grid
  complete_grid <- hrf_grid %>%
    select(.data$a1, .data$b1, .data$c) %>%
    distinct() %>%
    mutate(c = round(.data$c, 7))  # Round here too



  # Aggregate frequencies from rounded_params
  freq_df <- rounded_params %>%
    group_by(.data$a1, .data$b1, .data$c) %>%
    summarise(freq = dplyr::n(), .groups = "drop")

  # Left join to get zeros for missing combinations
  freq_df <- complete_grid %>%
    left_join(freq_df, by = c("a1", "b1", "c")) %>%
    mutate(freq = ifelse(is.na(.data$freq), 0, .data$freq))

  # Compute relative frequency for fill (legend only)
  freq_df <- freq_df %>%
    mutate(rel_freq = .data$freq / sum(.data$freq))

  # Create facet labels
  freq_df$c_label <- ifelse(freq_df$c == 0,
                            "No Undershoot (c=0)",
                            "With Undershoot (c=1/6)")

  # Gridline positions - now based on the complete grid
  a1_vals <- sort(unique(complete_grid$a1))
  b1_vals <- sort(unique(complete_grid$b1))
  a1_grid <- c(min(a1_vals) - 0.5, a1_vals + 0.5)
  b1_grid <- c(min(b1_vals) - 0.125, b1_vals + 0.125)

  p <- ggplot(freq_df, aes(x = .data$a1, y = .data$b1, fill = .data$rel_freq)) +
    geom_hline(yintercept = b1_grid, alpha = 0.1) +
    geom_vline(xintercept = a1_grid, alpha = 0.1) +
    geom_tile(alpha = 0.8) +
    geom_text(aes(label = ifelse(.data$freq == 0, "0", .data$freq))) +
    geom_rect(
      data = dplyr::filter(freq_df, .data$a1 == 6, .data$b1 == 1, abs(.data$c - 1/6) < 1e-5) %>%
            dplyr::mutate(c_label = "With Undershoot (c=1/6)"),
      aes(xmin = .data$a1 - 0.5, xmax = .data$a1 + 0.5,
          ymin = .data$b1 - 0.125, ymax = .data$b1 + 0.125),
      color = "black", size = 1.2, fill = NA
    ) +
    scale_fill_viridis_c(
      'Relative \nFrequency', option = 'A',
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, max(freq_df$rel_freq))
    ) +
    scale_x_continuous(breaks = a1_vals, expand = c(0, 0)) +
    scale_y_continuous(breaks = b1_vals, expand = c(0, 0)) +
    facet_grid(. ~ c_label) +
    theme_few() + list(...) +
    ggtitle(title)

  return(p)

 }
