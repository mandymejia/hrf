#' Plot method for `workingHRF` objects
#'
#' Dispatches to one of three plot types for visualizing results from a
#' \code{\link{fit_workingHRF}} object:
#'
#' \describe{
#'   \item{\code{type = "design"}}{Plots the design matrix for a specific subject.}
#'   \item{\code{type = "proportion"}}{Displays the cross-subject activation proportion map.}
#'   \item{\code{type = "binary"}}{Displays the binary activation mask after thresholding.}
#' }
#'
#' @param x A \code{workingHRF} object.
#' @param type Type of plot to display. One of \code{"design"}, \code{"proportion"}, or \code{"binary"}.
#' @param subject Integer. Subject index for the design matrix plot (only used if \code{type = "design"}).
#' @param ... Additional arguments passed to the specific plot function. If you wish to use you must set all positional arguments.
#'
#' @return A plot is rendered; invisibly returns the underlying object used for plotting.
#'
#' @examples
#' \dontrun{
#' result <- fit_workingHRF(...)
#' plot(result, type = "design", subject = 1)
#' plot(result, type = "proportion")
#' plot(result, type = "binary", threshold = 0.2)
#' }
#'
#' @export
plot.workingHRF <- function(x, type = c("design", "proportion", "binary"), subject = 1, ...) {
  type <- match.arg(type)
  cat("Plotting type:", type, "\n")

  switch(type,
    design = plot_design_fit(x, subject = subject, ...),
    proportion = plot_activation_proportion(x, ...),
    binary = plot_binary_map(x, ...)
  )
}

#' Plot HRF design matrix for a single subject
#'
#' This function generates line plots of the HRF regressors for a single subject's
#' design matrix, as estimated by \code{\link{fit_workingHRF}}. It distinguishes between:
#' \itemize{
#'   \item \strong{Main HRF}: Canonical regressors (no suffix)
#'   \item \strong{dHRF}: Temporal derivatives (columns ending in \code{_dHRF})
#'   \item \strong{ddHRF}: Dispersion derivatives (columns ending in \code{_ddHRF})
#' }
#' Each type is plotted separately using \pkg{ggplot2}, with lines grouped by task.
#'
#' @param x A \code{workingHRF} object as returned by \code{\link{fit_workingHRF}}.
#' @param subject Integer. Index of the subject whose design matrix should be plotted.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named list of \code{ggplot2} objects for each HRF type that was found.
#'         Also prints the plots to the graphics device. Components may include:
#'         \code{main}, \code{dhrf}, \code{ddhrf}.
#'
#' @import ggplot2
#' @keywords internal
plot_design_fit <- function(x, subject = 1, ...) {
  design_dHRFs <- x[["subject_results"]][[subject]][["design_matrix"]]
  plot_design_fit_core(design_dHRFs, ...)
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

#' Core Plotting Function for HRF Design Matrix
#'
#' Internal function to visualize the design matrix containing main HRF, dHRF, and ddHRF
#' regressors across timepoints and tasks. This function is typically used by higher-level
#' preview or diagnostic plotting tools, and is not intended for direct use by end-users.
#'
#' @param design_dHRFs A data frame or matrix containing the full design matrix, including
#'   HRF regressors. Column names must follow the convention: task names, optionally
#'   suffixed with \code{"_dHRF"} or \code{"_ddHRF"}.
#' @param title_prefix Title for the plot (default: \code{"Main HRF"}).
#' @param ... Additional arguments passed to internal \code{\link[ggplot2]{ggplot}} calls
#'   such as \code{theme()} or \code{facet_grid()} for further plot customization.
#'
#' @return Invisibly returns a named list of \code{ggplot2} objects: \code{main},
#'   \code{dhrf}, and \code{ddhrf}, if those regressors are present in the design.
#'   Plots are also rendered to the current device.
#'
#' @import ggplot2
#' @keywords internal
plot_design_fit_core <- function(design_dHRFs, title_prefix = "Main HRF", ...) {
  volume <- 1:nrow(design_dHRFs)

  # Define which HRF types to include and their suffixes
  hrf_types <- list(
    main = "",           # main regressors: no suffix
    dhrf = "_dHRF",      # dhrf regressors: ends with _dHRF
    ddhrf = "_ddHRF"
  )

  # Collect long-format rows for all HRF types
  df_list <- list()

  for (hrf_name in names(hrf_types)) {
    suffix <- hrf_types[[hrf_name]]

    # Get columns that match the suffix ("" means exclude _dHRF/_ddHRF)
    if (suffix == "") {
      regs <- colnames(design_dHRFs)[!grepl("_dHRF|_ddHRF", colnames(design_dHRFs))]
    } else {
      pattern <- paste0(suffix, "$")
      regs <- grep(pattern, colnames(design_dHRFs), value = TRUE)
    }

    df_list[[hrf_name]] <- do.call(rbind, lapply(regs, function(reg) {
      base_task <- sub(paste0(suffix, "$"), "", reg)  # Strip suffix if present
      data.frame(
        Volume = volume,
        Task = base_task,
        HRF = hrf_name,
        Value = design_dHRFs[, reg]
      )
    }))
  }

  df_long <- do.call(rbind, df_list)
  df_long$Task <- factor(df_long$Task, levels = unique(df_long$Task))

  plot_main <- ggplot(data = df_long[df_long$HRF == "main",]) +
    geom_line(aes(x = .data$Volume, y = .data$Value, color = .data$Task), linetype = "solid") +
    ylim(range(df_long$Value, na.rm = TRUE)) +
    labs(title = title_prefix, x = "Volume", y = "Value") +
    theme_minimal(base_size = 15) +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white")) +
    theme(legend.position = "right",
          strip.placement = "outside") +
    list(...) + 
    facet_grid(Task ~ .)

  plot_dhrf <- ggplot(data = df_long[df_long$HRF == "dhrf",]) +
    geom_line(aes(x = .data$Volume, y = .data$Value, color = .data$Task), linetype = "dashed") +
    ylim(range(df_long$Value, na.rm = TRUE)) +
    labs(title = "dHRF", x = "Volume", y = "Value") +
    theme_minimal(base_size = 15) +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white")) +
    theme(legend.position = "right",
          strip.placement = "outside") +
    list(...) +
    facet_grid(Task ~ .)

  plot_ddhrf <- ggplot(data = df_long[df_long$HRF == "ddhrf",]) +
    geom_line(aes(x = .data$Volume, y = .data$Value, color = .data$Task), linetype = "twodash") +
    ylim(range(df_long$Value, na.rm = TRUE)) +
    labs(title = "ddHRF", x = "Volume", y = "Value") +
    theme_minimal(base_size = 15) +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white")) +
    theme(legend.position = "right",
          strip.placement = "outside") +
    list(...) +
    facet_grid(Task ~ .)

  if ("main" %in% names(df_list)) print(plot_main)
  if ("dhrf" %in% names(df_list)) print(plot_dhrf)
  if ("ddhrf" %in% names(df_list)) print(plot_ddhrf)

  invisible(Filter(Negate(is.null), list(
    main = if ("main" %in% names(df_list)) plot_main else NULL,
    dhrf = if ("dhrf" %in% names(df_list)) plot_dhrf else NULL,
    ddhrf = if ("ddhrf" %in% names(df_list)) plot_ddhrf else NULL
  )))
}