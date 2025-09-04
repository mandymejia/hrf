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
#' @param ... Additional arguments passed to the specific plot function.
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
  cat("Plotting type:", type, "\n")
  type <- match.arg(type)

  switch(type,
    design = plot_design_fit(x, subject = subject, ...),
    proportion = plot_activation_proportion(x, ...),
    binary = plot_binary_mask(x, ...)
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
#' @examples
#' \dontrun{
#' result <- fit_workingHRF(...)
#' plots <- plot_design(result, subject = 2)
#' plots$main   # Access the main HRF plot
#' }
#'
#' @import ggplot2
#' @keywords internal
plot_design_fit <- function(x, subject = 1, ...) {
  design_dHRFs <- x[["subject_results"]][[subject]][["design_matrix"]]
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
    labs(title = "Main HRF", x = "Volume", y = "Value") +
    theme_minimal(base_size = 15) +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white")) +
    theme(legend.position = "right",
          strip.placement = "outside") +
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

#' @keywords internal
plot_activation_proportion <- function(x, ...) {
  stop("plot_activation_proportion() is not implemented yet.")
}

#' @keywords internal
plot_binary_mask <- function(x, ...) {
  stop("plot_binary_mask() is not implemented yet.")
}