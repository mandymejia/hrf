#' S3 method: use \code{\link[ciftiTools]{view_xifti}} to plot a \code{"BGLM"} object
#'
#' @param x An object of class "BfMRI_design".
#' @param ... Additional arguments to \code{\link{plot_design}}.
#' @method plot BfMRI_design
#' @export
#'
#' @return Result of the call to \code{\link{plot_design}}
#'
plot.BfMRI_design <- function(x, ...){
  stopifnot(inherits(x, "BfMRI_design"))
  plot_design(x$design, ...)
}

#' Plot design matrix
#'
#' @param design The timepoints by fields design matrix or data.frame.
#' @param method \code{"lineplot"} (default) or \code{"imageplot"}.
#' @param ... Additional arguments to \code{plot_design_line} or
#' \code{plot_design_image}.
#' @return A ggplot
#' @export
#'
plot_design <- function(design, method=c("lineplot", "imageplot"), ...){
  design <- as.data.frame(design)
  method <- match.arg(method, c("lineplot", "imageplot"))

  switch(method,
    lineplot=plot_design_line(design, ...),
    imageplot=plot_design_image(design, ...)
  )
}

#' Plot design with lineplot
#'
#' @rdname plot_design
#' @param design The timepoints by fields design matrix or data.frame.
#' @param colors The name of a ColorBrewer palette (see
#'  RColorBrewer::brewer.pal.info and colorbrewer2.org), the name of a
#'  viridisLite palette, or a character vector of colors. Default:
#'  \code{"Set1"}. A length-\code{nK} vector picks per-field colors.
#' @param linetype,linewidth,alpha Parameters for \code{ggplot2::geom_line}.
#'  Defaults: \code{"solid"} linetype, \code{0.7} linewidth and \code{0.8}
#'  alpha. \code{linetype} can also be a vector of options with length matching
#'  the number of fields in \code{design}.
#' @param style \code{"default"} (themed via \code{theme_bw}, legend right,
#'  grid on) or \code{"overlapping"} (themed via \code{theme_minimal}, legend
#'  at bottom, no grid, no legend title -- for the design-matrix overlap
#'  appendix figures where multiple task regressors are plotted on one panel
#'  using custom color/linetype per field).
#' @param drop_derivatives If \code{TRUE}, drop columns whose names match
#'  \code{_dHRF$} or \code{_ddHRF$} (the HRF derivative columns produced by
#'  \code{make_HRFs}). Default \code{FALSE}.
#' @param drop_onset_offset If \code{TRUE}, drop columns named \code{"onset"}
#'  and \code{"offset"} (the shared stimulus-onset and -offset regressors that
#'  \code{make_design} adds when \code{onset = TRUE} / \code{offset = TRUE}).
#'  Default \code{FALSE}.
#' @return A ggplot
#' @export
#' @importFrom ciftiTools make_color_pal
#'
plot_design_line <- function(design, colors="Set1", linetype="solid", linewidth=.7, alpha=.8,
                              style = c("default", "overlapping"),
                              drop_derivatives = FALSE,
                              drop_onset_offset = FALSE){
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please download the `ggplot2` package.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Please download the `tidyr` package.")
  }
  style <- match.arg(style)

  df <- as.data.frame(design)
  if (drop_derivatives) {
    df <- df[, !grepl("_d+HRF$", colnames(df)), drop = FALSE]
  }
  if (drop_onset_offset) {
    df <- df[, !colnames(df) %in% c("onset", "offset"), drop = FALSE]
  }
  nT <- nrow(df)
  nK <- ncol(df)

  # If `colors` is a vector matching nK, use it directly. Otherwise treat as palette name.
  if (length(colors) == nK && all(grepl("^#|^[A-Za-z]", colors))) {
    color_vec <- colors
  } else {
    color_vec <- suppressWarnings(
      make_color_pal(colors, "qualitative", zlim=nK)$color
    )
  }

  df <- cbind(df, data.frame(idx=seq(nT)))
  df <- tidyr::pivot_longer(df, seq(nK))
  colnames(df)[colnames(df)=="name"] <- "Field"
  df$Field <- factor(df$Field, levels=colnames(design)[colnames(design) %in% unique(df$Field)])

  if (length(linetype) == nK) {
    plt <- ggplot2::ggplot(df, ggplot2::aes_string(x="idx", y="value", col="Field", linetype="Field")) +
      ggplot2::geom_hline(yintercept=0, color="black", linetype="dashed") +
      ggplot2::geom_line(linewidth=linewidth, alpha=alpha) +
      ggplot2::scale_linetype_manual(values=linetype)
  } else {
    plt <- ggplot2::ggplot(df, ggplot2::aes_string(x="idx", y="value", col="Field")) +
      ggplot2::geom_hline(yintercept=0, color="black", linetype="dashed") +
      ggplot2::geom_line(linetype=linetype, linewidth=linewidth, alpha=alpha)
  }

  plt <- plt +
    ggplot2::xlab("Volume") + ggplot2::ylab("Value") +
    ggplot2::scale_color_manual(values=color_vec) +
    ggplot2::scale_x_continuous(limits=c(0-1e-8, nT+1e-8), expand=c(0,0))

  if (style == "overlapping") {
    plt + ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(legend.position  = "bottom",
                     panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.line        = ggplot2::element_line(color = "black"),
                     legend.title     = ggplot2::element_blank())
  } else {
    plt + ggplot2::theme_bw()
  }
}

#' Plot design with one facet-per-task panel per HRF type
#'
#' Renders the design matrix as up to three faceted ggplots -- one each for
#' the main HRF, dHRF, and ddHRF regressors -- where each plot is faceted by
#' task (one row per task). Intended for inspecting individual HRF components
#' in detail. Contrast with \code{plot_design_line(style = "overlapping")},
#' which overlays all regressors on a single panel.
#'
#' Column-name convention: task name with optional \code{_dHRF} or
#' \code{_ddHRF} suffix to flag derivative regressors. Onset/offset columns
#' are included as their own tasks.
#'
#' @rdname plot_design
#' @param design The timepoints by fields design matrix or data.frame.
#' @param title_prefix Title for the main-HRF plot (default \code{"Main HRF"}).
#' @return Invisibly, a named list of \code{ggplot} objects (\code{main},
#'   \code{dhrf}, \code{ddhrf}) for whichever HRF types are present.
#'   Plots are also printed to the current device.
#' @export
#' @import ggplot2
plot_design_unlapped <- function(design, title_prefix = "Main HRF") {
  design <- as.data.frame(design)
  volume <- seq_len(nrow(design))

  hrf_types <- list(main = "", dhrf = "_dHRF", ddhrf = "_ddHRF")

  df_list <- list()
  for (hrf_name in names(hrf_types)) {
    suffix <- hrf_types[[hrf_name]]
    regs <- if (suffix == "") {
      colnames(design)[!grepl("_dHRF|_ddHRF", colnames(design))]
    } else {
      grep(paste0(suffix, "$"), colnames(design), value = TRUE)
    }
    if (length(regs) == 0) next
    df_list[[hrf_name]] <- do.call(rbind, lapply(regs, function(reg) {
      base_task <- sub(paste0(suffix, "$"), "", reg)
      data.frame(Volume = volume, Task = base_task, HRF = hrf_name,
                 Value = design[, reg])
    }))
  }

  if (length(df_list) == 0) return(invisible(list()))

  df_long <- do.call(rbind, df_list)
  df_long$Task <- factor(df_long$Task, levels = unique(df_long$Task))
  yrange <- range(df_long$Value, na.rm = TRUE)

  make_panel <- function(hrf_key, title, ltype) {
    ggplot2::ggplot(df_long[df_long$HRF == hrf_key, ]) +
      ggplot2::geom_line(ggplot2::aes(x = .data$Volume, y = .data$Value,
                                      color = .data$Task), linetype = ltype) +
      ggplot2::ylim(yrange) +
      ggplot2::labs(title = title, x = "Volume", y = "Value") +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::theme(plot.background  = ggplot2::element_rect(fill = "white"),
                     panel.background = ggplot2::element_rect(fill = "white"),
                     legend.position  = "right",
                     strip.placement  = "outside") +
      ggplot2::facet_grid(Task ~ .)
  }

  plots <- list()
  if ("main"  %in% names(df_list)) plots$main  <- make_panel("main",  title_prefix, "solid")
  if ("dhrf"  %in% names(df_list)) plots$dhrf  <- make_panel("dhrf",  "dHRF",       "dashed")
  if ("ddhrf" %in% names(df_list)) plots$ddhrf <- make_panel("ddhrf", "ddHRF",      "twodash")

  for (p in plots) print(p)
  invisible(plots)
}

#' Plot design with imageplot
#'
#' @rdname plot_design
#' @param design The timepoints by fields design matrix or data.frame.
#' @return A ggplot
#' @export
#'
plot_design_image <- function(design){
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Please download the `grDevices` package.")
  }
  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop("Please download the `graphics` package.")
  }

  design <- as.matrix(design)
  nT <- nrow(design)
  nK <- ncol(design)

  graphics::image(
    t(design[seq(nT,1),]), # rev
    xlab="Field", ylab="Volume", col=grDevices::gray.colors(256), axes=FALSE
  )
  k_divs <- if (nK == 1) { 0 } else { seq(0, nK-1)/(nK-1) }
  graphics::axis(1, at = k_divs, labels=colnames(design))
  t_divs <-seq(0, 3)/3
  graphics::axis(2, at = t_divs, labels=rev(round(t_divs*nT))) # rev
}
