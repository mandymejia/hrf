#' Plot method for regularizeHRFs objects
#'
#' Diagnostic plots for the offset-based regularize output.
#' \itemize{
#'   \item \code{"pop_avg"} – population-average HRF parameter map (brain
#'     surface, grid-snapped values from \code{x$pop_avg}).
#'   \item \code{"mean_all"} – mean of raw per-subject best a1/b1/c across ALL
#'     subjects (no activation filter), from \code{x$best_params_df}.
#'   \item \code{"mean"} – mean of raw per-subject best a1/b1/c, restricted to
#'     subjects activated at each voxel (\code{best_params_df$mask == TRUE}).
#'   \item \code{"param_heatmap"} – frequency heatmap of per-subject best HRF
#'     parameters on the underlying HRF grid.
#' }
#'
#' @param x An object of class \code{"regularizeHRFs"} from
#'   \code{\link{regularize_allHRFs}}.
#' @param type Character. Plot type to generate.
#' @param ... Additional arguments passed to the type-specific helper
#'   (e.g. \code{param}, \code{mask}, \code{fname}).
#'
#' @return Invisibly returns the plot object.
#' @export
plot.regularizeHRFs <- function(x, type = c("pop_avg", "mean_all", "mean", "param_heatmap"), ...) {
  type <- match.arg(type)

  switch(type,
    pop_avg       = plot_pop_avg(x, ...),
    mean_all      = plot_mean_param(x, ...),
    mean          = plot_mean_param_filtered(x, ...),
    param_heatmap = plot_param_heatmap(x, ...)
  )
}


#' Plot Population Mean Parameter Maps (raw, all subjects)
#'
#' Brain-surface mean of raw per-subject best a1/b1/c from
#' \code{x$best_params_df}, averaged over all subjects with non-\code{NA} mask.
#'
#' @param x A \code{regularizeHRFs} object.
#' @param param One of \code{"a1"}, \code{"b1"}, \code{"c"}.
#' @param fname Optional path to save the rendered PNG.
#' @param title Optional plot title.
#' @param shadows,material,NA_color Standard \code{plot.xifti} aesthetics.
#' @param ... Additional args passed to \code{plot.xifti}.
#'
#' @return Invisibly returns the plot result.
#' @keywords internal
#' @importFrom ciftiTools newdata_xifti
#' @importFrom dplyr filter group_by summarize
plot_mean_param <- function(x,
                            param = c("a1", "b1", "c"),
                            fname = NULL,
                            title = NULL,
                            shadows = 1,
                            material = list(lit = TRUE, smooth = FALSE),
                            NA_color = "#505560",
                            ...) {
  param <- match.arg(param)

  xii_template <- attr(x, "xii")
  if (is.null(xii_template)) {
    stop("attr(x, \"xii\") is missing. Re-run regularize_allHRFs to attach it.")
  }

  best_params_df <- x$best_params_df
  mask_prop_NA   <- x$mask_prop_NA

  bp_avg <- best_params_df %>%
    dplyr::filter(!is.na(.data$mask)) %>%
    dplyr::group_by(.data$voxel) %>%
    dplyr::summarize(param_mean = mean(.data[[param]], na.rm = TRUE), .groups = "drop")

  full_vector <- rep(NA, length(mask_prop_NA))
  full_vector[bp_avg$voxel] <- bp_avg$param_mean * mask_prop_NA[bp_avg$voxel]

  zlim       <- list(a1 = c(4, 8), b1 = c(0.5, 1.5), c = c(0, 0.17))[[param]]
  color_mode <- if (param == "c") "sequential" else "diverging"

  if (is.null(title)) {
    title <- paste0("Mean of ", param, " (over all subjects)")
  }

  new_xifti <- ciftiTools::newdata_xifti(xii_template, full_vector)
  plot_result <- plot(
    new_xifti,
    zlim = zlim,
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


#' Plot Population Mean Parameter Maps (raw, activated subjects only)
#'
#' Brain-surface mean of raw per-subject best a1/b1/c from
#' \code{x$best_params_df}, restricted to subjects with
#' \code{mask == TRUE} at each voxel.
#'
#' @inheritParams plot_mean_param
#' @return Invisibly returns the plot result.
#' @keywords internal
#' @importFrom ciftiTools newdata_xifti
#' @importFrom dplyr filter group_by summarize
plot_mean_param_filtered <- function(x,
                                     param = c("a1", "b1", "c"),
                                     fname = NULL,
                                     title = NULL,
                                     shadows = 1,
                                     material = list(lit = TRUE, smooth = FALSE),
                                     NA_color = "#505560",
                                     ...) {
  param <- match.arg(param)

  xii_template <- attr(x, "xii")
  if (is.null(xii_template)) {
    stop("attr(x, \"xii\") is missing. Re-run regularize_allHRFs to attach it.")
  }

  best_params_df <- x$best_params_df
  mask_prop_NA   <- x$mask_prop_NA

  bp_avg <- best_params_df %>%
    dplyr::filter(.data$mask == TRUE) %>%
    dplyr::group_by(.data$voxel) %>%
    dplyr::summarize(param_mean = mean(.data[[param]], na.rm = TRUE), .groups = "drop")

  full_vector <- rep(NA, length(mask_prop_NA))
  full_vector[bp_avg$voxel] <- bp_avg$param_mean * mask_prop_NA[bp_avg$voxel]

  zlim       <- list(a1 = c(4, 8), b1 = c(0.5, 1.5), c = c(0, 0.17))[[param]]
  color_mode <- if (param == "c") "sequential" else "diverging"

  if (is.null(title)) {
    title <- paste0("Mean of ", param, " (over subjects with activation)")
  }

  new_xifti <- ciftiTools::newdata_xifti(xii_template, full_vector)
  plot_result <- plot(
    new_xifti,
    zlim = zlim,
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


#' Plot Population-Average HRF Parameter Map (Brain Surface)
#'
#' Renders the population-averaged HRF parameter from \code{x$pop_avg}
#' (snapped to the HRF grid) onto the cortex-surface xifti template attached
#' as \code{attr(x, "xii")} by \code{regularize_allHRFs()}.
#'
#' @param x A \code{regularizeHRFs} object.
#' @param param One of \code{"a1"}, \code{"b1"}, \code{"c"}.
#' @param fname Optional path to save the rendered PNG (or HTML on Mac without
#'   working OpenGL).
#' @param title Optional plot title. Defaults to \code{"Population avg <param> (snapped)"}.
#' @param shadows,material,NA_color Standard \code{plot.xifti} aesthetics.
#' @param zlim Optional length-2 numeric. If \code{NULL}, sensible per-param
#'   defaults are used (\code{a1: c(4,8)}, \code{b1: c(0.5,1.5)}, \code{c: c(0,0.17)}).
#' @param ... Additional args passed to \code{plot.xifti}.
#'
#' @return Invisibly returns the plot result.
#' @keywords internal
#' @importFrom ciftiTools newdata_xifti
plot_pop_avg <- function(x,
                         param = c("a1", "b1", "c"),
                         fname = NULL,
                         title = NULL,
                         shadows = 1,
                         material = list(lit = TRUE, smooth = FALSE),
                         NA_color = "#505560",
                         zlim = NULL,
                         ...) {
  param <- match.arg(param)

  xii_template <- attr(x, "xii")
  if (is.null(xii_template)) {
    stop("attr(x, \"xii\") is missing. This regularizeHRFs result was likely built ",
         "before the xii template was attached. Re-run regularize_allHRFs().")
  }

  pop_avg <- x$pop_avg
  mask_prop_NA <- x$mask_prop_NA
  snapped_col <- paste0(param, "_snapped")
  if (!snapped_col %in% colnames(pop_avg)) {
    stop("Column '", snapped_col, "' not found in x$pop_avg.")
  }

  full_vector <- rep(NA, length(mask_prop_NA))
  full_vector[pop_avg$voxel] <- pop_avg[[snapped_col]] * mask_prop_NA[pop_avg$voxel]

  if (is.null(zlim)) {
    zlim <- list(a1 = c(4, 8), b1 = c(0.5, 1.5), c = c(0, 0.17))[[param]]
  }
  if (is.null(title)) {
    title <- paste0("Population avg ", param, " (snapped)")
  }

  new_xifti <- ciftiTools::newdata_xifti(xii_template, full_vector)
  plot_result <- plot(
    new_xifti,
    zlim = zlim,
    color_mode = if (param == "c") "sequential" else "diverging",
    fname = fname,
    title = title,
    shadows = shadows,
    material = material,
    NA_color = NA_color,
    ...
  )
  invisible(plot_result)
}


#' Plot best-per-subject HRF parameter frequency on the HRF grid
#'
#' Heatmap of how often each \code{(a1, b1, c)} grid combination is the
#' per-voxel per-subject best HRF, faceted by \code{c}.
#'
#' @param x A \code{regularizeHRFs} object.
#' @param mask Logical. If \code{TRUE}, restrict to voxel-subject pairs where
#'   the working-HRF activation mask is \code{TRUE}.
#' @param title Optional custom title.
#' @param ... Additional ggplot \code{theme()} args (passed through via
#'   \code{list(...)}).
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom dplyr distinct summarise group_by select left_join filter mutate
#' @keywords internal
plot_param_heatmap <- function(x, mask = FALSE, title = NULL, ...) {
  if (!inherits(x, "regularizeHRFs")) {
    stop("Input object must be of class 'regularizeHRFs'.")
  }

  best_params_df <- x$best_params_df
  if (mask) {
    best_params_df <- best_params_df %>% dplyr::filter(.data$mask == TRUE)
  }

  if (is.null(title)) {
    title <- if (mask) "Best per-subject HRF params (masked)" else "Best per-subject HRF params"
  }

  hrf_grid <- x$hrf_grid
  best_params_df <- best_params_df %>% dplyr::mutate(c = round(.data$c, 7))

  complete_grid <- hrf_grid %>%
    dplyr::select(.data$a1, .data$b1, .data$c) %>%
    dplyr::distinct() %>%
    dplyr::mutate(c = round(.data$c, 7))

  freq_df <- best_params_df %>%
    dplyr::group_by(.data$a1, .data$b1, .data$c) %>%
    dplyr::summarise(freq = dplyr::n(), .groups = "drop")

  freq_df <- complete_grid %>%
    dplyr::left_join(freq_df, by = c("a1", "b1", "c")) %>%
    dplyr::mutate(freq = ifelse(is.na(.data$freq), 0, .data$freq)) %>%
    dplyr::mutate(rel_freq = .data$freq / sum(.data$freq))

  freq_df$c_label <- ifelse(freq_df$c == 0,
                            "No Undershoot (c=0)",
                            ifelse(abs(freq_df$c - 1/6) < 1e-5,
                                   "With Undershoot (c=1/6)",
                                   paste0("With Undershoot (c=", freq_df$c, ")")))

  a1_vals <- sort(unique(complete_grid$a1))
  b1_vals <- sort(unique(complete_grid$b1))
  a1_grid <- c(min(a1_vals) - 0.5, a1_vals + 0.5)
  b1_grid <- c(min(b1_vals) - 0.125, b1_vals + 0.125)

  fill_limit <- max(freq_df$rel_freq, na.rm = TRUE)

  ggplot(freq_df, aes(x = .data$a1, y = .data$b1, fill = .data$rel_freq)) +
    geom_hline(yintercept = b1_grid, alpha = 0.1) +
    geom_vline(xintercept = a1_grid, alpha = 0.1) +
    geom_tile(alpha = 0.8) +
    geom_text(aes(label = ifelse(.data$freq == 0 | .data$rel_freq * 100 < 0.1, "0%",
                                 ifelse(.data$rel_freq * 100 < 1,
                                        sprintf("%.1f%%", .data$rel_freq * 100),
                                        paste0(signif(.data$rel_freq * 100, 2), "%"))))) +
    geom_rect(
      data = dplyr::filter(freq_df, .data$a1 == 6, .data$b1 == 1, .data$c != 0),
      aes(xmin = .data$a1 - 0.5, xmax = .data$a1 + 0.5,
          ymin = .data$b1 - 0.125, ymax = .data$b1 + 0.125),
      color = "lightgrey", linewidth = 1.2, fill = NA
    ) +
    scale_fill_viridis_c(
      "Relative \nFrequency", option = "A",
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, fill_limit)
    ) +
    scale_x_continuous(breaks = a1_vals, expand = c(0, 0)) +
    scale_y_continuous(breaks = b1_vals, expand = c(0, 0)) +
    facet_grid(. ~ c_label) +
    theme_few() + list(...) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.spacing = unit(1.5, "lines"))
}


