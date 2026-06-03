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
plot.regularizeHRFs <- function(x, type = c("pop_avg", "pop_avg_continuous",
                                              "mean_all", "mean", "param_heatmap"), ...) {
  type <- match.arg(type)

  switch(type,
    pop_avg            = plot_pop_avg(x, ...),
    pop_avg_continuous = plot_pop_avg_continuous(x, ...),
    mean_all           = plot_mean_param(x, ...),
    mean               = plot_mean_param_filtered(x, ...),
    param_heatmap      = plot_param_heatmap(x, ...)
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
                            param = c("a1", "b1", "c", "time_to_peak", "FWHM"),
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

  # Defaults tuned for the standard HRF grid (a1 in {3..12}, b1 in {0.5..2},
  # c in {0, 1/6}). Pass `zlim` to override for any other grid.
  if (is.null(zlim)) {
    # Dynamic from the actual HRF grid (people can use different grids).
    zlim <- range(x$hrf_grid[[param]])
  }
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
                                     param = c("a1", "b1", "c", "time_to_peak", "FWHM"),
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

  # Defaults tuned for the standard HRF grid (a1 in {3..12}, b1 in {0.5..2},
  # c in {0, 1/6}). Pass `zlim` to override for any other grid.
  if (is.null(zlim)) {
    # Dynamic from the actual HRF grid (people can use different grids).
    zlim <- range(x$hrf_grid[[param]])
  }
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
#' @param zlim Optional length-2 numeric. If \code{NULL}, dynamically derived
#'   from \code{range(x$hrf_grid[[param]])} so the colorbar always spans the
#'   actual grid the user picked.
#' @param ... Additional args passed to \code{plot.xifti}.
#'
#' @return Invisibly returns the plot result.
#' @keywords internal
#' @importFrom ciftiTools newdata_xifti
plot_pop_avg <- function(x,
                         param = c("a1", "b1", "c", "time_to_peak", "FWHM"),
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
  # a1/b1/c use snapped grid values; t2p/FWHM use the imputed continuous means
  # produced by unmask_pop_avg (median-filled + surf_FWHM smoothed).
  col_name <- switch(param,
    a1           = "a1_snapped",
    b1           = "b1_snapped",
    c            = "c_snapped",
    time_to_peak = "t2p_mean",
    FWHM         = "fwhm_mean"
  )
  if (!col_name %in% colnames(pop_avg)) {
    stop("Column '", col_name, "' not found in x$pop_avg.")
  }

  # Show the imputed values directly. Multiplying by mask_prop_NA would re-NA
  # the voxels that unmask_pop_avg just filled (defeats the imputation step).
  full_vector <- rep(NA, length(mask_prop_NA))
  full_vector[pop_avg$voxel] <- pop_avg[[col_name]]

  if (is.null(zlim)) {
    # Dynamic from the actual HRF grid (people can use different grids).
    zlim <- range(x$hrf_grid[[param]])
  }
  if (is.null(title)) {
    title <- paste0("Population avg ", param,
                    if (param %in% c("a1","b1","c")) " (snapped)" else " (imputed)")
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


#' Plot continuous (imputed, non-snapped) per-voxel HRF parameter
#'
#' Computes the per-voxel mean of \code{best_params_df$<param>} across
#' subjects (mask-filtered), sprays it onto a full cortical xifti, then
#' fills + smooths it with \code{unmask_xifti} (same imputation used inside
#' \code{regularize_allHRFs} for t2p/FWHM). Plots the resulting continuous
#' map. Unlike \code{type = "pop_avg"} (which snaps to the HRF grid), this
#' preserves the continuous between-grid values.
#'
#' @param x A \code{regularizeHRFs} object.
#' @param param One of \code{"a1"} or \code{"b1"}.
#' @param fname,title,shadows,material,NA_color,zlim Plot args (see
#'   \code{plot_pop_avg}).
#' @param surf_FWHM Surface smoothing FWHM (mm) for the imputation step.
#'   Default \code{4} matches \code{unmask_pop_avg}.
#' @param impute_method Passed to \code{unmask_xifti}. Default \code{"median"}.
#' @param ... Additional args passed to \code{plot.xifti}.
#'
#' @return Invisibly returns the plot result.
#' @keywords internal
#' @importFrom ciftiTools newdata_xifti
plot_pop_avg_continuous <- function(x,
                                    param = c("a1", "b1"),
                                    fname = NULL,
                                    title = NULL,
                                    shadows = 1,
                                    material = list(lit = TRUE, smooth = FALSE),
                                    NA_color = "#505560",
                                    zlim = NULL,
                                    surf_FWHM = 4,
                                    impute_method = "median",
                                    ...) {
  param <- match.arg(param)

  xii_template <- attr(x, "xii")
  if (is.null(xii_template)) {
    stop("attr(x, \"xii\") is missing. Re-run regularize_allHRFs to attach it.")
  }

  # Replicate Step 3 of regularize_allHRFs (R/regularize_allHRFs.R:98-103) exactly,
  # but for the requested param instead of t2p/FWHM. Same mask + same voxel filter.
  bpd <- x$best_params_df
  mask_prop_NA <- x$mask_prop_NA
  pop_mask_voxels <- which(!is.na(mask_prop_NA))
  agg <- stats::aggregate(
    stats::as.formula(paste(param, "~ voxel")),
    data = bpd[bpd$mask & bpd$voxel %in% pop_mask_voxels, ],
    FUN = mean
  )

  N <- length(mask_prop_NA)
  v <- rep(NA_real_, N); v[agg$voxel] <- agg[[param]]
  xii <- ciftiTools::newdata_xifti(xii_template, v)
  xii <- unmask_xifti(xii, method = impute_method, surf_FWHM = surf_FWHM,
                      impute_FUN = function(z) mean(z, na.rm = TRUE))

  if (is.null(zlim)) {
    # Dynamic zlim from the actual HRF grid (people can use different grids).
    zlim <- range(x$hrf_grid[[param]])
  }
  if (is.null(title)) {
    title <- paste0("Population avg ", param, " (imputed, continuous)")
  }

  plot_result <- plot(
    xii,
    zlim = zlim,
    color_mode = "diverging",
    fname = fname,
    title = title,
    shadows = shadows,
    material = material,
    NA_color = NA_color,
    ...
  )
  invisible(plot_result)
}


#' Plot HRF parameter frequency on the HRF grid
#'
#' Heatmap of how often each \code{(a1, b1, c)} grid combination appears,
#' faceted by \code{c}.
#'
#' @param x A \code{regularizeHRFs} object.
#' @param source Character. \code{"subjects"} (default) counts per-voxel
#'   per-subject best HRFs from \code{x$best_params_df} (one count per
#'   (voxel, subject) pair). \code{"pop_avg"} counts per-voxel regularized
#'   population HRFs from \code{x$pop_avg} (one count per voxel; the
#'   "averaged" view used in Fig 3).
#' @param mask Logical. Only relevant when \code{source = "subjects"}.
#'   If \code{TRUE}, restrict to voxel-subject pairs where the working-HRF
#'   activation mask is \code{TRUE}.
#' @param title Optional custom title.
#'
#' @return A \code{ggplot} object. Compose additional ggplot layers at the
#'   call site (e.g. \code{plot(x, type = "param_heatmap") + theme(...)}).
#'
#' @importFrom dplyr distinct summarise group_by select left_join filter mutate
#' @keywords internal
plot_param_heatmap <- function(x,
                               source = c("subjects", "pop_avg"),
                               mask = FALSE,
                               title = NULL) {
  if (!inherits(x, "regularizeHRFs")) {
    stop("Input object must be of class 'regularizeHRFs'.")
  }
  source <- match.arg(source)

  if (source == "subjects") {
    src_df <- x$best_params_df
    if (mask) src_df <- src_df %>% dplyr::filter(.data$mask == TRUE)
    if (is.null(title)) {
      title <- if (mask) "Best per-subject HRF params (masked)" else "Best per-subject HRF params"
    }
  } else {
    src_df <- data.frame(a1 = x$pop_avg$a1_snapped,
                         b1 = x$pop_avg$b1_snapped,
                         c  = x$pop_avg$c_snapped)
    if (is.null(title)) {
      title <- "Population-average HRF params (per voxel)"
    }
  }

  hrf_grid <- x$hrf_grid
  src_df <- src_df %>% dplyr::mutate(c = round(.data$c, 7))

  complete_grid <- hrf_grid %>%
    dplyr::select(.data$a1, .data$b1, .data$c) %>%
    dplyr::distinct() %>%
    dplyr::mutate(c = round(.data$c, 7))

  freq_df <- src_df %>%
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
    theme_few() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.spacing = unit(1.5, "lines"))
}


