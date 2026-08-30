#' Plot method for regularizeHRFs objects
#'
#' Diagnostic plots for the offset-based regularize output.
#' \itemize{
#'   \item \code{"pop_best"} – population-average HRF parameter map (brain
#'     surface, grid-snapped values from \code{x$pop_best}).
#'   \item \code{"param_heatmap"} – frequency heatmap of per-subject best HRF
#'     parameters on the underlying HRF grid.
#'   \item \code{"param_grid"}, \code{"slices"}, \code{"hrfs"} – grid-shape
#'     methodology figures drawn from the grid this run actually fit
#'     (\code{x$hrf_grid}), delegated to \code{\link{plot.hrf_grid}}.
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
plot.regularizeHRFs <- function(x, type = c("pop_best", "param_heatmap",
                                            "param_grid", "slices", "hrfs"), ...) {
  type <- match.arg(type)

  switch(type,
    pop_best      = plot_pop_best(x, ...),
    param_heatmap = plot_param_heatmap(x, ...),
    # Grid-shape figures: plot the grid this run ACTUALLY fit (x$hrf_grid), so
    # param_grid / slices / hrfs always match the data (e.g. full-grid = 66,
    # standard = 48) instead of a hardcoded regenerated grid.
    param_grid = ,
    slices     = ,
    hrfs       = {
      g <- x$hrf_grid
      class(g) <- c("hrf_grid", "data.frame")
      plot(g, type = type, ...)
    }
  )
}


#' Plot Population-Average HRF Parameter Map (Brain Surface)
#'
#' Renders the population-averaged HRF parameter from \code{x$pop_best}
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
#' @param mask Logical. \code{FALSE} (default) shows the modal-filled (unmasked)
#'   pop_best; \code{TRUE} re-applies \code{mask_prop_NA} so non-activated voxels
#'   stay NA (the masked view).
#' @param ... Additional args passed to \code{plot.xifti}.
#'
#' @return Invisibly returns the plot result.
#' @keywords internal
#' @importFrom ciftiTools newdata_xifti
plot_pop_best <- function(x,
                         param = c("a1", "b1", "c", "time_to_peak", "FWHM"),
                         fname = NULL,
                         title = NULL,
                         shadows = 1,
                         material = list(lit = TRUE, smooth = FALSE),
                         NA_color = "#505560",
                         zlim = NULL,
                         mask = FALSE,
                         ...) {
  param <- match.arg(param)

  xii_template <- attr(x, "xii")
  if (is.null(xii_template)) {
    stop("attr(x, \"xii\") is missing. This regularizeHRFs result was likely built ",
         "before the xii template was attached. Re-run regularize_allHRFs().")
  }

  pop_best <- x$pop_best
  mask_prop_NA <- x$mask_prop_NA
  # All columns are the winning grid point's values (pop_rss argmin); t2p_mean /
  # fwhm_mean are the winner's t2p / FWHM (legacy "_mean" names), not subject means.
  col_name <- switch(param,
    a1           = "a1",
    b1           = "b1",
    c            = "c",
    time_to_peak = "t2p_mean",
    FWHM         = "fwhm_mean"
  )
  if (!col_name %in% colnames(pop_best)) {
    stop("Column '", col_name, "' not found in x$pop_best.")
  }

  # pop_best is already modal-filled (unmasked). Default shows it as-is; mask=TRUE
  # re-NAs the off-pop-mask voxels -- a pure mask, never a weighting (non-weighted RSS).
  full_vector <- rep(NA, length(mask_prop_NA))
  full_vector[pop_best$voxel] <- pop_best[[col_name]]
  if (mask) full_vector[is.na(mask_prop_NA)] <- NA

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


#' Plot HRF parameter frequency on the HRF grid
#'
#' Heatmap of how often each \code{(a1, b1, c)} grid combination appears,
#' faceted by \code{c}.
#'
#' @param x A \code{regularizeHRFs} object.
#' @param source Character. \code{"subjects"} (default) counts per-voxel
#'   per-subject best HRFs from \code{x$best_params_df} (one count per
#'   (voxel, subject) pair). \code{"pop_best"} counts per-voxel regularized
#'   population HRFs from \code{x$pop_best} (one count per voxel; the
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
                               source = c("subjects", "pop_best"),
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
    # Count only ACTIVATED voxels -- the real per-voxel RSS picks. Non-activated
    # voxels are modal-filled by modal_unmask (imputation, not real picks), so
    # counting them dumps a large artificial spike on the single modal grid point.
    keep   <- !is.na(x$mask_prop_NA[x$pop_best$voxel])
    src_df <- data.frame(a1 = x$pop_best$a1[keep],
                         b1 = x$pop_best$b1[keep],
                         c  = x$pop_best$c[keep])
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


