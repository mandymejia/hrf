
#' Plot method for hrf_grid objects
#'
#' Creates diagnostic plots for HRF parameter grid results. Supports visualization
#' of raw HRFs, tapered HRFs, parameter grid heatmaps, individual HRFs, and multiple HRF overlays.
#'
#' @param x An object of class \code{"hrf_grid"} from \code{\link{generate_hrf_grid}}.
#' @param type Character. Type of plot to generate. Options are:
#'   \itemize{
#'     \item \code{"hrfs"} – plot all HRFs from the parameter grid (raw or tapered based on \code{tapered} arg).
#'     \item \code{"param_grid"} – plot parameter grid heatmaps of time-to-peak and FWHM.
#'     \item \code{"single_hrf"} – plot a single HRF specified by \code{hrf_idx}.
#'     \item \code{"multiple_hrf"} – plot multiple overlapping HRFs specified by \code{hrf_idx} vector.
#'     \item \code{"slices"} – two-panel didactic view: vary a1 (with b1 fixed at canonical),
#'       and vary b1 (with a1 fixed at canonical). Canonical HRF drawn thick + black.
#'   }
#' @param hrf_idx Integer or integer vector. Row index(es) for plotting HRF(s).
#'   Single value for \code{type = "single_hrf"}, vector for \code{type = "multiple_hrf"}.
#' @param tapered Logical. Whether to plot tapered version. Applies to \code{type = "hrfs"}, 
#'   \code{type = "single_hrf"}, and \code{type = "multiple_hrf"}. Default is TRUE.
#' @param colors Character vector. Colors for each HRF (required for \code{type = "multiple_hrf"}).
#'   Must have same length as \code{hrf_idx}.
#' @param ... Additional arguments passed to the specific plotting function.
#'
#' @return Invisibly returns the result of the specific plotting function.
#'   Plots are rendered to the current graphics device.
#'
#' @seealso \code{\link{generate_hrf_grid}}, \code{\link{compute_hrf_metrics}}
#'
#' @examples
#' \dontrun{
#' # Plot all raw HRFs
#' plot(my_hrf_grid, type = "hrfs", tapered = FALSE)
#'
#' # Plot all tapered HRFs (default)
#' plot(my_hrf_grid, type = "hrfs")
#' plot(my_hrf_grid, type = "hrfs", tapered = TRUE)
#'
#' # Plot a single HRF
#' plot(my_hrf_grid, type = "single_hrf", hrf_idx = 5)
#' plot(my_hrf_grid, type = "single_hrf", hrf_idx = 5, tapered = TRUE)
#'
#' # Plot multiple overlapping HRFs
#' plot(my_hrf_grid, type = "multiple_hrf", hrf_idx = c(1, 5, 10),
#'      colors = c("#2c7fb8", "#d95f02", "#7570b3"))
#'
#' # Plot parameter grid metrics
#' plot(my_hrf_grid, type = "param_grid")
#' }
#'
#' @export
#' @param show_t2p_limits Logical. For \code{type = "hrfs"}, whether to draw the
#'   dashed vertical lines at 2 s and 10 s marking the physiologically plausible
#'   time-to-peak range. Default is \code{TRUE}.
#'
plot.hrf_grid <- function(x, type = c("hrfs", "param_grid", "single_hrf", "multiple_hrf", "slices"),
                          hrf_idx = 1, tapered = TRUE, colors = NULL,
                          show_t2p_limits = TRUE, ...) {
  type <- match.arg(type)
  switch(type,
         hrfs         = plot_hrfs(x, tapered = tapered, show_t2p_limits = show_t2p_limits, ...),
         param_grid   = plot_param_grid_metrics(x, ...),
         single_hrf   = plot_hrf_single(x, hrf_idx = hrf_idx, tapered = tapered, ...),
         multiple_hrf = plot_hrf_multiple(x, hrf_idx = hrf_idx, colors = colors, tapered = tapered, ...),
         slices       = plot_hrfs_slices(x, ...)
  )
}

#' Plot all HRFs from parameter grid - wrapper (internal)
#'
#' Internal wrapper function that calls either the raw or tapered HRF plotting
#' function based on the \code{tapered} argument.
#'
#' @param hrf_grid Data frame with HRF parameters and call_params attributes
#' @param tapered Logical. If TRUE, plots tapered HRFs; if FALSE, plots raw HRFs.
#'   Default is TRUE.
#' @param show_t2p_limits Logical. Whether to draw dashed vertical lines at 2 s
#'   and 10 s marking the physiologically plausible time-to-peak range. Default
#'   is TRUE.
#'
#' @return A ggplot object
#'
#' @keywords internal
plot_hrfs <- function(hrf_grid, tapered = TRUE, show_t2p_limits = TRUE, ...) {
  if (tapered) {
    plot_hrfs_all_tapered(hrf_grid, show_t2p_limits = show_t2p_limits, ...)
  } else {
    plot_hrfs_all(hrf_grid, show_t2p_limits = show_t2p_limits, ...)
  }
}


#' Plot multiple HRFs from parameter grid (internal)
#'
#' Internal helper for \code{plot.hrf_grid}. Plots multiple overlapping HRFs
#' specified by index vector with custom colors.
#'
#' @param hrf_grid Data frame with columns a1, b1, a2, b2, c and call_params attributes
#' @param hrf_idx Integer vector. Row indices of the HRFs to plot
#' @param colors Character vector. Colors for each HRF (must match length of hrf_idx)
#' @param tapered Logical. Whether to plot tapered versions
#'
#' @return A ggplot object
#'
#' @keywords internal
plot_hrf_multiple <- function(hrf_grid, hrf_idx, colors, tapered = TRUE) {
  stopifnot(all(hrf_idx >= 1), all(hrf_idx <= nrow(hrf_grid)))
  if (length(hrf_idx) != length(colors)) stop("Length of 'hrf_idx' must equal length of 'colors'")
  
  # Grab TR from attributes
  TR <- 2
  
  # Make time vector
  inds <- seq(1/100, 30, 1/100) * TR
  
  # Initialize empty data frame
  plot_df <- data.frame()
  
  # Loop through each HRF index
  for (i in seq_along(hrf_idx)) {
    idx <- hrf_idx[i]
    params <- hrf_grid[idx, ]
    
    # Call HRF_calc with those params
    hrf_vals <- hrf::HRF_calc(
      t = inds, 
      deriv = 0,
      a1 = params$a1,
      b1 = params$b1,
      a2 = params$a2,
      b2 = params$b2,
      c = params$c
    )
    
    # If tapered, recalculate with taper params
    if (tapered && params$c > 0) {
      peak2_time <- inds[which.min(hrf_vals)]
      taper_start <- min(peak2_time, 25)
      
      if (abs(hrf_vals[which.min(abs(inds - 30))]) > .HRF_THRESHOLD) {
        hrf_vals <- hrf::HRF_calc(
          t = inds, 
          deriv = 0,
          a1 = params$a1,
          b1 = params$b1,
          a2 = params$a2,
          b2 = params$b2,
          c = params$c,
          taper_start = taper_start,
          taper_end = 30,
          taper_power = 1
        )
      }
    }
    
    # Add to data frame with HRF label
    temp_df <- data.frame(
      sec = inds, 
      HRF = hrf_vals,
      hrf_label = sprintf("HRF #%d (a1=%.1f, b1=%.2f, c=%.3f)", 
                          idx, params$a1, params$b1, params$c)
    )
    plot_df <- rbind(plot_df, temp_df)
  }
  # end_time <- calculate_hrf_endpoint(hrf_grid, hrf_idx, TR, tapered)
  end_time <- 30
  # Multiple HRF overlay plot
  p <- ggplot(plot_df, aes(x = .data$sec, y = .data$HRF, 
                           color = .data$hrf_label, group = .data$hrf_label)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +
    scale_color_manual(values = colors) +
    labs(title = sprintf("HRF Comparison%s", 
                         ifelse(tapered, "", " (non-tapered)")),
         x = "Time (seconds)", 
         y = "HRF Response",
         color = "HRF") +
    theme_minimal() +
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c("right", "top"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
    coord_cartesian(xlim = c(0, end_time))
  return(p)
}


#' Calculate appropriate x-axis endpoint for HRF plots
#'
#' @param hrf_grid Data frame with HRF parameters
#' @param hrf_idx Integer vector of row indices to check
#' @param TR Time repetition value
#' @param tapered Logical, whether HRFs are tapered
#' @param threshold Numeric, threshold for considering HRF "zero"
#' @param buffer Numeric, seconds to add after last non-zero value
#' @param max_time Numeric, maximum allowed time
#'
#' @return Numeric endpoint for x-axis
#' @keywords internal
calculate_hrf_endpoint <- function(hrf_grid, hrf_idx, TR, tapered = TRUE, 
                                   threshold = .HRF_THRESHOLD, buffer = 2, max_time = 50) {
  if (tapered) return(30)
  
  inds <- seq(1/100, 30, 1/100) * TR
  max_end <- 0
  
  for (idx in hrf_idx) {
    params <- hrf_grid[idx, ]
    hrf_vals <- hrf::HRF_calc(
      t = inds, 
      deriv = 0,
      a1 = params$a1,
      b1 = params$b1,
      a2 = params$a2,
      b2 = params$b2,
      c = params$c
    )
    last_nonzero <- max(which(abs(hrf_vals) > threshold))
    max_end <- max(max_end, inds[last_nonzero])
  }
  
  end_time <- min(max_end + buffer, max_time)
  message(sprintf("Longest HRF ends at %.1f seconds, plot ends at %.1f seconds", 
                  max_end, end_time))
  return(end_time)
}

#' Plot a single HRF from parameter grid (internal)
#'
#' Internal helper for \code{plot.hrf_grid}. Plots a single HRF specified by index.
#'
#' @param hrf_grid Data frame with columns a1, b1, a2, b2, c and call_params attributes
#' @param hrf_idx Integer. Row index of the HRF to plot
#' @param tapered Logical. Whether to plot tapered version
#'
#' @return A ggplot object
#'
#' @keywords internal
plot_hrf_single <- function(hrf_grid, hrf_idx = 1, tapered = TRUE) {
  stopifnot(hrf_idx <= nrow(hrf_grid), hrf_idx >= 1)
  # Grab the row
  params <- hrf_grid[hrf_idx, ]
  
  # Grab TR from attributes
  TR <- 2
  
  # Make time vector
  inds <- seq(1/100, 30, 1/100) * TR
  
  # Call HRF_calc with those params
  hrf_vals <- hrf::HRF_calc(
    t = inds, 
    deriv = 0,
    a1 = params$a1,
    b1 = params$b1,
    a2 = params$a2,
    b2 = params$b2,
    c = params$c
  )
  
  # If tapered, recalculate with taper params
  if (tapered && params$c > 0) {
    peak2_time <- inds[which.min(hrf_vals)]
    taper_start <- min(peak2_time, 25)
    
    if (abs(hrf_vals[which.min(abs(inds - 30))]) > .HRF_THRESHOLD) {
      hrf_vals <- hrf::HRF_calc(
        t = inds, 
        deriv = 0,
        a1 = params$a1,
        b1 = params$b1,
        a2 = params$a2,
        b2 = params$b2,
        c = params$c,
        taper_start = taper_start,
        taper_end = 30,
        taper_power = 1
      )
    }
  }

  end_time <- 30  # default for tapered
  # if (!tapered) {
  #   threshold <- .HRF_THRESHOLD
  #   last_nonzero <- max(which(abs(hrf_vals) > threshold))
  #   end_time <- min(inds[last_nonzero] + 2, 60)  # Add 2 sec buffer, max 60
  # }
  
  # Plot it with ggplot
  plot_df <- data.frame(sec = inds, HRF = hrf_vals)
  
  p <- ggplot(plot_df, aes(x = .data$sec, y = .data$HRF)) +
    geom_line(color = '#2c7fb8', linewidth = 1) +
    geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +
    labs(title = sprintf("HRF #%d%s: a1=%.1f, b1=%.2f, c=%.3f", 
                         hrf_idx, 
                         ifelse(tapered && params$c > 0, " (tapered)", ""),
                         params$a1, params$b1, params$c),
         x = "Time (seconds)", 
         y = "HRF Response") +
    theme_minimal() +
    coord_cartesian(xlim = c(0, end_time))

    return(p)
}


#' Compute HRF metrics from parameter grid
#'
#' Takes an HRF parameter grid and computes all metrics including time-to-peak,
#' FWHM, peak times, and the actual HRF values at super-resolution.
#' All parameters are extracted from the grid's attributes.
#'
#' @param hrf_grid Data frame with columns: a1, b1, c, a2, b2 and call_params attributes
#'
#' @return List with two elements:
#'   \item{hrf_params}{Data frame with all HRF parameters and computed metrics}
#'   \item{hrf_df}{Data frame with super-resolution HRF values for plotting}
#'
#' @keywords internal
compute_hrf_metrics <- function(hrf_grid) {
  # Extract parameters from attributes
  call_params <- attr(hrf_grid, "call_params")

  # Get parameters from attributes
  TR <- 2
  sr_factor <- call_params$sr_factor
  max_time <- 20  # Not in attributes, so using default

  # Initialize hrf_params with the grid data
  hrf_params <- hrf_grid
  hrf_params$TR <- TR

  # Calculate analytical metrics
  hrf_params$shape1 <- hrf_params$a1 / hrf_params$b1
  hrf_params$rate1 <- hrf_params$TR / hrf_params$b1
  hrf_params$TRs_to_peak <- (hrf_params$shape1 - 1) / hrf_params$rate1
  hrf_params$time_to_peak <- hrf_params$TRs_to_peak * hrf_params$TR

  # Initialize storage for computed metrics
  hrf_params$FWHM <- NA
  hrf_params$time_to_end <- NA
  hrf_params$peak2_time <- NA

  # Create super-resolution time vector
  inds <- seq(1/sr_factor, max_time, 1/sr_factor) * TR
  nt <- length(inds)

  # Pre-allocate hrf_df components
  n_params <- nrow(hrf_params)
  total_rows <- n_params * nt

  # Initialize hrf_df with repeated parameters
  hrf_df <- data.frame(
    a1 = rep(hrf_params$a1, each = nt),
    b1 = rep(hrf_params$b1, each = nt),
    TR = TR,
    sec = rep(inds, n_params),
    a2 = rep(hrf_params$a2, each = nt),
    b2 = rep(hrf_params$b2, each = nt),
    c = rep(hrf_params$c, each = nt)
  )

  # Initialize HRF storage
  hrf_df$HRF <- NA
  hrf_df$HRF_tapered <- NA
  hrf_df$FWHM <- NA
  hrf_df$time_to_end <- NA
  hrf_df$peak2_time <- NA

  for (ii in 1:n_params) {

    a1_ii <- hrf_params$a1[ii]
    b1_ii <- hrf_params$b1[ii]
    a2_ii <- hrf_params$a2[ii]
    b2_ii <- hrf_params$b2[ii]
    c_ii <- hrf_params$c[ii]

    # Calculate HRF at super-resolution
    hrf_ii <- hrf::HRF_calc(t = inds, deriv = 0, a1 = a1_ii,b1 = b1_ii, a2 = a2_ii, b2 = b2_ii, c = c_ii)

    # Find peak value and location
    peak_val <- max(hrf_ii)
    peak_idx <- which.max(hrf_ii)

    # Calculate FWHM
    vals_left <- hrf_ii[1:peak_idx]
    vals_right <- hrf_ii[peak_idx:length(hrf_ii)]

    # Find half-maximum crossings
    x1 <- min(which(vals_left > 0.5 * peak_val))
    x2 <- max(which(vals_right > 0.5 * peak_val))

    time1 <- inds[x1]
    time2 <- inds[peak_idx + x2 - 1]
    FWHM_ii <- time2 - time1

    # Find undershoot peak time
    if (c_ii > 0) {
      peak2_idx <- which.min(hrf_ii)
    } else {
      peak2_idx <- peak_idx  # If no undershoot, use main peak
    }
    peak2_time_ii <- inds[peak2_idx]

    # Find time when HRF resolves (drops below threshold)
    time_to_end_ii <- inds[max(which(abs(hrf_ii) > .HRF_THRESHOLD))]

    # Store in hrf_params
    hrf_params$FWHM[ii] <- FWHM_ii
    hrf_params$time_to_end[ii] <- time_to_end_ii
    hrf_params$peak2_time[ii] <- peak2_time_ii

    # Calculate tapered version if needed
    taper_start <- min(peak2_time_ii, 25)
    if (c_ii > 0 && abs(hrf_ii[which.min(abs(inds - 30))]) > .HRF_THRESHOLD) {
      hrf_tapered_ii <- hrf::HRF_calc(
        t = inds, deriv = 0,  a1 = a1_ii, b1 = b1_ii,a2 = a2_ii, b2 = b2_ii,
        c = c_ii,taper_start = taper_start,taper_end = 30,taper_power = 1)
    } else {
      hrf_tapered_ii <- hrf_ii
    }

    # Store in hrf_df
    rows_ii <- ((ii - 1) * nt + 1):(ii * nt)
    hrf_df$HRF[rows_ii] <- hrf_ii
    hrf_df$HRF_tapered[rows_ii] <- hrf_tapered_ii
    hrf_df$FWHM[rows_ii] <- FWHM_ii
    hrf_df$time_to_end[rows_ii] <- time_to_end_ii
    hrf_df$peak2_time[rows_ii] <- peak2_time_ii
  }

  # Add factor versions for plotting
  hrf_df$c_num <- hrf_df$c

  # Create labels for c values - handle 0 and 1/6 specially, others generically
  c_unique <- sort(unique(hrf_df$c))
  c_labels <- sapply(c_unique, function(c_val) {
    if (c_val == 0) {
      'No Undershoot (c=0)'
    } else if (abs(c_val - 1/6) < 1e-10) {
      'With Undershoot (c=1/6)'
    } else {
      paste0('c = ', round(c_val, 3))
    }
  })

  hrf_df$c <- factor(
    hrf_df$c,
    levels = c_unique,
    labels = c_labels
  )

  # Format b1 for plotting
  hrf_df$b1 <- paste0('b1 = ', hrf_df$b1)

  return(list(
    hrf_params = hrf_params,
    hrf_df = hrf_df
  ))
}

#' Plot raw HRFs from parameter grid (internal)
#'
#' Internal helper for \code{plot.hrf_grid}. Generates a faceted plot of all raw
#' HRFs from the parameter grid.
#'
#' @param hrf_grid Data frame with HRF parameters (created by
#'   \code{generate_hrf_grid}).
#'
#' @importFrom ggthemes theme_few
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
plot_hrfs_all <- function(hrf_grid, show_t2p_limits = TRUE,
                           canonical_a1 = 6, canonical_b1 = 1,
                           dashed_a1 = c(3, 5, 7, 9, 11)) {

  hrf_results <- compute_hrf_metrics(hrf_grid)
  hrf_df <- hrf_results$hrf_df

  n_c <- length(unique(hrf_df$c))

  a1_vals  <- sort(unique(hrf_df$a1))
  roy      <- ciftiTools::ROY_BIG_BL(min(a1_vals), max(a1_vals), mid = canonical_a1)
  roy      <- roy[order(roy$value), ]
  swatches <- grDevices::colorRampPalette(roy$col)(length(a1_vals))
  ltypes   <- ifelse(a1_vals %in% dashed_a1, "dashed", "solid")
  names(swatches) <- as.character(a1_vals)
  names(ltypes)   <- as.character(a1_vals)

  canonical_grid <- hrf_grid[hrf_grid$a1 == canonical_a1 &
                             hrf_grid$b1 == canonical_b1, ]
  class(canonical_grid) <- c("hrf_grid", "data.frame")
  canonical_df <- if (nrow(canonical_grid) > 0) {
    compute_hrf_metrics(canonical_grid)$hrf_df
  } else NULL

  p <- ggplot(hrf_df, aes(x = .data$sec, y = .data$HRF,
                          color    = factor(.data$a1),
                          linetype = factor(.data$a1),
                          group    = factor(.data$a1))) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_line(linewidth = 0.8) +
    xlim(0, 30) +
    scale_color_manual(values = swatches, name = "a1") +
    scale_linetype_manual(values = ltypes, guide = "none") +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(y = "HRF", x = "Time (s)") +
    theme_few() +
    theme(legend.position = 'right',
          panel.border = element_blank(),
          axis.line    = element_line(color = "black"))

  if (!is.null(canonical_df)) {
    p <- p + geom_line(data = canonical_df, color = "black",
                       linewidth = 1.4, linetype = "solid",
                       inherit.aes = FALSE,
                       mapping = aes(x = .data$sec, y = .data$HRF))
  }

  if (show_t2p_limits) {
    p <- p + geom_vline(xintercept = c(2, 10), linetype = 2, color = 'gray')
  }

  p <- p + if (n_c == 1) {
    facet_wrap(~ b1, ncol = 1, strip.position = "right")
  } else {
    facet_grid(b1 ~ c)
  }

  return(p)
}

#' Plot HRF parameter grid with computed metrics (internal)
#'
#' Internal helper for \code{plot.hrf_grid}. Generates a two-panel plot showing
#' the HRF parameter space with time-to-peak (top) and FWHM (bottom) values
#' displayed as heatmaps with text labels, faceted by undershoot parameter.
#'
#' @param hrf_grid A data frame of HRF parameters with \code{call_params}
#'   attributes, typically created by \code{generate_hrf_grid()}.
#'
#' @importFrom ggthemes theme_few
#'
#' @return A grob object from \code{gridExtra::grid.arrange} combining the two plots.
#'
#' @keywords internal
plot_param_grid_metrics <- function(hrf_grid) {
  # Compute HRF metrics
  hrf_results <- compute_hrf_metrics(hrf_grid)
  hrf_params <- hrf_results$hrf_params

  # Create c_label for faceting (used only when multiple c values present)
  hrf_params$c_label <- ifelse(hrf_params$c == 0,
                               'No Undershoot (c=0)',
                               paste0('With Undershoot (c=', round(hrf_params$c, 3), ')'))

  n_c <- length(unique(hrf_params$c))

  # Get unique a1 and b1 values for grid lines
  a1_vals <- sort(unique(hrf_params$a1))
  b1_vals <- sort(unique(hrf_params$b1))

  # Calculate grid line positions (between values)
  a1_grid <- c(min(a1_vals) - 0.5, a1_vals + 0.5)
  b1_grid <- c(min(b1_vals) - 0.125, b1_vals + 0.125)

  # SPM canonical marker: white hollow square at (a1=6, b1=1) when present.
  # geom_rect (not geom_tile) so the fill scale doesn't paint an NA ghost.
  step_a1 <- if (length(a1_vals) > 1) min(diff(a1_vals)) else 1
  step_b1 <- if (length(b1_vals) > 1) min(diff(b1_vals)) else 0.25
  show_canonical <- 6 %in% a1_vals && 1 %in% b1_vals
  canonical_layer <- if (show_canonical) {
    geom_rect(data = data.frame(a1 = 6, b1 = 1),
              aes(xmin = .data$a1 - step_a1 / 2, xmax = .data$a1 + step_a1 / 2,
                  ymin = .data$b1 - step_b1 / 2, ymax = .data$b1 + step_b1 / 2),
              fill = NA, color = "white", linewidth = 1.2, inherit.aes = FALSE)
  } else NULL

  base_theme <- theme_few() +
    theme(panel.spacing = unit(1.5, "lines"),
          legend.position = "right",
          plot.margin = margin(5, 0, 5, 0))

  # Conditional facet: drop facet when only one c value present.
  facet_layer <- if (n_c == 1) NULL else facet_grid(. ~ c_label)

  # Time-to-peak heatmap
  p1 <- ggplot(hrf_params, aes(x = .data$a1, y = .data$b1, fill = .data$time_to_peak)) +
    geom_hline(yintercept = b1_grid, alpha = 0.1) +
    geom_vline(xintercept = a1_grid, alpha = 0.1) +
    geom_tile(alpha = 0.8) +
    geom_text(aes(label = round(.data$time_to_peak, 1))) +
    canonical_layer +
    scale_fill_viridis_c('Time to Peak  ', option = 'C') +
    scale_x_continuous(breaks = a1_vals, expand = c(0, 0)) +
    scale_y_continuous(breaks = b1_vals, expand = c(0, 0)) +
    facet_layer +
    base_theme

  # FWHM heatmap
  p2 <- ggplot(hrf_params, aes(x = .data$a1, y = .data$b1, fill = .data$FWHM)) +
    geom_hline(yintercept = b1_grid, alpha = 0.1) +
    geom_vline(xintercept = a1_grid, alpha = 0.1) +
    geom_tile(alpha = 0.8) +
    geom_text(aes(label = round(.data$FWHM, 1))) +
    canonical_layer +
    scale_fill_viridis_c('Width (FWHM)', option = 'D') +
    scale_x_continuous(breaks = a1_vals, expand = c(0, 0)) +
    scale_y_continuous(breaks = b1_vals, expand = c(0, 0)) +
    facet_layer +
    base_theme

  # Side-by-side layout
  combined_plot <- gridExtra::grid.arrange(p1, p2, ncol = 2)

  return(combined_plot)
}

#' Plot tapered HRFs from parameter grid (internal)
#'
#' Two-panel "slices" view of the HRF parameter grid (internal)
#'
#' Vary a1 with b1 fixed at canonical, then vary b1 with a1 fixed at canonical.
#' Canonical HRF drawn thick + black in both panels. Discrete diverging
#' (ROY_BIG_BL) palette centered at the canonical value.
#'
#' @param hrf_grid Data frame with HRF parameters from \code{generate_hrf_grid}.
#' @param canonical_a1,canonical_b1 The fixed canonical parameter values.
#'   Defaults: SPM canonical (\code{a1 = 6}, \code{b1 = 1}).
#' @param c_val Single c value to use. If \code{hrf_grid} has multiple c values,
#'   the grid is filtered to this one. Default \code{1/6}.
#' @param dashed_a1,dashed_b1 Values to draw as dashed (alternating to help
#'   distinguish neighboring swatches in the discrete palette).
#'
#' @keywords internal
#' @return A \code{grid.arrange} grob with two stacked panels.
plot_hrfs_slices <- function(hrf_grid, canonical_a1 = 6, canonical_b1 = 1,
                              c_val = 1/6,
                              dashed_a1 = c(3, 5, 7, 9, 11),
                              dashed_b1 = c(0.75, 1.25, 1.75)) {

  c_vals <- unique(hrf_grid$c)
  if (length(c_vals) > 1) {
    if (!(c_val %in% c_vals)) {
      stop(sprintf("c_val=%s not in grid (c values: %s). Pass a valid c_val.",
                   c_val, paste(c_vals, collapse = ", ")))
    }
    hrf_grid <- hrf_grid[hrf_grid$c == c_val, ]
    class(hrf_grid) <- c("hrf_grid", "data.frame")
  }

  hrf_df <- compute_hrf_metrics(hrf_grid)$hrf_df
  # plot_hrfs_all_tapered re-prefixes b1 strip values with "b1 = "; we use the
  # numeric value here.
  if (is.character(hrf_df$b1)) {
    hrf_df$b1 <- as.numeric(sub("b1 = ", "", hrf_df$b1))
  }

  make_panel <- function(df, color_var, fixed_label, mid_val, dashed_vals) {
    vals <- sort(unique(df[[color_var]]))
    roy  <- ciftiTools::ROY_BIG_BL(min(vals), max(vals), mid = mid_val)
    roy  <- roy[order(roy$value), ]
    swatches <- grDevices::colorRampPalette(roy$col)(length(vals))

    ltypes <- ifelse(vals %in% dashed_vals, "dashed", "solid")
    names(ltypes) <- as.character(vals)

    canonical <- df[df$a1 == canonical_a1 & df$b1 == canonical_b1, ]

    ggplot(df, aes(x = .data$sec, y = .data$HRF_tapered,
                   color    = factor(.data[[color_var]]),
                   linetype = factor(.data[[color_var]]),
                   group    = factor(.data[[color_var]]))) +
      geom_hline(yintercept = 0, color = 'gray') +
      geom_line(linewidth = 0.8) +
      geom_line(data = canonical, color = "black", linewidth = 1.4,
                linetype = "solid", inherit.aes = FALSE,
                mapping = aes(x = .data$sec, y = .data$HRF_tapered)) +
      scale_color_manual(values = swatches, name = color_var) +
      scale_linetype_manual(values = ltypes, guide = "none") +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      xlim(0, 30) +
      labs(x = "Time (s)", y = "HRF",
           title = sprintf("Effect of varying %s (%s)", color_var, fixed_label)) +
      theme_few(base_size = 11) +
      theme(legend.position = "right",
            panel.border = element_blank(),
            axis.line    = element_line(color = "black"))
  }

  pA <- make_panel(hrf_df[hrf_df$b1 == canonical_b1, ], "a1",
                   sprintf("b1 = %g", canonical_b1), canonical_a1, dashed_a1)
  pB <- make_panel(hrf_df[hrf_df$a1 == canonical_a1, ], "b1",
                   sprintf("a1 = %g", canonical_a1), canonical_b1, dashed_b1)

  gridExtra::grid.arrange(pA, pB, nrow = 2)
}

#' Plot all tapered HRFs from parameter grid (internal)
#'
#' Internal helper for \code{plot.hrf_grid}. Generates a faceted plot of all
#' tapered hemodynamic response functions (HRFs) from the parameter grid.
#' Tapered HRFs are used when the response does not resolve to baseline
#' by 30 seconds.
#'
#' @param hrf_grid A data frame of HRF parameters with \code{call_params}
#'   attributes, typically created by \code{generate_hrf_grid()}.
#' @param show_t2p_limits Logical. Whether to draw dashed vertical lines at 2 s
#'   and 10 s marking the physiologically plausible time-to-peak range.
#'
#' @importFrom ggthemes theme_few
#'
#' @return A \code{ggplot} object.
#'
#' @keywords internal
plot_hrfs_all_tapered <- function(hrf_grid, show_t2p_limits = TRUE,
                                   canonical_a1 = 6, canonical_b1 = 1,
                                   dashed_a1 = c(3, 5, 7, 9, 11)) {

  hrf_results <- compute_hrf_metrics(hrf_grid)
  hrf_df <- hrf_results$hrf_df

  n_c <- length(unique(hrf_df$c))

  # Discrete a1 scale: one ROY_BIG_BL swatch per integer a1, centered at
  # canonical. Odd-valued a1 drawn dashed so neighboring swatches are easier
  # to distinguish.
  a1_vals  <- sort(unique(hrf_df$a1))
  roy      <- ciftiTools::ROY_BIG_BL(min(a1_vals), max(a1_vals), mid = canonical_a1)
  roy      <- roy[order(roy$value), ]
  swatches <- grDevices::colorRampPalette(roy$col)(length(a1_vals))
  ltypes   <- ifelse(a1_vals %in% dashed_a1, "dashed", "solid")
  names(swatches) <- as.character(a1_vals)
  names(ltypes)   <- as.character(a1_vals)

  # Canonical HRF (a1==canonical_a1, b1==canonical_b1) for thick-black overlay.
  canonical_grid <- hrf_grid[hrf_grid$a1 == canonical_a1 &
                             hrf_grid$b1 == canonical_b1, ]
  class(canonical_grid) <- c("hrf_grid", "data.frame")
  canonical_df <- if (nrow(canonical_grid) > 0) {
    compute_hrf_metrics(canonical_grid)$hrf_df
  } else NULL

  p <- ggplot(hrf_df, aes(x = .data$sec, y = .data$HRF_tapered,
                          color    = factor(.data$a1),
                          linetype = factor(.data$a1),
                          group    = factor(.data$a1))) +
    geom_hline(yintercept = 0, color = 'gray') +
    geom_line(linewidth = 0.8) +
    xlim(0, 30) +
    scale_color_manual(values = swatches, name = "a1") +
    scale_linetype_manual(values = ltypes, guide = "none") +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(y = "HRF", x = "Time (s)") +
    theme_few() +
    theme(legend.position = 'right',
          panel.border = element_blank(),
          axis.line    = element_line(color = "black"))

  if (!is.null(canonical_df)) {
    p <- p + geom_line(data = canonical_df, color = "black",
                       linewidth = 1.4, linetype = "solid",
                       inherit.aes = FALSE,
                       mapping = aes(x = .data$sec, y = .data$HRF_tapered))
  }

  if (show_t2p_limits) {
    p <- p + geom_vline(xintercept = c(2, 10), linetype = 2, color = 'gray')
  }

  # Single c value → drop the c facet; multiple c → keep both.
  p <- p + if (n_c == 1) {
    facet_wrap(~ b1, ncol = 1, strip.position = "right")
  } else {
    facet_grid(b1 ~ c)
  }

  return(p)
}
