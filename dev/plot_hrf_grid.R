library(here)
library(qs)
library(ggplot2)
library(ggthemes)
library(gridExtra)


ciftiTools::ciftiTools.setOption('wb_path','/Applications/workbench/bin_macosxub/wb_command')

devtools::load_all("~/Documents/Github/hrf-z") # hrf-z is Zeshawn's branch with latest hrf package changes.
# devtools::load_all("~/Documents/Github/hrf-HRFcalc-mods") # Temporary while cleaning up hrf-z branch

# ------------------------------------------------------------------------------
#                              HRF GRID GENERATION
# ------------------------------------------------------------------------------
hrf_grid <- default_hrf_grid
compute_hrf_metrics <- function(hrf_grid) {
  call_params <- attr(hrf_grid, "call_params")

  TR <- call_params$TR

  a1 <- seq(call_params$a1_range[[1]], call_params$a1_range[[2]], by=call_params$a1_step)
  b1 <- seq(call_params$b1_range[[1]], call_params$b1_range[[2]], by=call_params$b1_step)
  hrf_params <- expand.grid(a1 = a1,
                            b1 = b1,
                            TR=TR)

  #calculate time to peak
  hrf_params$shape1 <- hrf_params$a1/hrf_params$b1
  hrf_params$rate1 <- hrf_params$TR/hrf_params$b1
  hrf_params$TRs_to_peak <- (hrf_params$shape1 - 1)/(hrf_params$rate1) #mode of the first Gamma, in terms of seconds
  hrf_params$time_to_peak <- hrf_params$TRs_to_peak * hrf_params$TR #convert to seconds

  hrf_params$a2 <- (16/sqrt(6))*sqrt(hrf_params$a1)*sqrt(hrf_params$b1) # Hardcoded A2
  hrf_params$b2 <- hrf_params$b1 # Hardcoded B2

  #exclude HRFs that peak faster than x
  hrf_params <- hrf_params[hrf_params$time_to_peak >= call_params$time_to_peak_min,] #min timing of response

  a1_grid <- c(min(a1) - call_params$a1_step, a1) + call_params$a1_step / 2
  b1_grid <- c(min(b1) - call_params$b1_step, b1) + call_params$b1_step / 2

  duration_TRs <- 20 # Hardcoded duration_TRs
  inds <- seq(1 / call_params$sr_factor, duration_TRs, by = 1 / call_params$sr_factor) * call_params$TR
  nt <- length(inds)


  hrf_df <- data.frame(a1 = rep(hrf_params$a1, each=nt),
                       b1 = rep(hrf_params$b1, each=nt),
                       TR = TR,
                       sec = inds) #will be recycled


  hrf_df$a2 <- (16/sqrt(6))*sqrt(hrf_df$a1)*sqrt(hrf_df$b1)
  hrf_df$b2 <- hrf_df$b1
  hrf_df$c <- hrf_grid$c
}


my_hrf_grid <- generate_hrf_grid(peak2_time_max = 50, c_vals = c(0, 1/6, 2))
plot_obj <- plot_hrf_param_grid(my_hrf_grid)
ggsave("dev/test_plots/hrf_grid/param_grid_all.pdf", plot_obj, width=10, height=9)




call_params <- attr(default_hrf_grid, "call_params")


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
#' @export
compute_hrf_metrics <- function(hrf_grid) {
  # Extract parameters from attributes
  call_params <- attr(hrf_grid, "call_params")

  # Get parameters from attributes
  TR <- call_params$TR
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
    time_to_end_ii <- inds[max(which(abs(hrf_ii) > 1e-3))]

    # Store in hrf_params
    hrf_params$FWHM[ii] <- FWHM_ii
    hrf_params$time_to_end[ii] <- time_to_end_ii
    hrf_params$peak2_time[ii] <- peak2_time_ii

    # Calculate tapered version if needed
    taper_start <- min(peak2_time_ii, 25)
    if (c_ii > 0 && abs(hrf_ii[which.min(abs(inds - 30))]) > 1e-3) {
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
result <-compute_hrf_metrics(default_hrf_grid)

plot_hrfs_all <- function(hrf_grid) {

  # Compute HRF metrics
  hrf_results <- compute_hrf_metrics(hrf_grid)
  hrf_df <- hrf_results$hrf_df
  hrf_params <- hrf_results$hrf_params

  # Create the plot
  p <- ggplot(hrf_df, aes(x = sec, y = HRF, color = a1, group = a1)) +
    geom_vline(xintercept = c(2, 10), linetype = 2, color = 'gray') +  # limits on time to peak
    geom_hline(yintercept = 0, color = 'gray') +
    geom_line() +
    xlim(0, 30) +
    scale_color_viridis_c(breaks = seq(3, 12, 3)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    facet_grid(b1 ~ c) +
    theme_few() +
    theme(legend.position = 'bottom')

  return(p)
}
my_hrf_grid <- generate_hrf_grid(peak2_time_max = 50)
plot_hrfs_all(my_hrf_grid)


p <- plot_hrfs_all(default_hrf_grid)
ggsave("HRFs_all.pdf", p, width = 8, height = 10.5)
print(p)


#' Plot all tapered HRFs from parameter grid
#'
#' Creates a faceted plot showing all tapered HRF shapes from the parameter grid,
#' organized by b1 values (rows) and undershoot presence (columns).
#' Tapered HRFs are used when the response doesn't resolve to baseline by 30 seconds.
#'
#' @param hrf_grid Data frame with HRF parameters (must have call_params attributes)
#'
#' @return A ggplot object that can be saved or further customized
#'
#' @expor
plot_hrfs_all_tapered <- function(hrf_grid) {

  # Compute HRF metrics
  hrf_results <- compute_hrf_metrics(hrf_grid)
  hrf_df <- hrf_results$hrf_df
  hrf_params <- hrf_results$hrf_params

  # Create the plot
  p <- ggplot(hrf_df, aes(x = sec, y = HRF_tapered, color = a1, group = a1)) +
    geom_vline(xintercept = c(2, 10), linetype = 2, color = 'gray') +  # limits on time to peak
    geom_hline(yintercept = 0, color = 'gray') +
    geom_line() +
    xlim(0, 30) +
    scale_color_viridis_c(breaks = seq(3, 12, 3)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    facet_grid(b1 ~ c) +
    theme_few() +
    theme(legend.position = 'bottom')

  return(p)
}

p <- plot_hrfs_all_tapered(default_hrf_grid)
ggsave("HRFs_tapered.pdf", p, width = 8, height = 10.5)
print(p)

#' Plot HRF parameter grid with computed metrics
#'
#' Creates a two-panel plot showing the HRF parameter space with time-to-peak
#' (top) and FWHM (bottom) values displayed as heatmaps with text labels.
#'
#' @param hrf_grid Data frame with HRF parameters (must have call_params attributes)
#'
#' @return A gridded plot object from grid.arrange that can be saved
#'
#' @export
plot_param_grid_metrics <- function(hrf_grid) {

  # Compute HRF metrics
  hrf_results <- compute_hrf_metrics(hrf_grid)
  hrf_params <- hrf_results$hrf_params

  # Create c_label for faceting
  hrf_params$c_label <- ifelse(hrf_params$c == 0,
                               'No Undershoot (c=0)',
                               paste0('With Undershoot (c=', round(hrf_params$c, 3), ')'))

  # Get unique a1 and b1 values for grid lines
  a1_vals <- sort(unique(hrf_params$a1))
  b1_vals <- sort(unique(hrf_params$b1))

  # Calculate grid line positions (between values)
  a1_grid <- c(min(a1_vals) - 0.5, a1_vals + 0.5)
  b1_grid <- c(min(b1_vals) - 0.125, b1_vals + 0.125)

  # Create time-to-peak plot
  p1 <- ggplot(hrf_params, aes(x = a1, y = b1, fill = time_to_peak)) +
    geom_hline(yintercept = b1_grid, alpha = 0.1) +
    geom_vline(xintercept = a1_grid, alpha = 0.1) +
    geom_tile(alpha = 0.8) +
    geom_text(aes(label = round(time_to_peak, 1))) +
    scale_fill_viridis_c('Time to Peak  ', option = 'C') +
    scale_x_continuous(breaks = a1_vals) +
    scale_y_continuous(breaks = b1_vals) +
    facet_grid(. ~ c_label) +
    theme_few()

  # Create FWHM plot
  p2 <- ggplot(hrf_params, aes(x = a1, y = b1, fill = FWHM)) +
    geom_hline(yintercept = b1_grid, alpha = 0.1) +
    geom_vline(xintercept = a1_grid, alpha = 0.1) +
    geom_tile(alpha = 0.8) +
    geom_text(aes(label = round(FWHM, 1))) +
    scale_fill_viridis_c('Width (FWHM)', option = 'C') +
    scale_x_continuous(breaks = a1_vals) +
    scale_y_continuous(breaks = b1_vals) +
    facet_grid(. ~ c_label) +
    theme_few()

  # Combine plots
  combined_plot <- gridExtra::grid.arrange(p1, p2, nrow = 2)

  return(combined_plot)
}

p <- plot_param_grid_metrics(my_hrf_grid)
ggsave("dev/test_plots/hrf_grid/param_grid_all.pdf", p, width = 10, height = 9)
print(p)
