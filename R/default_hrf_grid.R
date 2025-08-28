#' Default HRF Parameter Grid
#'
#' A pre-computed grid of hemodynamic response function (HRF) parameters
#' for systematic HRF modeling in fMRI data analysis.
#'
#' @format A data frame with 101 rows and 13 columns:
#' \describe{
#'   \item{a1}{Delay of response (seconds). Range: 3-12}
#'   \item{b1}{Response dispersion (seconds). Range: 0.5-2}  
#'   \item{TR}{Temporal resolution (seconds). Value: 2}
#'   \item{shape1}{Shape parameter for first gamma (a1/b1). Range: 2-24}
#'   \item{rate1}{Rate parameter for first gamma (TR/b1). Range: 1-4}
#'   \item{TRs_to_peak}{Time to peak in TR units. Range: 1-5.75}
#'   \item{time_to_peak}{Time to HRF peak (seconds). Range: 2-11.5}
#'   \item{a2}{Delay of undershoot (seconds). Range: 5.66-45.25}
#'   \item{b2}{Dispersion of undershoot (seconds). Equal to b1. Range: 0.5-2}
#'   \item{c}{Scale of undershoot. Values: 0 (no undershoot) or 0.167 (1/6, canonical)}
#'   \item{time_to_end}{Time when HRF resolves (seconds). Range: 6.32-29.34}
#'   \item{FWHM}{Full width at half maximum (seconds). Range: 2.28-10.6}
#'   \item{c_label}{Factor describing undershoot: "No Undershoot (c=0)" or "With Undershoot (c=1/6)"}
#' }
#'
#' @details
#' This grid was constructed by systematically varying HRF parameters around
#' the SPM canonical HRF (a1=6, b1=1, c=1/6) and filtering for physiologically
#' plausible responses. The grid includes 66 HRFs without undershoot (c=0) 
#' and 35 HRFs with canonical undershoot (c=1/6).
#'
#' All calculations assume TR=2 seconds. Parameters were derived using the 
#' double-gamma HRF model and filtered to exclude responses with unrealistic
#' timing characteristics.
#'
#' @source Generated from systematic parameter exploration based on 
#' double-gamma HRF model with physiological constraints
#'
#' @examples
#' # View the grid
#' head(default_hrf_grid)
#' 
#' # Basic statistics
#' summary(default_hrf_grid)
#' 
#' # Check distribution of undershoot types
#' table(default_hrf_grid$c_label)
#' 
#' # Get canonical HRF parameters (a1=6, b1=1, c=1/6)
#' canonical <- default_hrf_grid[default_hrf_grid$a1 == 6 & 
#'                               default_hrf_grid$b1 == 1 & 
#'                               abs(default_hrf_grid$c - 1/6) < 1e-6, ]
#' print(canonical)
#' 
#' # Get fast responses (peak time < 5 seconds)
#' fast_hrfs <- default_hrf_grid[default_hrf_grid$time_to_peak < 5, ]
#' nrow(fast_hrfs)
#' 
#' # Get responses without undershoot only
#' no_undershoot <- default_hrf_grid[default_hrf_grid$c == 0, ]
#' nrow(no_undershoot)  # Should be 66
#' 
#' # Plot parameter space (if ggplot2 available)
#' \dontrun{
#' if(require(ggplot2)) {
#'   ggplot(default_hrf_grid, aes(x = a1, y = b1, color = c_label)) +
#'     geom_point(size = 2) +
#'     labs(title = "HRF Parameter Grid", 
#'          x = "a1 (delay)", y = "b1 (dispersion)", color = "Undershoot") +
#'     theme_minimal()
#' }
#' }
#' 
#' # Use with fit_allHRFs (when available)
#' \dontrun{
#' result <- fit_allHRFs(BOLD, EVs, TR = 0.72, hrf_grid = default_hrf_grid)
#' }
"default_hrf_grid"

#' Generate HRF parameter grid
#'
#' Creates a grid of HRF parameter combinations following the methodology
#' from the research pipeline. Filters out physiologically implausible
#' combinations based on timing constraints.
#'
#' @inheritParams TR_Param
#' @param a1_range Numeric vector of length 2. Min/max values for a1 parameter.
#'   Default spans canonical value ±2x: c(3, 12).
#' @param a1_step Numeric. Step size for a1 sequence. Default: 1.
#' @param b1_range Numeric vector of length 2. Min/max values for b1 parameter.
#'   Default spans canonical value ±2x: c(0.5, 2.0).
#' @param b1_step Numeric. Step size for b1 sequence. Default: 0.25.
#' @param c_vals Numeric vector. Undershoot scaling values to include.
#'   Default: c(1/6, 0) for with/without undershoot.
#' @param sr_factor Integer. Super-resolution factor for HRF timing calculations.
#'   Default: 100.
#' @param time_to_peak_min Numeric. Minimum allowed time-to-peak in seconds.
#'   Default: 2.
#' @param peak2_time_max Numeric. Maximum allowed undershoot peak time in seconds.
#'   Default: 22.74.
#'
#' @return Data frame with columns a1, b1, c, a2, b2 containing valid
#'   HRF parameter combinations.
#'
#' @export
generate_default_hrf_grid <- function(
    TR               = 2,
    a1_range         = c(a1_spm/2, a1_spm*2),
    a1_step          = 1,
    b1_range         = c(b1_spm/2, b1_spm*2),
    b1_step          = 0.25,
    c_vals           = c(1/6, 0),
    sr_factor        = 100,
    time_to_peak_min = 2,
    peak2_time_max   = 22.74
) {
  # canonical SPM defaults
  a1_spm <- 6
  b1_spm <- 1

  # 1) raw grid
  a1_seq <- seq(a1_range[1], a1_range[2], by = a1_step)
  b1_seq <- seq(b1_range[1], b1_range[2], by = b1_step)
  g <- expand.grid(
    a1 = a1_seq,
    b1 = b1_seq,
    TR = TR,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  # 2) analytic time_to_peak
  g$time_to_peak <- ((g$a1 / g$b1) - 1) / (g$TR / g$b1) * g$TR

  # 3) filter too-fast
  g <- subset(g, time_to_peak >= time_to_peak_min)

  # 4) keep only TR == 2 (as before)
  g <- subset(g, TR == 2)

  # 5) compute a2, b2
  g$a2 <- (16 / sqrt(6)) * sqrt(g$a1) * sqrt(g$b1)
  g$b2 <- g$b1

  # 6) add c variants
  g <- do.call(rbind, lapply(c_vals, function(cc) transform(g, c = cc)))

  # 7) super-res time axis
  inds <- seq(1/sr_factor, 20, length.out = 20*sr_factor) * unique(g$TR)

  # 8) compute peak2_time
  g$peak2_time <- NA_real_
  for (i in seq_len(nrow(g))) {
    hr <- hrf::HRF_calc(
      t     = inds,
      deriv = 0,
      a1    = g$a1[i],
      b1    = g$b1[i],
      a2    = g$a2[i],
      c     = g$c[i]
    )
    g$peak2_time[i] <- if (g$c[i] > 0) inds[which.min(hr)] else inds[which.max(hr)]
  }

  # 9) final undershoot-peak filter
  g <- subset(g, peak2_time <= peak2_time_max)

  # 10) return only your five columns, reset rownames
  rownames(g) <- NULL
  g[, c("a1","b1","c","a2","b2"), drop = FALSE]
}

#' Resolve HRF grid specification to actual grid
#'
#' Handles both function and data frame inputs for HRF grid specification.
#' Validates that result contains required columns with appropriate data types.
#'
#' @param hrf_grid Function returning HRF grid data frame, or data frame directly.
#' @param ... Additional arguments passed to hrf_grid if it's a function.
#'
#' @return Data frame with exactly 5 columns: a1, b1, c, a2, b2.
#'
#' @keywords internal
set_hrf_grid <- function(hrf_grid, ...) {
  # If a function, call it. If a data.frame, just use it.
  if (is.function(hrf_grid)) {
    g <- hrf_grid(...)
  } else if (is.data.frame(hrf_grid)) {
    g <- hrf_grid
  } else {
    stop("hrf_grid must be either a data.frame or a function returning a data.frame")
  }

  stopifnot("`hrf_grid` must have at least 5 columns (a1,b1,c,a2,b2)" = ncol(g) >= 5)

  # extract first five columns, force names
  result <- as.data.frame(g)[ , seq_len(5), drop = FALSE]
  names(result) <- c("a1","b1","c","a2","b2")

  stopifnot("columns a1, b1, c, a2, b2 must all be numeric" = all(vapply(result, is.numeric, TRUE)))
  stopifnot("column c (undershoot) must be >= 0" = all(result$c >= 0, na.rm = TRUE))

  result
}
