#' Default HRF Parameter Grid
#'
#' A pre-computed grid of hemodynamic response function (HRF) parameters
#' for systematic HRF modeling in fMRI data analysis.
#'
#' @format A data frame with 101 rows and 5 columns:
#' \describe{
#'   \item{a1}{Delay of response (seconds). Range: 3–12}
#'   \item{b1}{Response dispersion (seconds). Range: 0.5–2}
#'   \item{c}{Undershoot scale. Values: 0 (no undershoot) or 1/6 (canonical)}
#'   \item{a2}{Delay of undershoot (seconds). Calculated from a1 and b1}
#'   \item{b2}{Dispersion of undershoot (seconds). Equal to b1}
#' }
#'
#' @details
#' This grid was constructed by systematically varying HRF parameters around
#' the SPM canonical HRF (a1=6, b1=1, c=1/6) and filtering for physiologically
#' plausible responses. All calculations assume TR = 2 seconds.
#'
#' @source Generated from systematic parameter exploration based on
#' double-gamma HRF model with physiological constraints
#'
#' @examples
#' head(default_hrf_grid)
#' summary(default_hrf_grid)
#' table(default_hrf_grid$c)
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
generate_hrf_grid <- function(
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
  result <- g[, c("a1","b1","c","a2","b2"), drop = FALSE]

  # attach call params so downstream functions can use them
  attr(result, "call_params") <- list(
    TR = TR,
    a1_range = a1_range,
    a1_step = a1_step,
    b1_range = b1_range,
    b1_step = b1_step,
    c_vals = c_vals,
    sr_factor = sr_factor,
    time_to_peak_min = time_to_peak_min,
    peak2_time_max = peak2_time_max
  )

  class(result) <- c("hrf_grid", "data.frame")
  return(result)
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
