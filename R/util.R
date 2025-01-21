#' Mask out invalid data
#'
#' Mask out data locations that are invalid (missing data, low mean, or low
#'  variance) for any session.
#'
#' @param BOLD A session-length list of \eqn{T \times V} numeric BOLD data.
#' @param meanTol,varTol Tolerance for mean and variance of each data location.
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Defaults: \code{1e-6}.
#' @param verbose Print messages counting how many locations are removed?
#'  Default: \code{TRUE}.
#'
#' @importFrom matrixStats colVars
#' @return A logical vector indicating locations that are valid across all sessions.
#'
#' @examples
#' nT <- 30
#' nV <- 400
#' BOLD1 <- matrix(rnorm(nT*nV), nrow=nT)
#' BOLD1[,seq(30,50)] <- NA
#' BOLD2 <- matrix(rnorm(nT*nV), nrow=nT)
#' BOLD2[,65] <- BOLD2[,65] / 1e10
#' BOLD <- list(sess1=BOLD1, sess2=BOLD2)
#' do_QC(BOLD)
#'
#' @export
do_QC <- function(BOLD, meanTol=1e-6, varTol=1e-6, verbose=TRUE){

  nS <- length(BOLD)
  nV <- ncol(BOLD[[1]])

  # For each BOLD data matrix,
  mask_na <- mask_mean <- mask_var <- mask_snr <- rep(TRUE, nV)
  for (ss in seq(nS)) {
    # Mark columns with any NA or NaN values for removal.
    na_ss <- is.na(BOLD[[ss]]) | is.nan(BOLD[[ss]])
    mask_na[apply(na_ss, 2, any)] <- FALSE
    # Calculate means and variances of columns, except those with any NA or NaN.
    # Mark columns with mean/var falling under the thresholds for removal.
    means_ss <- colMeans(BOLD[[ss]][,mask_na,drop=FALSE])
    vars_ss <- matrixStats::colVars(BOLD[[ss]][,mask_na,drop=FALSE])
    snr_ss <- rep(NA, nV) # use `NA` values for `NA` columns.
    snr_ss[mask_na] <- means_ss/sqrt(vars_ss)
    mask_mean[mask_na][means_ss < meanTol] <- FALSE
    mask_var[mask_na][vars_ss < varTol] <- FALSE
    #mask_snr[mask_na][snr_ss < snrTol] <- FALSE
  }

  # Print counts of locations removed, for each reason.
  if (verbose) {
    warn_part1 <- " locations"
    warn_part2 <- if (nS > 1) { " in at least one session.\n" } else { ".\n" }
    if (any(!mask_na)) {
      cat(paste0(
        "\t", sum(!mask_na), warn_part1,
        " removed due to NA/NaN values", warn_part2
      ))
      warn_part1 <- " additional locations"
    }
    # Do not include NA locations in count.
    mask_mean2 <- mask_mean | (!mask_na)
    if (any(!mask_mean2)) {
      cat(paste0(
        "\t", sum(!mask_mean2), warn_part1,
        " removed due to low mean", warn_part2
      ))
      warn_part1 <- " additional locations"
    }
    # Do not include NA or low-mean locations in count.
    mask_var2 <- mask_var | (!mask_mean) | (!mask_na)
    if (any(!mask_var2)) {
      cat(paste0(
        "\t", sum(!mask_var2), warn_part1,
        " removed due to low variance", warn_part2
      ))
    }
    # # Do not include NA or low-mean or low-var locations in count.
    # mask_snr2 <- mask_snr | (!mask_mean) | (!mask_var) | (!mask_na)
    # if (any(!mask_snr2)) {
    #   cat(paste0(
    #     "\t", sum(!mask_snr2), warn_part1,
    #     " removed due to low SNR", warn_part2
    #   ))
    # }
  }

  # Return composite mask and other masks
  list(
    mask =  mask_na & mask_mean & mask_var, # & mask_snr,
    mask_na = mask_na,
    mask_mean = mask_mean,
    mask_var = mask_var,
    mask_snr = snr_ss #just return the actual SNR values for now
  )
}

#' Is a matrix or data.frame?
#'
#' Is this a matrix or data.frame?
#'
#' @param x The object
#' @return Length-one logical.
#' @keywords internal
is_matrix_or_df <- function(x){
  is.matrix(x) || is.data.frame(x)
}

#' `cbind` if first argument might be \code{NULL}
#'
#' `cbind`, but return the second argument if the first is \code{NULL}
#' @param mat_or_NULL \code{NULL} or a numeric matrix
#' @param to_add A numeric matrix with the same number of rows as \code{mat_or_NULL}
#' @return \code{cbind(mat_or_NULL, to_add)}, or just \code{to_add} if the first argument is NULL.
#' @keywords internal
cbind2 <- function(mat_or_NULL, to_add) {
  if (!is.null(mat_or_NULL)) {
    cbind(mat_or_NULL, to_add)
  } else {
    to_add
  }
}
