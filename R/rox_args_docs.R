#' Connectome Workbench
#'
#' @section Connectome Workbench Requirement:
#'
#'  This function uses a system wrapper for the 'wb_command' executable. The
#'  user must first download and install the Connectome Workbench, available
#'  from https://www.humanconnectome.org/software/get-connectome-workbench .
#'
#' @name Connectome_Workbench_Description
NULL

#' aic
#'
#' @param aic (For prewhitening) Use the Akaike information criterion (AIC) to
#'  select AR model orders between \code{0} and \code{ar_order}? Default:
#'  \code{FALSE}.
#'
#' @name aic_Param
NULL

#' ar_order
#'
#' @param ar_order (For prewhitening) The order of the autoregressive (AR) model
#'  to use for prewhitening. If \code{0}, do not prewhiten. Default: \code{6}.
#'
#'  For multi-session modeling, note that a single AR model is used; its
#'  coefficients will be the average estimate from each session.
#'
#' @name ar_order_Param
NULL

#' ar_smooth
#'
#' @param ar_smooth (For prewhitening) The FWHM parameter for spatially
#'  smoothing the coefficient estimates for the AR model to use for
#'  prewhitening. Recall that
#'  \eqn{\sigma = \frac{FWHM}{2*sqrt(2*log(2)}}. Set to \code{0} to not smooth
#'  the estimates. Default: \code{5}.
#'
# [TO DO] vol vs surf?
#' @name ar_smooth_Param
NULL

#'  faces
#'
#' @param faces An \eqn{F \times 3} matrix, where each row contains the vertex
#'  indices for a given triangular face in the mesh. \eqn{F} is the number of
#'  faces in the mesh.
#'
#' @name faces_Param
NULL

#' mask: vertices
#'
#' @param mask  A length \eqn{V} logical vector indicating if each vertex is
#'  within the input mask.
#'
#' @name mask_Param_vertices
NULL

#' mean and variance tolerance
#'
#' @param meanTol,varTol Tolerance for mean and variance of each data location.
#'  Locations which do not meet these thresholds are masked out of the analysis.
#'  Default: \code{1e-6} for both.
#'
#' @name mean_var_Tol_Param
NULL

#' scale_BOLD
#'
#' @param scale_BOLD Controls scaling the BOLD response at each location.
#'  \describe{
#'    \item{"mean":}{   Scale the data to percent local signal change.}
#'    \item{"sd":}{   Scale the data by local standard deviation.}
#'    \item{"none":}{   Center the data but do not scale it.}
#' }
#'
#' @name scale_BOLD_Param
NULL

#' session_names
#'
#' @param session_names The names of the task-fMRI \code{BOLD} sessions, for
#'  multi-session analysis. If not provided here, will be inferred from
#'  \code{names(BOLD)}, inferred from \code{names(design)}, or generated
#'  automatically, in that order.
#'
#' @name session_names_Param
NULL

#' field_names
#'
#' @param field_names (Optional) Names of fields represented in design matrix.
#'
#' @name field_names_Param
NULL

#' verbose
#'
#' @param verbose \code{1} (default) to print occasional updates during model
#'  computation; \code{2} for occasional updates as well as running INLA in
#'  verbose mode (if \code{Bayes}), or \code{0} for no printed updates.
#'
#' @name verbose_Param
NULL
