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

#' scrub
#'
#' @param scrub (Optional) A \eqn{T \times N_{scrub}} matrix of spike regressors
#'  (one 1 value at the timepoint to scrub, and 0 for all other values), or a
#'  logical vector indicating the timepoints to scrub (\code{TRUE} to scrub, and
#'  \code{FALSE} to keep).
#'
#'  The spike regressors will be included in the nuisance
#'  regression, and afterwards the timepoints indicated in \code{scrub} will be
#'  removed from the BOLD data and design matrix.
#'
#' @name scrub_Param
NULL


#' nT
#'
#' @param nT The number of timepoints in the dataset. Used to validate that
#'  scrub vectors, spike regressors, and nuisance regressors have consistent
#'  dimensions.
#'
#' @name nT_Param
NULL

#' spike_matrix
#'
#' @param spike_matrix A binary \eqn{T \times N} spike regressor matrix, where
#'  \eqn{T} is the number of timepoints and \eqn{N} is the number of spike
#'  regressors (one per scrubbed timepoint). Each column must contain exactly
#'  one `1` and \eqn{T-1} zeros.
#'
#' @name spike_matrix_Param
NULL

#' onset_skip
#'
#' @param onset_skip Minimum event duration (in seconds) for tasks to be
#'  included in onset/offset regressors. Tasks where ANY event has duration less than
#'  the threshold will be excluded from the respective onset/offset regressors. This is
#'  useful for excluding short events (e.g., cues) from block-based onset/offset modeling.
#'  Default: 5 seconds. Set to \code{NULL} or 0 to include all tasks.
#'
#' @name onset_skip_Param
NULL

#' offset_skip
#'
#' @param offset_skip Minimum event duration (in seconds) for tasks to be
#'  included in onset/offset regressors. Tasks where ANY event has duration less than
#'  the threshold will be excluded from the respective onset/offset regressors. This is
#'  useful for excluding short events (e.g., cues) from block-based onset/offset modeling.
#'  Default: 5 seconds. Set to \code{NULL} or 0 to include all tasks.
#'
#' @name offset_skip_Param
NULL

#' onsets_sep
#'
#' @param onsets_sep Model the onsets separately for each task? Default: \code{FALSE}, to
#'  model all onsets together as a single field in the design.
#'
#' @name onsets_sep_Param
NULL

#' offsets_sep
#'
#' @param offsets_sep Model the offsets separately for each task? Default: \code{FALSE}, to
#'  model all offsets together as a single field in the design.
#'
#' @name offsets_sep_Param
NULL

#' onset
#'
#' @param onset Add task regressors indicating the onset of each event block? 
#'  Provide the names of the tasks as a character vector. All onsets across the 
#'  specified tasks will be represented by one additional column in the design matrix. 
#'  The task names must match the names of \code{EVs}. Can also be \code{TRUE} or 
#'  \code{"all"} to use all tasks.
#'
#'  Onsets/offset modeling is only compatible with a block design experiment.
#'  An error will be raised if the events in \code{EVs} do not have duration
#'  greater than one second.
#'
#' @name onset_Param
NULL

#' offset
#'
#' @param offset Add task regressors indicating the offset of each event block? 
#'  Provide the names of the tasks as a character vector. All offsets across the 
#'  specified tasks will be represented by one additional column in the design matrix. 
#'  The task names must match the names of \code{EVs}. Can also be \code{TRUE} or 
#'  \code{"all"} to use all tasks.
#'
#'  Onsets/offset modeling is only compatible with a block design experiment.
#'  An error will be raised if the events in \code{EVs} do not have duration
#'  greater than one second.
#'
#' @name offset_Param
NULL