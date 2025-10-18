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

#' BOLD Data Input
#'
#' @param BOLD Character vector of file paths to CIFTI files; 
#'   one entry per subject/session. CIFTI files should
#'   contain preprocessed BOLD time-series data.
#'
#' @name BOLD_Param
NULL

#' Event Definitions
#'
#' @param EVs List, same length as \code{BOLD}. Each element is a data frame
#'   (or similar) giving event onsets and durations for that subject. Each
#'   data frame should have columns for onset times and durations in seconds.
#'
#' @name EVs_Param
NULL

#' Nuisance Regressors
#'
#' @param nuisance Character vector of file paths (one per subject) or
#'   \code{NULL}. Each file must contain a design matrix of nuisance
#'   regressors (e.g., motion parameters, physiological noise). If \code{NULL},
#'   no nuisance regression is performed.
#'
#' @name nuisance_Param
NULL

#' Scrubbing Parameters
#'
#' @param scrub List (one logical vector per subject) or \code{NULL}.
#'   Frames set to \code{TRUE} are excluded from the GLM. Used to remove
#'   high-motion or artifact-contaminated timepoints.
#'
#' @name scrub_Param
NULL

#' Temporal Resolution
#'
#' @param TR Numeric. Repetition time in seconds. The temporal sampling
#'   rate of the BOLD acquisition.
#'
#' @name TR_Param
NULL

#' Brain Structures
#'
#' @param brainstructures Character vector passed to
#'   \code{ciftiTools::read_cifti()}, specifying which brain structures
#'   to analyze. Typical values: \code{c("left", "right")} for cortical
#'   hemispheres, or \code{c("left", "right", "subcortical")}.
#'
#' @name brainstructures_Param
NULL

#' Surface Resampling
#'
#' @param resamp_res Integer or \code{NULL}. Surface resampling resolution
#'   (e.g., 10000 for 10k vertices per hemisphere). If \code{NULL}, uses
#'   native resolution of the input data.
#'
#' @name resamp_res_Param
NULL

#' High-Pass Filtering
#'
#' @param hpf Numeric or \code{NULL}. High-pass filter cut-off frequency
#'   in Hz. Set \code{NULL} to disable filtering. Recommended for removing
#'   low-frequency drift in fMRI data.
#'
#' @name hpf_Param
NULL

#' HRF Parameters
#'
#' @param hrf_params Named list with at least \code{a1}, \code{b1}, \code{c}.
#'   Additional parameters \code{a2}, \code{b2} are optional and default to
#'   SPM values. Controls the shape of the hemodynamic response function:
#'   \itemize{
#'     \item \code{a1}: Delay of main response (default: 6)
#'     \item \code{b1}: Dispersion of main response (default: 1)  
#'     \item \code{c}: Scale of undershoot (default: 1/6)
#'     \item \code{a2}: Delay of undershoot (default: 16)
#'     \item \code{b2}: Dispersion of undershoot (default: 1)
#'   }
#'
#' @name hrf_params_Param
NULL

#' HRF Derivatives
#'
#' @param derivatives Logical. Include temporal and dispersion derivatives
#'   of the HRF in the design matrix? This allows for flexibility in HRF
#'   timing and shape across brain regions.
#'
#' @name derivatives_Param
NULL

#' Onset Modeling
#'
#' @param onsets Logical or character vector. Add separate onset regressors
#'   per trial type? If \code{TRUE}, adds onset regressors for all tasks.
#'   If character vector, specifies which tasks to model onsets for.
#'   Only compatible with block designs (events with duration > 1s).
#'
#' @name onsets_Param  
NULL

#' Offset Modeling
#'
#' @param offsets Logical or character vector. Add separate offset regressors
#'   per trial type? If \code{TRUE}, adds offset regressors for all tasks.
#'   If character vector, specifies which tasks to model offsets for.
#'   Only compatible with block designs (events with duration > 1s).
#'
#' @name offsets_Param
NULL

#' Statistical Threshold
#'
#' @param alpha Numeric. P-value threshold used when generating activation
#'   masks across subjects. Typical values: 0.05, 0.01, or 0.001.
#'
#' @name alpha_Param
NULL

#' Verbosity Level
#'
#' @param verbose Integer controlling output verbosity:
#'   \itemize{
#'     \item 0 = silent
#'     \item 1 = high-level progress messages  
#'     \item 2 = detailed per-subject information
#'   }
#'
#' @name verbose_Param
NULL

#' Parallel Processing
#'
#' @param n_cores Positive integer. Number of physical CPU cores to allocate
#'   for parallel processing. If \code{n_cores > 1}, subjects are processed
#'   in parallel using \code{parLapplyLB()}. Each core may spawn multiple
#'   workers to mask I/O latency.
#'
#' @name n_cores_Param
NULL

#' Subject Index
#'
#' @param subject_idx Integer. Index identifying the current subject being
#'   processed (used internally for progress reporting and error tracking).
#'
#' @name subject_idx_Param
NULL

#' Design Matrix Data
#'
#' @param design_matrix Numeric matrix containing the GLM design matrix
#'   (timepoints × regressors). Created from event definitions and HRF
#'   convolution.
#'
#' @name design_matrix_Param
NULL

#' Design Array Format
#'
#' @param design_array 3-dimensional array version of design matrix,
#'   formatted for compatibility with \code{multiGLM()}. Third dimension
#'   contains multiple design variations.
#'
#' @name design_array_Param
NULL

#' Subject Results
#'
#' @param subject_results List containing GLM results for all processed
#'   subjects. Each element contains model fits, statistics, and metadata
#'   for one subject.
#'
#' @name subject_results_Param
NULL

#' BOLD Data Object
#'
#' @param BOLD_xii A \code{xifti} object containing loaded BOLD time-series
#'   data, typically after resampling and preprocessing.
#'
#' @name BOLD_xii_Param
NULL

#' Nuisance File Path
#'
#' @param nuisance_file Character or NULL. File path to nuisance regressors
#'   text file. Expected format: tab/space-separated values, no header, 
#'   one timepoint per row. If NULL, no nuisance regression is performed.
#'
#' @name nuisance_file_Param
NULL

#' Log Directory Path
#'
#' @param log_dir Character or NULL. Directory path for saving log files
#'   and temporary outputs during processing. If NULL, uses current working
#'   directory or system temp directory.
#'
#' @name log_dir_Param
NULL

#' HRF Parameter Grid
#'
#' @param hrf_grid Data frame or function returning HRF parameter combinations.
#'   Must contain columns: a1, b1, c, a2, b2. If function, called with ... arguments.
#'
#' @name hrf_grid_Param
NULL

#' work_dir
#'
#' @param work_dir Character. Directory path where temporary result files will be
#'  saved during processing. If relative, interpreted relative to current working
#'  directory. If the directory doesn't exist, it will be created. Default: "work".
#'  Large intermediate results are saved here to manage memory usage, especially
#'  in parallel processing mode.
#'
#' @name work_dir_Param
NULL

#' Minimum Active Subjects Threshold
#'
#' @param min_active_subjects Integer. Minimum number of subjects that
#'   must show significant activation (p < alpha) at a voxel for it to be
#'   included in the group mask used for regularization. Default is 20.
#'   If the number of successful subjects is less than this value, it will
#'   be automatically adjusted with a warning.
#'
#' @name min_active_subjects_Param
NULL