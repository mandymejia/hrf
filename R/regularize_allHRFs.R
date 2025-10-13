#' Regularize HRF Parameters Across Subjects
#'
#' Performs spatial regularization of hemodynamic response function (HRF)
#' parameters by combining activation masks from working HRF analysis with
#' best-fitting parameters from all-HRFs analysis. The function applies both
#' ordinary least squares (OLS) and weighted least squares (WLS) methods to
#' estimate subject-level deviations from population averages.
#'
#' @param workingHRF_results Results object from \code{fit_workingHRF()},
#'   containing activation masks that indicate which voxels showed significant
#'   activation for each subject using the working HRF model.
#' @param allHRF_results Results object from \code{fit_allHRFs()}, containing
#'   best-fitting HRF parameters (a1, b1, c) for each voxel-subject combination
#'   from fitting multiple HRF variants.
#' @param rounding Logical; if TRUE (default), round the regularized HRF parameters
#'   to the nearest values in the provided HRF grid for better interpretability.
#'
#' @return An object of class \code{"regularizeHRFs"} containing:
#'   \describe{
#'     \item{best_params_df}{Combined dataframe with HRF parameters and
#'       activation masks for all subjects and voxels.}
#'     \item{regularized_params}{Nested list structure with regularization
#'       results for each parameter (a1, b1, c), including population averages,
#'       subject-level estimates, residual variances, and brain visualization
#'       objects for both OLS and WLS methods.}
#'     \item{rounded_params}{Rounded regularized HRF parameters (if
#'       \code{rounding = TRUE}), aligned to the nearest values in the HRF grid.}
#'     \item{mask_prop_NA}{Logical vector mask indicating voxels with sufficient
#'       activation consistency across subjects (values outside threshold are set
#'       to \code{NA}). Useful for filtering voxels in downstream analyses.}
#'   }
#'
#' @details
#' The regularization process consists of several steps:
#' \enumerate{
#'   \item Validates input dimensions and extracts brain template
#'   \item Merges activation masks with parameter estimates
#'   \item For each parameter (a1, b1, c):
#'     \itemize{
#'       \item Computes population-averaged parameter estimates
#'       \item Estimates subject-level slopes and intercepts using OLS
#'       \item Calculates residual variances for uncertainty quantification
#'       \item Re-estimates parameters using WLS with inverse variance weights
#'       \item Creates brain visualization objects for all results
#'     }
#' }
#'
#' The function handles missing data through activation masking and provides
#' both unweighted (OLS) and precision-weighted (WLS) estimates to account
#' for heteroscedastic variance across brain regions.
#'
#' @seealso \code{\link{fit_workingHRF}}, \code{\link{fit_allHRFs}},
#'   \code{\link{HRF_regularize}}
#'
#' @examples
#' \dontrun{
#' # Assuming you have results from previous steps
#' working_results <- fit_workingHRF(BOLD, EVs, TR = 0.72)
#' all_results <- fit_allHRFs(BOLD, EVs, TR = 0.72, hrf_grid = my_grid)
#'
#' # Regularize parameters
#' regularized <- regularize_allHRFs(working_results, all_results)
#'
#' }
#'
#' @export
regularize_allHRFs <- function(workingHRF_results, allHRF_results, rounding = TRUE) {
  cat("****************Version 0.1.11********************\n")

  xii <- validate_previous_results(workingHRF_results, allHRF_results)
  result <- create_best_params_df(workingHRF_results, allHRF_results)

  best_params_df <- result$best_params_df
  mask_prop_NA <- result$mask_prop_NA

  hrf_grid <- allHRF_results$hrf_grid
  rm(workingHRF_results, allHRF_results); gc()
  regularized_params <- HRF_regularize(best_params_df, mask_prop_NA, xii, mask = TRUE,
                                       log = TRUE, truncate = TRUE, WLS = TRUE, 
                                       use_subject_activation_mask = TRUE)

  if (rounding) {
    rounded <-round_regularized_params(regularized_params,
                                       hrf_grid, model = "A", inside_grid = TRUE)
  }

  result <- list (
    best_params_df = best_params_df,
    regularized_params = regularized_params,
    rounded_params = rounded,
    mask_prop_NA = mask_prop_NA
  )

  class(result) <- "regularizeHRFs"
  attr(result, "xii") <- xii
  return(result)

}

#' Create Best Parameters Dataframe for HRF Regularization
#'
#' This function merges activation masks from working HRF analysis with best-fitting
#' HRF parameters from all-HRFs analysis to create a consolidated dataframe for
#' regularization. It prepares the data structure needed for the regularization
#' step by combining subject-specific activation information with overfitted
#' parameter estimates.
#'
#' @param workingHRF_results Results object from \code{fit_workingHRF()}, containing
#'   activation masks that indicate which voxels showed significant activation
#'   for each subject using the working HRF model.
#' @param allHRF_results Results object from \code{fit_allHRFs()}, containing
#'   best-fitting HRF parameters (a1, b1, c) for each voxel-subject combination
#'   from fitting multiple HRF variants.
#'
#' @return A dataframe with the following columns:
#'   \describe{
#'     \item{a1, b1, c}{HRF parameters (numeric) from the best-fitting model}
#'     \item{voxel}{Voxel/location identifier (integer)}
#'     \item{subject}{Subject identifier (character)}
#'     \item{mask}{Logical indicating significant activation for this voxel-subject pair}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts activation masks from working HRF results
#'   \item Reshapes masks from wide to long format using \code{reshape2::melt()}
#'   \item Extracts best-fitting parameters from all-HRFs results
#'   \item Converts parameters to factors then back to numeric (matching research pipeline)
#'   \item Merges masks with parameters using \code{dplyr::left_join()}
#' }
#'
#' The resulting dataframe contains both the overfitted parameter estimates and
#' activation information needed for spatial regularization of HRF parameters.
#'
#' @keywords internal
create_best_params_df <- function(workingHRF_results, allHRF_results ) {
  masks <- workingHRF_results[["activation_masks"]][["masks"]]
  masks_df <- masks_df <- as.data.frame(masks) # Reshape dataframe
  masks_df$voxel <- 1:nrow(masks_df)
  masks_df <- reshape2::melt(masks_df, id.vars = 'voxel', variable.name = 'subject', value.name = 'mask')

  masks_df_prop <- masks_df %>% group_by(voxel) %>% summarize(prop = mean(mask, na.rm=TRUE))
  
  subj_prop <- get_subj_prop(workingHRF_results)
  mask_prop <- (masks_df_prop$prop > subj_prop)   # mask_prop <- (masks_df_prop$prop > 0.1)
  
  mask_prop_NA <- mask_prop; mask_prop_NA[!mask_prop_NA] <- NA


  # At this point research code saves mask_prop_NA
  #################################################

  # best_params_df creation
  best_params_df <- allHRF_results[["best_params_results"]]
  a1_vals <- sort(unique(best_params_df$a1))
  b1_vals <- sort(unique(best_params_df$b1))
  c_vals <- sort(unique(best_params_df$c))
  best_params_df$a1 <- factor(best_params_df$a1, levels=a1_vals)
  best_params_df$b1 <- factor(best_params_df$b1, levels=b1_vals)
  best_params_df$c <- factor(best_params_df$c, levels=c_vals)

  best_params_df$subject <- as.character(best_params_df$subject) # Convert to left join
  masks_df$subject <- as.character(masks_df$subject)


  #1. merge masks into best_params_df
  best_params_df <- left_join(best_params_df, masks_df)

  best_params_df$a1 <- as.numeric(as.character(best_params_df$a1))
  best_params_df$b1 <- as.numeric(as.character(best_params_df$b1))
  best_params_df$c <- as.numeric(as.character(best_params_df$c))

  return(list(
    best_params_df = best_params_df,
    mask_prop_NA = mask_prop_NA
  ))

}

#' Calculate the Proportion of Active Subjects
#'
#' Computes the proportion of subjects considered "active" based on the
#' `min_active_subjects` parameter in the `workingHRF_results` object.
#' If missing, defaults to the greater of 2 or 10% of total subjects.
#'
#' @param workingHRF_results HRF analysis results list.
#' @return Numeric proportion of active subjects (rounded to two decimals).
#' @keywords internal
get_subj_prop <- function(workingHRF_results) {
  min_active_subjects <- workingHRF_results[["call_info"]][["min_active_subjects"]]
  if (is.null(min_active_subjects)) {
    total_subjects <- sum(vapply(workingHRF_results[["subject_results"]], \(x) identical(x[["status"]], "success"), logical(1)))
    min_active_subjects <- max(2, ceiling(total_subjects * 0.10))
    message("min_active_subjects not found in workingHRF_results (old version), defaulting to max(10% of total, 2): ", min_active_subjects, " subjects")
  } 
  subj_prop <- round(min_active_subjects / length(workingHRF_results[["subject_results"]]), 2)
  message("Using group mask threshold: ", subj_prop, " (", min_active_subjects, " out of ", length(workingHRF_results[["subject_results"]]), " subjects)")
  return(subj_prop)
}

#' Regularize HRF Parameters Using Hierarchical Modeling
#'
#' Performs hierarchical regularization of HRF parameters by modeling
#' subject-specific deviations from population averages. The function
#' implements both ordinary least squares (OLS) and weighted least squares
#' (WLS) approaches to estimate subject-level parameters while accounting
#' for spatial heterogeneity in estimation uncertainty.
#'
#' @param best_params_df Dataframe containing HRF parameter estimates with
#'   columns for parameters (a1, b1, c), voxel identifiers, subject
#'   identifiers, and activation masks.
#' @param xii Template xifti object used for creating brain visualizations
#'   with proper anatomical structure and dimensions.
#' @param mask Logical. If \code{TRUE} (default), only analyze voxels that
#'   showed significant activation in the working HRF analysis. If \code{FALSE},
#'   analyze all voxels.
#' @param log Logical. Reserved for future log-transformation of parameters.
#'   Currently not implemented.
#' @param truncate Logical. If \code{TRUE} (default), subjects with negative
#'   slopes are identified and their slopes are set to zero, with intercepts
#'   re-estimated accordingly.
#' @param WLS Logical. If \code{TRUE} (default), performs weighted least
#'   squares estimation using inverse variance weights after initial OLS
#'   estimation. If \code{FALSE}, only performs OLS estimation.
#' @param use_subject_activation_mask Logical. If \code{TRUE} (default), only
#'   includes subjects with activation at each voxel when computing population
#'   averages (i.e., filters by \code{best_params_df$mask == TRUE}). If
#'   \code{FALSE}, includes all subjects regardless of individual activation
#'   status.
#'
#' @return A nested list with one entry per parameter (a1, b1, c), where each
#'   parameter's results contain:
#'   \describe{
#'     \item{pop_mean_xii}{Brain map of population-averaged parameter values}
#'     \item{best_params_df_est_OLS}{Subject-level OLS estimates (slopes,
#'       intercepts, offsets)}
#'     \item{results_OLS}{Full dataframe with OLS predictions, residuals,
#'       and variances}
#'     \item{residual_variance_xii_A_OLS, residual_variance_xii_B_OLS}{Brain
#'       maps of OLS residual variances for Model A and B}
#'     \item{best_params_df_est_WLS}{Subject-level WLS estimates (if WLS=TRUE)}
#'     \item{results_WLS}{Full dataframe with WLS predictions, residuals,
#'       and variances (if WLS=TRUE)}
#'     \item{residual_variance_xii_A_WLS, residual_variance_xii_B_WLS}{Brain
#'       maps of WLS residual variances (if WLS=TRUE)}
#'   }
#'
#' @details
#' The regularization process for each parameter involves:
#' \enumerate{
#'   \item \strong{Population Estimation}: Computing voxel-wise averages
#'     across subjects for activated regions
#'   \item \strong{OLS Subject Modeling}: For each subject, fitting linear
#'     models relating individual estimates to population averages:
#'     \itemize{
#'       \item Model A: \eqn{param_{ij} = intercept_i + slope_i \times pop\_avg_j}
#'       \item Model B: \eqn{param_{ij} = offset_i + pop\_avg_j}
#'     }
#'   \item \strong{Variance Estimation}: Computing residual variances to
#'     quantify estimation uncertainty
#'   \item \strong{WLS Re-estimation}: Using inverse variance weights to
#'     improve parameter estimates in heteroscedastic regions
#' }
#'
#' Negative slope truncation prevents unrealistic parameter relationships where
#' subjects systematically deviate in the opposite direction from population trends.
#'
#' @keywords internal
HRF_regularize <- function(best_params_df, mask_prop_NA, xii, mask = TRUE,
                           log = TRUE, truncate = TRUE, WLS = TRUE, 
                           use_subject_activation_mask = TRUE  ) {
  param_names <-  extract_param_names(best_params_df)
  model_outputs <- list() # To store the results for each param

  best_params_df$use <- if(!use_subject_activation_mask) TRUE else best_params_df$mask

  for(pp in param_names) {
    print(paste("Running analysis for param", pp, "..."))

    pop_results <- estimate_population_parameters(best_params_df, mask_prop_NA, pp, xii, mask)
    best_params_df <- pop_results$best_params_df
    mean0_xii <- pop_results$mean0_xii

    # Initialize an empty list and save pop_mean_xii as function output
    model_outputs[[pp]] <- list()
    model_outputs[[pp]]$pop_mean_xii <- mean0_xii

    # OLS estimation and variance calculation
    ols_results <- estimate_ols_and_variance(best_params_df, truncate, xii)
    best_params_df <- ols_results$best_params_df
    model_outputs[[pp]]$best_params_df_est_OLS <- ols_results$best_params_df_est
    model_outputs[[pp]]$results_OLS <- ols_results$best_params_df
    model_outputs[[pp]]$residual_variance_xii_A_OLS <- ols_results$variance_xii_A
    model_outputs[[pp]]$residual_variance_xii_B_OLS <- ols_results$variance_xii_B

    print("Completed OLS")

    if(WLS) {

      wls_results <- estimate_wls_and_variance(best_params_df, truncate, mean0_xii)
      model_outputs[[pp]]$best_params_df_est_WLS <- wls_results$best_params_df_est
      model_outputs[[pp]]$results_WLS <- wls_results$best_params_df
      model_outputs[[pp]]$residual_variance_xii_A_WLS <- wls_results$variance_xii_A
      model_outputs[[pp]]$residual_variance_xii_B_WLS <- wls_results$variance_xii_B

      print("Completed WLS")
    }
  }
  return(model_outputs)
}

#' Estimate Population-Level HRF Parameters
#'
#' Computes population-averaged HRF parameter estimates by averaging individual
#' subject estimates across voxels. This represents the first step in the HRF
#' regularization process, creating group-level parameter maps and preparing
#' data for subsequent subject-level modeling.
#'
#' @param best_params_df Dataframe containing HRF parameter estimates with columns
#'   for parameters (a1, b1, c), voxel identifiers, subject identifiers, and
#'   activation masks. Modified in place with additional columns.
#' @param pp Character string specifying which parameter column to process
#'   (e.g., "a1", "b1", "c").
#' @param xii Template xifti object used for creating brain visualizations with
#'   proper anatomical structure and dimensions.
#'
#' @return A list containing:
#'   \describe{
#'     \item{best_params_df}{Updated dataframe with new columns: \code{param},
#'       \code{param0}, and \code{param_mean0} joined from population averaging.
#'       Previous iteration columns are cleared.}
#'     \item{mean0_xii}{xifti object containing population-averaged parameter
#'       values properly formatted for brain visualization, with NA padding
#'       for voxels not meeting activation criteria.}
#'   }
#'
#' @details
#' The function performs the following operations:
#' \enumerate{
#'   \item Selects the current parameter and clears previous iteration columns
#'   \item Filters to voxels meeting activation criteria (\code{use} column)
#'   \item Computes voxel-wise averages across subjects
#'   \item Joins population averages back to the main dataframe
#'   \item Creates an xifti visualization object with NA padding for missing voxels
#' }
#'
#' This function modifies the input dataframe by adding parameter-specific columns
#' and clearing columns from previous parameter iterations to ensure a clean state.
#'
#' @keywords internal
estimate_population_parameters <- function(best_params_df, mask_prop_NA, pp, xii, mask) {
  best_params_df$param <- best_params_df[,pp] #select the current parameter
  best_params_df$param0 <- best_params_df$param

  best_params_df$param_mean0 <- NULL # Initiliazation
  best_params_df$param_offset <- best_params_df$param_int <- best_params_df$param_slope <- NULL
  best_params_df$param_pred_A <- best_params_df$param_pred_B <- NULL
  print("Running Ordinary Least Squares Method")

  if(mask) {
    # Use population mask + activation mask
    best_params_df_avg0 <- best_params_df %>%
      filter(use & voxel %in% which(!is.na(mask_prop_NA))) %>%
      group_by(voxel) %>%
      summarize(param_mean0 = mean(param0, na.rm=TRUE))
  } else {
    # Original behavior - all voxels, all subjects
    best_params_df_avg0 <- best_params_df %>%
      filter(use) %>%  # When mask=FALSE, use=TRUE for everyone
      group_by(voxel) %>%
      summarize(param_mean0 = mean(param0, na.rm=TRUE))
  }

  best_params_df <- best_params_df %>% left_join(best_params_df_avg0, by = "voxel") # Bing population means back into full df

  best_params_df <- best_params_df %>%
    filter(!is.na(param_mean0)) # Because mask_prop_NA is applied, so the left join can create NAs

  mean0_xii <- create_xifti_with_padding(best_params_df_avg0, "param_mean0", xii)
  # mean0_xii <- ciftiTools::newdata_xifti(xii, best_params_df_avg0$param_mean0)

  return(list(
    best_params_df = best_params_df,
    mean0_xii = mean0_xii
  ))
}

#' Estimate OLS Parameters and Calculate Variance
#'
#' Performs ordinary least squares estimation for subject-level parameters
#' and calculates residual variances for Model A and Model B predictions.
#'
#' @param best_params_df Dataframe with parameter estimates and population means
#' @param truncate Logical, whether to truncate negative slopes
#' @param xii Template xifti object for creating variance visualizations
#'
#' @return List containing updated dataframe, estimates, and variance objects
#' @keywords internal
estimate_ols_and_variance <- function(best_params_df, truncate, xii) {
  ### NORMAL METHOD IMPLEMENTATION (NON-WLS)
  # This section contains the original non-weighted implementation

  ### Step 2: for each subject, fit a model relating the over-fitted estimates to the population average
  ### Model A: intercept and slope
  ### Model B: offset-only (an intercept-only model, which might be preferable when model A results in slope = 0 or negative slopes)

  best_params_df_est <- best_params_df %>%
    filter(use) %>%
    group_by(subject) %>%
    summarize(param_offset = mean(param - param_mean0),
              param_slope = slope_fun(param, param_mean0, lm=FALSE),
              param_int = mean(param) - param_slope*mean(param_mean0))


  ### Step 3: truncate negative slopes, re-estimate intercepts as the straight mean
  if(truncate) best_params_df_est <- truncate_negative_slopes(best_params_df,
                                                              best_params_df_est,
                                                              method = "OLS")

  # Keep subjects in numerical order
  best_params_df_est <- arrange(best_params_df_est, subject)


  # Bring subject-level slopes and intercepts back into full df
  best_params_df <- best_params_df %>% left_join(best_params_df_est, by = 'subject')

  ### Step 4: generate predictions & update \gamma_v using predictions (for iteration)
  best_params_df$param_pred_A <- best_params_df$param_int + best_params_df$param_slope * best_params_df$param_mean0
  best_params_df$param_pred_B <- best_params_df$param_offset + best_params_df$param_mean0

  # Step 4.a: Compute residuals
  best_params_df$residual_A <- best_params_df$param - best_params_df$param_pred_A
  best_params_df$residual_B <- best_params_df$param - best_params_df$param_pred_B


  # Step 4.b: Estimate variance over subjects (grouping by voxel)
  residual_variance <- best_params_df %>%
    filter(use) %>%
    group_by(voxel) %>%
    # summarize(
    #   n_obs = dplyr::n(),
    #   # Return 0 variance if only one subjects for voxel
    #   variance_A = if(n_obs >= 2) var(residual_A, na.rm = TRUE) else 0,
    #   variance_B = if(n_obs >= 2) var(residual_B, na.rm = TRUE) else 0
    # )
    summarize(variance_A = var(residual_A, na.rm = TRUE), # TODO: Do we need above fix?
              variance_B = var(residual_B, na.rm = TRUE))

  # Create xifti objects for variance datasets for OLS
  variance_xii_A <- create_xifti_with_padding(residual_variance, "variance_A", xii)
  variance_xii_B <- create_xifti_with_padding(residual_variance, "variance_B", xii)

  # Bring variances back to main dataframe
  best_params_df <- best_params_df %>%
    select(-any_of(c("variance_A", "variance_B"))) %>%
    left_join(residual_variance, by = "voxel")

  print("Saving predictions, residuals, variances, variance xiftis for Ordinary Least Squares Method")
  return(list(
    best_params_df = best_params_df,
    best_params_df_est = best_params_df_est,
    variance_xii_A = variance_xii_A,
    variance_xii_B = variance_xii_B
  ))
}

#' Estimate WLS Parameters and Calculate Final Variance
#'
#' Performs weighted least squares estimation using inverse variance weights
#' and calculates final residual variances for Model A and Model B predictions.
#'
#' @param best_params_df Dataframe with parameter estimates, population means, and variances
#' @param truncate Logical, whether to truncate negative slopes
#' @param mean0_xii Template xifti object for creating variance visualizations
#'
#' @return List containing estimates and variance objects
#' @keywords internal
estimate_wls_and_variance <- function(best_params_df, truncate, mean0_xii) {
  # Now that we have variances, we will use WLS method to re-estimate the slope and intercept
  # gamma_iv ->param , gamma_v -> param_mean0
  # weight_A -> inverse square root of variances for model A
  # weight_B -> inverse square root of variances for model B

  print("Running weighted Least squares method")

  best_params_df$weight_A <- 1 / sqrt(best_params_df$variance_A)
  best_params_df$weight_B <- 1 / sqrt(best_params_df$variance_B) # TODO: Do we need below fix?

  # best_params_df$weight_A <- ifelse(best_params_df$variance_A > 1e-10,
  #                                   1 / sqrt(best_params_df$variance_A), 0)
  # best_params_df$weight_B <- ifelse(best_params_df$variance_B > 1e-10,
  #                                   1 / sqrt(best_params_df$variance_B), 0)

  # Now calculate weighted estimates with proper weights
  # best_params_df$weight_A[ !is.finite(best_params_df$weight_A) ] <- 0  # zero-out NA / Inf
  # best_params_df$weight_B[ !is.finite(best_params_df$weight_B) ] <- 0 # FIX FOR ROUND
  best_params_df_est <- best_params_df %>%
    filter(use) %>%
    group_by(subject) %>%
    summarize(param_offset = weighted.mean(param - param_mean0, weight_B, na.rm = TRUE),
              param_slope = slope_fun_WLS(param, param_mean0, weight_A, lm=FALSE),
              param_int = weighted.mean(param, weight_A, na.rm = TRUE) -
                param_slope * weighted.mean(param_mean0, weight_A, na.rm = TRUE))

  if(truncate) best_params_df_est <- truncate_negative_slopes(best_params_df,
                                                              best_params_df_est,
                                                              method = "WLS")


  # Keep subjects in numerical order
  best_params_df_est <- arrange(best_params_df_est, subject) #new param_offset, param_slope, param_int


  best_params_df <- best_params_df %>%
    select(-any_of(c("param_int", "param_slope", "param_offset"))) %>%
    left_join(best_params_df_est, by = "subject")

  # "variance A","variance_B", "param_pred_A","param_pred_B", "residual_A", "residual_B"

  best_params_df$param_pred_A <- best_params_df$param_int + best_params_df$param_slope * best_params_df$param_mean0
  best_params_df$param_pred_B <- best_params_df$param_offset + best_params_df$param_mean0

  # Compute residuals for final variance calculation
  best_params_df$residual_A <- best_params_df$weight_A * (best_params_df$param - best_params_df$param_pred_A)
  best_params_df$residual_B <- best_params_df$weight_B * (best_params_df$param - best_params_df$param_pred_B)

  # Calculate final variance
  residual_variance <- best_params_df %>%
    filter(use) %>%
    group_by(voxel) %>%
    # summarize(
    #   n_obs = dplyr::n(),
    #   variance_A = if(n_obs >= 2) var(residual_A, na.rm = TRUE) else 0,
    #   variance_B = if(n_obs >= 2) var(residual_B, na.rm = TRUE) else 0
    # )
    summarize(variance_A = var(residual_A, na.rm = TRUE), # TODO: Do we need above fix?
              variance_B = var(residual_B, na.rm = TRUE))

  # Create xifti objects for variance datasets
  variance_xii_A <- create_xifti_with_padding(residual_variance, "variance_A", mean0_xii)
  variance_xii_B <- create_xifti_with_padding(residual_variance, "variance_B", mean0_xii)

  # Bring variances back to main dataframe
  best_params_df <- best_params_df %>%
    select(-any_of(c("variance_A", "variance_B"))) %>%
    left_join(residual_variance, by = "voxel")

  print('Saving predictions, residuals, variances, weights, variance xiftis for Weighted Least Squares Method')

  return(list(
    best_params_df = best_params_df,
    best_params_df_est = best_params_df_est,
    variance_xii_A = variance_xii_A,
    variance_xii_B = variance_xii_B
  ))
}



#' Truncate Negative Slopes and Re-estimate Intercepts
#'
#' This function finds subjects with negative slopes in the parameter estimates
#' and sets their slopes to zero. It then re-estimates their intercepts as the
#' simple mean of the parameter values. The resulting estimates are merged back
#' into the original set of estimates.
#'
#' @param best_params_df The full dataframe containing parameter values,
#'   population means, and subject IDs.
#' @param best_params_df_est A summarized dataframe with columns `subject`,
#'   `param_offset`, `param_slope`, and `param_int`.
#'
#' @return A modified version of `best_params_df_est` with negative slopes
#'   truncated to 0 and corresponding intercepts re-estimated.
#'
#' @keywords internal
truncate_negative_slopes <- function(best_params_df, best_params_df_est, method = "OLS") {
  print("In truncataion*****************")
  if (sum(best_params_df_est$param_slope < 0) > 0) {
    if (method == "OLS") {
      cat("OLS: Some subjects had negative slopes, truncating...\n")
      subj_neg_slope <- best_params_df_est$subject[best_params_df_est$param_slope < 0]

      best_params_df_est2 <- best_params_df %>%
        filter(subject %in% subj_neg_slope) %>%
        filter(use) %>%
        group_by(subject) %>%
        summarize(param_offset = mean(param - param_mean0),
                  param_slope = 0,
                  param_int = mean(param))

      best_params_df_est <- rbind(filter(best_params_df_est, param_slope >= 0), best_params_df_est2)
    } else if (method == "WLS") {
      cat("WLS: Some subjects had negative slopes, truncating...\n")
      subj_neg_slope <- best_params_df_est$subject[best_params_df_est$param_slope < 0]

      best_params_df_est2 <- best_params_df %>%
        filter(subject %in% subj_neg_slope) %>%
        filter(use) %>%
        group_by(subject) %>%
        summarize(param_offset = weighted.mean(param - param_mean0, weight_B, na.rm = TRUE),
                  param_slope = 0,
                  param_int = weighted.mean(param, weight_B, na.rm = TRUE))

      best_params_df_est <- rbind(filter(best_params_df_est, param_slope >= 0), best_params_df_est2)
    } else { stop("Attempted to truncate slopes but passed through WLS and OLS methods.")}
  } else {print("***********No need for truncation.")}

  best_params_df_est
}

#' Compute Subject-Level Slope Relative to Population HRF Estimate
#'
#' Computes the slope of subject-level parameter estimates (`gamma_iv`)
#' against population-averaged estimates (`gamma_v`) either via direct formula
#' or linear regression.
#'
#' @param gamma_iv Numeric vector of subject-specific HRF parameter estimates
#'   (length = number of vertices used for a specific subject).
#' @param gamma_v Numeric vector of voxel-wise population-averaged HRF parameters.
#' @param lm Logical; if `TRUE`, uses linear regression. If `FALSE`, computes slope manually.
#'
#' @return Numeric scalar representing the slope estimate.
#'
#' @keywords internal
slope_fun <- function(gamma_iv, gamma_v, lm=TRUE){
  # All arguments must be of length V, where V is the number of vertices
  # included in a particular model (which may differ by subject).
  # gamma_v contains the voxel-wise parameter estimates, averaged across subjects.
  # gamma_iv contains the subject-level, voxel-wise parameter estimates.

  if(!lm) {
    gamma_v_ctr <- gamma_v - mean(gamma_v)
    return(sum(gamma_v_ctr * (gamma_iv - mean(gamma_iv))) / sum(gamma_v_ctr^2))
  }

  if(lm) coefficients(lm(gamma_iv ~ gamma_v))[2]
}

#' Compute Subject-Level Slope Using Weighted Least Squares
#'
#' Computes the weighted slope of subject-level parameter estimates (`gamma_iv`)
#' against population-averaged estimates (`gamma_v`) using inverse variance
#' weighting for improved estimation accuracy.
#'
#' @param gamma_iv Numeric vector of subject-specific HRF parameter estimates
#'   (length = number of vertices used for a specific subject).
#' @param gamma_v Numeric vector of voxel-wise population-averaged HRF parameters.
#' @param weight_A Numeric vector of weights (typically inverse variance) for
#'   weighted least squares estimation.
#' @param lm Logical; currently unused, kept for compatibility. The function
#'   uses manual weighted calculation rather than `lm()`.
#'
#' @return Numeric scalar representing the weighted slope estimate, or NA if
#'   weight normalization fails.
#'
#' @details
#' The function performs weighted least squares slope estimation by:
#' \enumerate{
#'   \item Computing weighted means for centering both variables
#'   \item Normalizing weights to sum to 1 for numerical stability
#'   \item Computing the weighted slope using the manual formula
#' }
#'
#' The manual approach is used instead of `lm()` to maintain control over
#' the weight normalization and reuse normalized weights efficiently.
#'
#' @keywords internal
slope_fun_WLS <- function(gamma_iv, gamma_v, weight_A, lm=FALSE) {
  # weight_A[!is.finite(weight_A)] <- 0
  stopifnot("All weights are zero or invalid - cannot estimate slope" = sum(weight_A) > 0)

  # Center gamma_v using weighted mean
  gamma_v_mean <- weighted.mean(gamma_v, weight_A, na.rm = TRUE)
  gamma_iv_mean <- weighted.mean(gamma_iv, weight_A, na.rm = TRUE)
  gamma_v_ctr <- gamma_v - gamma_v_mean

  # Recompute weights sum to 1
  weight_A <- weight_A/sum(weight_A) # weighted.mean()
  if (abs(sum(weight_A) - 1) > 1e-10) {
    print(paste("ERROR: Sum of weights not normalized to 1; instead:", sum(weight_A)))
    return(NA)
  }

  return(sum(weight_A * gamma_v_ctr * (gamma_iv - gamma_iv_mean)) /
           sum(weight_A * gamma_v_ctr^2))

  # if (lm) {
  #   return(coefficients(lm(gamma_iv ~ gamma_v, weights = weight_A))[2])
  # }
}


#' Extract Parameter Column Names
#'
#' Given a data frame of HRF parameters along with columns \code{"voxel"},
#' \code{"subject"}, and \code{"mask"}, this function returns the names of
#' all remaining numeric columns.
#'
#' @param params_df A data frame containing numeric parameter columns and the
#'   identifier columns \code{"voxel"}, \code{"subject"}, and \code{"mask"}.
#' @return A character vector of column names in \code{params_df} that are numeric
#'   and not one of \code{"voxel"}, \code{"subject"}, or \code{"mask"}.
#' @keywords internal
extract_param_names <- function(params_df) {
  param_names <- names(params_df)[
    vapply(params_df, is.numeric, logical(1)) &
      !(names(params_df) %in% c("voxel", "subject", "mask"))
  ]
  param_names
}


#' Validate Previous Results and Extract Brain Template
#'
#' Validates that \code{workingHRF} and \code{allHRF} results have matching
#' dimensions and extracts a brain template (\code{xii}) for visualization.
#'
#' @param workingHRF_results Results from \code{fit_workingHRF()}.
#' @param allHRF_results Results from \code{fit_allHRFs()}.
#'
#' @return A \code{xii} template object (converted to dscalar format).
#' @keywords internal
validate_previous_results <- function(workingHRF_results, allHRF_results) {

  stopifnot(length(workingHRF_results[["subject_results"]]) == length(allHRF_results[["subject_results"]]))

  # Extract xii objects from both results
  working_xii <- workingHRF_results[["subject_results"]][[1]][["glm_results"]][["bestmodel_xii"]]
  all_xii <- allHRF_results[["subject_results"]][[1]][["glm_result"]][["bestmodel_xii"]]

  # Get dimensions
  working_dims <- dim(as.matrix(working_xii))
  all_dims <- dim(as.matrix(all_xii))

  # Validate dimensions match
  if (!identical(working_dims, all_dims)) {
    stop("Dimension mismatch between workingHRF and allHRF results!\n",
         "  workingHRF dimensions: ", paste(working_dims, collapse=" x "), "\n",
         "  allHRF dimensions: ", paste(all_dims, collapse=" x "))
  }

  cat("Results validation passed - dimensions match:", paste(working_dims, collapse=" x "), "\n")

  # Extract template (following research code pattern)
  xii <- ciftiTools::convert_xifti(working_xii, "dscalar")

  cat("Brain template extracted successfully\n")

  return(xii)
}

#' Create xifti Object from Filtered Data with NA Padding
#'
#' Takes a dataframe with voxel-wise data and creates an xifti object
#' by padding missing voxels with NA to match template dimensions.
#'
#' @param summary_df Dataframe containing 'voxel' column and data column
#' @param data_col Name of the column containing the values to extract
#' @param template_xii Template xifti object to match dimensions and structure
#' @return New xifti object with data from summary_df and NA padding
#' @keywords internal
create_xifti_with_padding <- function(summary_df, data_col, template_xii) {
  # Create full-length vector with NAs
  full_vector <- rep(NA, nrow(as.matrix(template_xii)))

  # Fill in available values
  if (data_col %in% names(summary_df) && "voxel" %in% names(summary_df)) {
    full_vector[summary_df$voxel] <- summary_df[[data_col]]
  } else {
    stop("summary_df must contain 'voxel' column and specified data_col: ", data_col)
  }

  # Create new xifti object
  new_xii <- ciftiTools::newdata_xifti(template_xii, full_vector)

  return(new_xii)
}
