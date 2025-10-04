#' Data Validation Functions for Surplus Production Model
#'
#' This file contains comprehensive validation functions for CPUE and catch data
#' used in the Pella-Tomlinson surplus production model, following the data schema
#' specifications and quality assessment requirements.
#'
#' @name validation
NULL

#' Detect Outliers in Numeric Data
#'
#' Helper function to detect outliers using various statistical methods.
#'
#' @param x Numeric vector to test for outliers
#' @param method Character string specifying the method: "iqr", "zscore"
#' @param threshold Numeric threshold for outlier detection
#'
#' @return Integer vector of outlier indices
#'
#' @keywords internal
detect_outliers <- function(x, method = c("iqr", "zscore"), threshold = NULL) {
  method <- match.arg(method)

  if (length(x) < 4) {
    return(integer(0))
  }

  outlier_indices <- integer(0)

  if (method == "iqr") {
    if (is.null(threshold)) threshold <- 1.5

    q1 <- quantile(x, 0.25, na.rm = TRUE)
    q3 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- q3 - q1

    lower_bound <- q1 - threshold * iqr
    upper_bound <- q3 + threshold * iqr

    outlier_indices <- which(x < lower_bound | x > upper_bound)
  } else if (method == "zscore") {
    if (is.null(threshold)) threshold <- 3.0

    mean_x <- mean(x, na.rm = TRUE)
    sd_x <- sd(x, na.rm = TRUE)

    if (sd_x > 0) {
      z_scores <- abs((x - mean_x) / sd_x)
      outlier_indices <- which(z_scores > threshold)
    }
  }

  return(outlier_indices)
}

#' Validate CPUE Data
#'
#' Comprehensive validation of CPUE (Catch Per Unit Effort) data for Antarctic
#' toothfish stock assessment, including data quality checks, outlier detection,
#' and consistency verification.
#'
#' @param cpue Numeric vector of CPUE values
#' @param years Integer vector of corresponding years (optional)
#' @param effort Numeric vector of effort data (optional, for consistency checks)
#' @param catch Numeric vector of catch data (optional, for consistency checks)
#' @param min_cpue Minimum acceptable CPUE value (default: 0.001)
#' @param max_cpue Maximum acceptable CPUE value (default: 100)
#' @param outlier_method Method for outlier detection: "iqr", "zscore", or "none"
#' @param outlier_threshold Threshold for outlier detection (default: 3.0 for zscore, 1.5 for IQR)
#' @param return_details Logical, whether to return detailed validation results
#'
#' @return If return_details=FALSE (default), returns TRUE if valid or stops with error.
#'   If return_details=TRUE, returns a list with validation results and diagnostics.
#'
#' @details
#' Validation checks performed:
#' \itemize{
#'   \item Basic data type and structure validation
#'   \item Range checks (positive, finite values within acceptable bounds)
#'   \item Missing value detection
#'   \item Outlier detection using specified method
#'   \item Temporal consistency (if years provided)
#'   \item Cross-validation with effort and catch data (if provided)
#'   \item Data quality assessment metrics
#' }
#'
#' CPUE data requirements per data schema:
#' \itemize{
#'   \item Must be positive numeric values
#'   \item No missing values allowed
#'   \item Values should be within reasonable biological bounds
#'   \item Temporal series should show realistic variation
#' }
#'
#' @examples
#' \dontrun{
#' # Basic validation
#' cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
#' validate_cpue_data(cpue)
#'
#' # Detailed validation with outlier detection
#' years <- 2000:2005
#' effort <- c(500, 550, 600, 450, 500, 650)
#' catch <- c(1000, 1100, 1200, 950, 1050, 1300)
#'
#' result <- validate_cpue_data(
#'   cpue = cpue,
#'   years = years,
#'   effort = effort,
#'   catch = catch,
#'   outlier_method = "zscore",
#'   return_details = TRUE
#' )
#' }
#'
#' @export
validate_cpue_data <- function(cpue, years = NULL, effort = NULL, catch = NULL,
                               min_cpue = 0.001, max_cpue = 100,
                               outlier_method = c("iqr", "zscore", "none"),
                               outlier_threshold = NULL,
                               return_details = FALSE) {
  outlier_method <- match.arg(outlier_method)

  # Set default outlier thresholds
  if (is.null(outlier_threshold)) {
    outlier_threshold <- switch(outlier_method,
      "zscore" = 3.0,
      "iqr" = 1.5,
      "none" = NA
    )
  }

  # Initialize validation results
  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    outliers = integer(0),
    quality_metrics = list(),
    summary = list()
  )

  # Basic data type validation
  if (!is.numeric(cpue)) {
    validation_results$errors <- c(
      validation_results$errors,
      "CPUE data must be numeric"
    )
    validation_results$valid <- FALSE
  }

  if (length(cpue) == 0) {
    validation_results$errors <- c(
      validation_results$errors,
      "CPUE data cannot be empty"
    )
    validation_results$valid <- FALSE
  }

  # If basic validation failed, return early
  if (!validation_results$valid) {
    if (return_details) {
      return(validation_results)
    }
    stop(paste(validation_results$errors, collapse = "; "))
  }

  # Missing value check
  if (any(is.na(cpue))) {
    na_indices <- which(is.na(cpue))
    validation_results$errors <- c(
      validation_results$errors,
      paste(
        "CPUE data contains missing values at positions:",
        paste(na_indices, collapse = ", ")
      )
    )
    validation_results$valid <- FALSE
  }

  # Infinite value check
  if (any(!is.finite(cpue[!is.na(cpue)]))) {
    inf_indices <- which(!is.finite(cpue) & !is.na(cpue))
    validation_results$errors <- c(
      validation_results$errors,
      paste(
        "CPUE data contains infinite values at positions:",
        paste(inf_indices, collapse = ", ")
      )
    )
    validation_results$valid <- FALSE
  }

  # Positive value check
  finite_cpue <- cpue[is.finite(cpue)]
  if (length(finite_cpue) > 0 && any(finite_cpue <= 0)) {
    zero_neg_indices <- which(cpue <= 0 & is.finite(cpue))
    validation_results$errors <- c(
      validation_results$errors,
      paste(
        "CPUE data must be positive. Non-positive values at positions:",
        paste(zero_neg_indices, collapse = ", ")
      )
    )
    validation_results$valid <- FALSE
  }

  # Range validation
  if (length(finite_cpue) > 0) {
    if (any(finite_cpue < min_cpue)) {
      below_min <- which(cpue < min_cpue & is.finite(cpue))
      validation_results$warnings <- c(
        validation_results$warnings,
        paste(
          "CPUE values below minimum threshold", min_cpue,
          "at positions:", paste(below_min, collapse = ", ")
        )
      )
    }

    if (any(finite_cpue > max_cpue)) {
      above_max <- which(cpue > max_cpue & is.finite(cpue))
      validation_results$warnings <- c(
        validation_results$warnings,
        paste(
          "CPUE values above maximum threshold", max_cpue,
          "at positions:", paste(above_max, collapse = ", ")
        )
      )
    }
  }

  # Years validation (if provided)
  if (!is.null(years)) {
    if (length(years) != length(cpue)) {
      validation_results$errors <- c(
        validation_results$errors,
        paste(
          "Years length (", length(years),
          ") must match CPUE length (", length(cpue), ")"
        )
      )
      validation_results$valid <- FALSE
    } else {
      # Check for duplicate years
      if (any(duplicated(years))) {
        dup_years <- years[duplicated(years)]
        validation_results$warnings <- c(
          validation_results$warnings,
          paste(
            "Duplicate years found:",
            paste(unique(dup_years), collapse = ", ")
          )
        )
      }

      # Check for temporal ordering
      if (any(diff(years) <= 0)) {
        validation_results$warnings <- c(
          validation_results$warnings,
          "Years are not in strictly increasing order"
        )
      }
    }
  }

  # Cross-validation with effort data
  if (!is.null(effort)) {
    if (length(effort) != length(cpue)) {
      validation_results$errors <- c(
        validation_results$errors,
        paste(
          "Effort length (", length(effort),
          ") must match CPUE length (", length(cpue), ")"
        )
      )
      validation_results$valid <- FALSE
    } else if (!is.null(catch)) {
      # Check CPUE = catch/effort relationship
      calculated_cpue <- catch / effort
      if (any(is.finite(calculated_cpue) & is.finite(cpue))) {
        relative_diff <- abs(cpue - calculated_cpue) / calculated_cpue
        large_diff_indices <- which(relative_diff > 0.01 & is.finite(relative_diff))
        if (length(large_diff_indices) > 0) {
          validation_results$warnings <- c(
            validation_results$warnings,
            paste(
              "CPUE values inconsistent with catch/effort at positions:",
              paste(large_diff_indices, collapse = ", ")
            )
          )
        }
      }
    }
  }

  # Outlier detection
  if (outlier_method != "none" && length(finite_cpue) > 3) {
    outliers <- detect_outliers(finite_cpue, method = outlier_method, threshold = outlier_threshold)
    if (length(outliers) > 0) {
      # Map back to original indices
      finite_indices <- which(is.finite(cpue))
      validation_results$outliers <- finite_indices[outliers]
      validation_results$warnings <- c(
        validation_results$warnings,
        paste(
          "Potential outliers detected at positions:",
          paste(validation_results$outliers, collapse = ", ")
        )
      )
    }
  }

  # Calculate quality metrics
  if (length(finite_cpue) > 0) {
    validation_results$quality_metrics <- list(
      n_observations = length(cpue),
      n_valid = length(finite_cpue),
      n_missing = sum(is.na(cpue)),
      n_infinite = sum(!is.finite(cpue) & !is.na(cpue)),
      n_outliers = length(validation_results$outliers),
      mean_cpue = mean(finite_cpue),
      median_cpue = median(finite_cpue),
      sd_cpue = sd(finite_cpue),
      cv_cpue = sd(finite_cpue) / mean(finite_cpue),
      min_cpue = min(finite_cpue),
      max_cpue = max(finite_cpue),
      range_cpue = max(finite_cpue) - min(finite_cpue)
    )
  }

  # Create summary
  validation_results$summary <- list(
    data_quality = ifelse(validation_results$valid, "PASS", "FAIL"),
    n_errors = length(validation_results$errors),
    n_warnings = length(validation_results$warnings),
    completeness = ifelse(length(finite_cpue) > 0,
      length(finite_cpue) / length(cpue) * 100, 0
    )
  )

  # Return results
  if (return_details) {
    return(validation_results)
  } else {
    if (!validation_results$valid) {
      stop(paste(
        "CPUE data validation failed:",
        paste(validation_results$errors, collapse = "; ")
      ))
    }
    if (length(validation_results$warnings) > 0) {
      warning(paste(
        "CPUE data validation warnings:",
        paste(validation_results$warnings, collapse = "; ")
      ))
    }
    return(TRUE)
  }
}

#' Validate Catch Data
#'
#' Comprehensive validation of catch data for Antarctic toothfish stock assessment,
#' including biological plausibility checks, temporal consistency, and data quality
#' assessment.
#'
#' @param catch Numeric vector of catch values (tonnes)
#' @param years Integer vector of corresponding years (optional)
#' @param cpue Numeric vector of CPUE data (optional, for consistency checks)
#' @param effort Numeric vector of effort data (optional, for consistency checks)
#' @param min_catch Minimum acceptable catch value (default: 0)
#' @param max_catch Maximum acceptable catch value (default: 50000 tonnes)
#' @param outlier_method Method for outlier detection: "iqr", "zscore", or "none"
#' @param outlier_threshold Threshold for outlier detection
#' @param return_details Logical, whether to return detailed validation results
#'
#' @return If return_details=FALSE (default), returns TRUE if valid or stops with error.
#'   If return_details=TRUE, returns a list with validation results and diagnostics.
#'
#' @details
#' Validation checks performed:
#' \itemize{
#'   \item Basic data type and structure validation
#'   \item Range checks (non-negative values within biological bounds)
#'   \item Missing value detection
#'   \item Outlier detection for unrealistic catch values
#'   \item Temporal consistency checks
#'   \item Cross-validation with CPUE and effort data
#'   \item Biological plausibility assessment
#' }
#'
#' Catch data requirements per data schema:
#' \itemize{
#'   \item Must be non-negative numeric values (tonnes)
#'   \item No missing values in the time series
#'   \item Values should be within reasonable biological/operational bounds
#'   \item Should show realistic temporal variation patterns
#' }
#'
#' @examples
#' \dontrun{
#' # Basic validation
#' catch <- c(1000, 1100, 1200, 950, 1050, 1300)
#' validate_catch_data(catch)
#'
#' # Detailed validation with cross-checks
#' years <- 2000:2005
#' cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
#' effort <- c(500, 550, 600, 450, 500, 650)
#'
#' result <- validate_catch_data(
#'   catch = catch,
#'   years = years,
#'   cpue = cpue,
#'   effort = effort,
#'   outlier_method = "iqr",
#'   return_details = TRUE
#' )
#' }
#'
#' @export
validate_catch_data <- function(catch, years = NULL, cpue = NULL, effort = NULL,
                                min_catch = 0, max_catch = 50000,
                                outlier_method = c("iqr", "zscore", "none"),
                                outlier_threshold = NULL,
                                return_details = FALSE) {
  outlier_method <- match.arg(outlier_method)

  # Set default outlier thresholds
  if (is.null(outlier_threshold)) {
    outlier_threshold <- switch(outlier_method,
      "zscore" = 3.0,
      "iqr" = 1.5,
      "none" = NA
    )
  }

  # Initialize validation results
  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    outliers = integer(0),
    quality_metrics = list(),
    summary = list()
  )

  # Basic data type validation
  if (!is.numeric(catch)) {
    validation_results$errors <- c(
      validation_results$errors,
      "Catch data must be numeric"
    )
    validation_results$valid <- FALSE
  }

  if (length(catch) == 0) {
    validation_results$errors <- c(
      validation_results$errors,
      "Catch data cannot be empty"
    )
    validation_results$valid <- FALSE
  }

  # If basic validation failed, return early
  if (!validation_results$valid) {
    if (return_details) {
      return(validation_results)
    }
    stop(paste(validation_results$errors, collapse = "; "))
  }

  # Missing value check
  if (any(is.na(catch))) {
    na_indices <- which(is.na(catch))
    validation_results$errors <- c(
      validation_results$errors,
      paste(
        "Catch data contains missing values at positions:",
        paste(na_indices, collapse = ", ")
      )
    )
    validation_results$valid <- FALSE
  }

  # Infinite value check
  if (any(!is.finite(catch[!is.na(catch)]))) {
    inf_indices <- which(!is.finite(catch) & !is.na(catch))
    validation_results$errors <- c(
      validation_results$errors,
      paste(
        "Catch data contains infinite values at positions:",
        paste(inf_indices, collapse = ", ")
      )
    )
    validation_results$valid <- FALSE
  }

  # Non-negative value check
  finite_catch <- catch[is.finite(catch)]
  if (length(finite_catch) > 0 && any(finite_catch < 0)) {
    neg_indices <- which(catch < 0 & is.finite(catch))
    validation_results$errors <- c(
      validation_results$errors,
      paste(
        "Catch data must be non-negative. Negative values at positions:",
        paste(neg_indices, collapse = ", ")
      )
    )
    validation_results$valid <- FALSE
  }

  # Range validation
  if (length(finite_catch) > 0) {
    if (any(finite_catch < min_catch)) {
      below_min <- which(catch < min_catch & is.finite(catch))
      validation_results$warnings <- c(
        validation_results$warnings,
        paste(
          "Catch values below minimum threshold", min_catch,
          "at positions:", paste(below_min, collapse = ", ")
        )
      )
    }

    if (any(finite_catch > max_catch)) {
      above_max <- which(catch > max_catch & is.finite(catch))
      validation_results$warnings <- c(
        validation_results$warnings,
        paste(
          "Catch values above maximum threshold", max_catch,
          "at positions:", paste(above_max, collapse = ", ")
        )
      )
    }
  }

  # Years validation (if provided)
  if (!is.null(years)) {
    if (length(years) != length(catch)) {
      validation_results$errors <- c(
        validation_results$errors,
        paste(
          "Years length (", length(years),
          ") must match catch length (", length(catch), ")"
        )
      )
      validation_results$valid <- FALSE
    } else {
      # Check for duplicate years
      if (any(duplicated(years))) {
        dup_years <- years[duplicated(years)]
        validation_results$warnings <- c(
          validation_results$warnings,
          paste(
            "Duplicate years found:",
            paste(unique(dup_years), collapse = ", ")
          )
        )
      }

      # Check for temporal gaps
      if (length(years) > 1) {
        year_diffs <- diff(sort(years))
        if (any(year_diffs > 1)) {
          validation_results$warnings <- c(
            validation_results$warnings,
            "Temporal gaps detected in catch time series"
          )
        }
      }
    }
  }

  # Cross-validation with CPUE and effort
  if (!is.null(cpue) && !is.null(effort)) {
    if (length(cpue) != length(catch)) {
      validation_results$errors <- c(
        validation_results$errors,
        paste(
          "CPUE length (", length(cpue),
          ") must match catch length (", length(catch), ")"
        )
      )
      validation_results$valid <- FALSE
    } else if (length(effort) != length(catch)) {
      validation_results$errors <- c(
        validation_results$errors,
        paste(
          "Effort length (", length(effort),
          ") must match catch length (", length(catch), ")"
        )
      )
      validation_results$valid <- FALSE
    } else {
      # Check catch = cpue * effort relationship
      calculated_catch <- cpue * effort
      if (any(is.finite(calculated_catch) & is.finite(catch))) {
        relative_diff <- abs(catch - calculated_catch) / pmax(calculated_catch, 1e-6)
        large_diff_indices <- which(relative_diff > 0.01 & is.finite(relative_diff))
        if (length(large_diff_indices) > 0) {
          validation_results$warnings <- c(
            validation_results$warnings,
            paste(
              "Catch values inconsistent with CPUE*effort at positions:",
              paste(large_diff_indices, collapse = ", ")
            )
          )
        }
      }
    }
  }

  # Outlier detection
  if (outlier_method != "none" && length(finite_catch) > 3) {
    outliers <- detect_outliers(finite_catch, method = outlier_method, threshold = outlier_threshold)
    if (length(outliers) > 0) {
      # Map back to original indices
      finite_indices <- which(is.finite(catch))
      validation_results$outliers <- finite_indices[outliers]
      validation_results$warnings <- c(
        validation_results$warnings,
        paste(
          "Potential outliers detected at positions:",
          paste(validation_results$outliers, collapse = ", ")
        )
      )
    }
  }

  # Biological plausibility checks
  if (length(finite_catch) > 1) {
    # Check for extreme year-to-year variations
    catch_changes <- abs(diff(finite_catch))
    median_catch <- median(finite_catch)
    extreme_changes <- which(catch_changes > 2 * median_catch)
    if (length(extreme_changes) > 0) {
      validation_results$warnings <- c(
        validation_results$warnings,
        paste(
          "Extreme year-to-year catch changes detected between years:",
          paste(extreme_changes, extreme_changes + 1, sep = "-", collapse = ", ")
        )
      )
    }
  }

  # Calculate quality metrics
  if (length(finite_catch) > 0) {
    validation_results$quality_metrics <- list(
      n_observations = length(catch),
      n_valid = length(finite_catch),
      n_missing = sum(is.na(catch)),
      n_infinite = sum(!is.finite(catch) & !is.na(catch)),
      n_zero_catch = sum(finite_catch == 0),
      n_outliers = length(validation_results$outliers),
      total_catch = sum(finite_catch),
      mean_catch = mean(finite_catch),
      median_catch = median(finite_catch),
      sd_catch = sd(finite_catch),
      cv_catch = ifelse(mean(finite_catch) > 0, sd(finite_catch) / mean(finite_catch), NA),
      min_catch = min(finite_catch),
      max_catch = max(finite_catch),
      range_catch = max(finite_catch) - min(finite_catch)
    )
  }

  # Create summary
  validation_results$summary <- list(
    data_quality = ifelse(validation_results$valid, "PASS", "FAIL"),
    n_errors = length(validation_results$errors),
    n_warnings = length(validation_results$warnings),
    completeness = ifelse(length(finite_catch) > 0,
      length(finite_catch) / length(catch) * 100, 0
    )
  )

  # Return results
  if (return_details) {
    return(validation_results)
  } else {
    if (!validation_results$valid) {
      stop(paste(
        "Catch data validation failed:",
        paste(validation_results$errors, collapse = "; ")
      ))
    }
    if (length(validation_results$warnings) > 0) {
      warning(paste(
        "Catch data validation warnings:",
        paste(validation_results$warnings, collapse = "; ")
      ))
    }
    return(TRUE)
  }
}

#' Validate Pella-Tomlinson Model Parameters
#'
#' Comprehensive validation of parameters for the Pella-Tomlinson surplus production
#' model, ensuring biological plausibility and mathematical constraints.
#'
#' @param parameters Named numeric vector of model parameters containing:
#'   r (intrinsic growth rate), K (carrying capacity), m (shape parameter),
#'   q (catchability), sigma_proc (process error), sigma_obs (observation error)
#' @param return_details Logical, whether to return detailed validation results
#'
#' @return If return_details=FALSE (default), returns TRUE if valid or stops with error.
#'   If return_details=TRUE, returns a list with validation results and diagnostics.
#'
#' @details
#' Parameter validation checks:
#' \itemize{
#'   \item r (intrinsic growth rate): 0.001 < r < 2.0 per year
#'   \item K (carrying capacity): 1000 < K < 1,000,000 tonnes
#'   \item m (Pella-Tomlinson shape): 0.1 < m < 10.0
#'   \item q (catchability coefficient): 1e-6 < q < 1.0
#'   \item sigma_proc (process error): 0.001 < sigma_proc < 1.0
#'   \item sigma_obs (observation error): 0.001 < sigma_obs < 2.0
#' }
#'
#' Biological constraints:
#' \itemize{
#'   \item MSY = r*K*m^(m/(m-1))/(m+1)^((m+1)/(m-1)) must be > 0
#'   \item BMSY = K/(m+1)^(1/(m-1)) must be > 0 and < K
#'   \item Parameters must result in stable population dynamics
#' }
#'
#' @examples
#' \dontrun{
#' # Validate default parameters
#' params <- c(
#'   r = 0.1, K = 50000, m = 2.0, q = 0.001,
#'   sigma_proc = 0.1, sigma_obs = 0.2
#' )
#' validate_pt_parameters(params)
#'
#' # Detailed validation
#' result <- validate_pt_parameters(params, return_details = TRUE)
#' print(result$biological_metrics)
#' }
#'
#' @export
validate_pt_parameters <- function(parameters, return_details = FALSE) {
  # Initialize validation results
  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    biological_metrics = list(),
    summary = list()
  )

  # Required parameter names
  required_params <- c("r", "K", "m", "q", "sigma_proc", "sigma_obs")

  # Basic structure validation
  if (!is.numeric(parameters)) {
    validation_results$errors <- c(
      validation_results$errors,
      "Parameters must be numeric"
    )
    validation_results$valid <- FALSE
  }

  if (is.null(names(parameters))) {
    validation_results$errors <- c(
      validation_results$errors,
      "Parameters must have names"
    )
    validation_results$valid <- FALSE
  }

  # Check for required parameters
  if (!is.null(names(parameters))) {
    missing_params <- setdiff(required_params, names(parameters))
    if (length(missing_params) > 0) {
      validation_results$errors <- c(
        validation_results$errors,
        paste(
          "Missing required parameters:",
          paste(missing_params, collapse = ", ")
        )
      )
      validation_results$valid <- FALSE
    }

    extra_params <- setdiff(names(parameters), required_params)
    if (length(extra_params) > 0) {
      validation_results$warnings <- c(
        validation_results$warnings,
        paste(
          "Unknown parameters detected:",
          paste(extra_params, collapse = ", ")
        )
      )
    }
  }

  # If basic validation failed, return early
  if (!validation_results$valid) {
    if (return_details) {
      return(validation_results)
    }
    stop(paste(validation_results$errors, collapse = "; "))
  }

  # Individual parameter validation
  param_bounds <- list(
    r = c(0.001, 2.0),
    K = c(1000, 1000000),
    m = c(0.1, 10.0),
    q = c(1e-6, 1.0),
    sigma_proc = c(0.001, 1.0),
    sigma_obs = c(0.001, 2.0)
  )

  for (param_name in required_params) {
    if (param_name %in% names(parameters)) {
      param_value <- parameters[param_name]
      bounds <- param_bounds[[param_name]]

      # Check for finite values
      if (!is.finite(param_value)) {
        validation_results$errors <- c(
          validation_results$errors,
          paste("Parameter", param_name, "must be finite")
        )
        validation_results$valid <- FALSE
        next
      }

      # Check for positive values
      if (param_value <= 0) {
        validation_results$errors <- c(
          validation_results$errors,
          paste("Parameter", param_name, "must be positive")
        )
        validation_results$valid <- FALSE
        next
      }

      # Check bounds
      if (param_value < bounds[1]) {
        validation_results$warnings <- c(
          validation_results$warnings,
          paste(
            "Parameter", param_name, "is below recommended minimum:",
            bounds[1]
          )
        )
      }

      if (param_value > bounds[2]) {
        validation_results$warnings <- c(
          validation_results$warnings,
          paste(
            "Parameter", param_name, "is above recommended maximum:",
            bounds[2]
          )
        )
      }
    }
  }

  # Biological plausibility checks
  if (all(required_params %in% names(parameters)) && validation_results$valid) {
    r <- parameters["r"]
    K <- parameters["K"]
    m <- parameters["m"]
    q <- parameters["q"]

    tryCatch(
      {
        # Calculate MSY using correct Pella-Tomlinson formula
        if (m > 1.001) { # Avoid numerical issues near m=1
          # MSY = r*K*m^(m/(m-1)) / (m+1)^((m+1)/(m-1))
          term1 <- m^(m / (m - 1))
          term2 <- (m + 1)^((m + 1) / (m - 1))
          msy <- r * K * term1 / term2
        } else {
          # For m ≈ 1, use Fox model approximation: MSY = r*K/e
          msy <- r * K / exp(1)
        }

        # Calculate BMSY using correct formula
        if (m > 0.1) {
          # BMSY = K / (m+1)^(1/(m-1))
          bmsy <- K / (m + 1)^(1 / (m - 1))
        } else {
          # For very small m, use limiting case
          bmsy <- K * exp(-1 / m)
        }

        # Calculate reference points
        validation_results$biological_metrics <- list(
          msy = as.numeric(msy),
          bmsy = as.numeric(bmsy),
          fmsy = as.numeric(msy / bmsy),
          bmsy_k_ratio = as.numeric(bmsy / K),
          productivity = as.numeric(r * K / 4) # Approximate productivity metric
        )

        # Validate reference points
        if (!is.finite(msy) || msy <= 0) {
          validation_results$errors <- c(
            validation_results$errors,
            "Parameters result in invalid MSY"
          )
          validation_results$valid <- FALSE
        }

        if (!is.finite(bmsy) || bmsy <= 0 || bmsy >= K) {
          validation_results$errors <- c(
            validation_results$errors,
            "Parameters result in invalid BMSY"
          )
          validation_results$valid <- FALSE
        }

        # Check for reasonable productivity
        if (is.finite(msy) && msy > K) {
          validation_results$warnings <- c(
            validation_results$warnings,
            "MSY is greater than carrying capacity (K)"
          )
        }

        # Check shape parameter plausibility
        if (m < 0.5) {
          validation_results$warnings <- c(
            validation_results$warnings,
            "Shape parameter m < 0.5 may indicate overly convex production"
          )
        }

        if (m > 5) {
          validation_results$warnings <- c(
            validation_results$warnings,
            "Shape parameter m > 5 may indicate unrealistic production dynamics"
          )
        }
      },
      error = function(e) {
        validation_results$errors <- c(
          validation_results$errors,
          paste("Error calculating biological metrics:", e$message)
        )
        validation_results$valid <- FALSE
      }
    )
  }

  # Create summary
  validation_results$summary <- list(
    parameter_check = ifelse(validation_results$valid, "PASS", "FAIL"),
    n_errors = length(validation_results$errors),
    n_warnings = length(validation_results$warnings),
    biological_plausibility = ifelse(length(validation_results$biological_metrics) > 0, "CALCULATED", "FAILED")
  )

  # Return results
  if (return_details) {
    return(validation_results)
  } else {
    if (!validation_results$valid) {
      stop(paste(
        "Parameter validation failed:",
        paste(validation_results$errors, collapse = "; ")
      ))
    }
    if (length(validation_results$warnings) > 0) {
      warning(paste(
        "Parameter validation warnings:",
        paste(validation_results$warnings, collapse = "; ")
      ))
    }
    return(TRUE)
  }
}
