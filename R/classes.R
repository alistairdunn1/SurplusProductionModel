#' Validation for ProductionModel Objects
#'
#' This file contains the validation function for ProductionModel S3 objects.
#' A ProductionModel is a named list with class \code{"ProductionModel"}
#' containing the following elements:
#'
#' \describe{
#'   \item{parameters}{A named numeric vector of model parameters (r, K, m, q, sigma_proc, sigma_obs, and optionally B0).
#'     For multi-area/index models, parameters may include q.<area>, q.<area>.<label>, or B0.<area>.}
#'   \item{data}{A list containing model data (years, catch, cpue, effort).}
#'   \item{results}{A list containing model fitting results (empty until fitted).}
#'   \item{fitted}{Logical indicating whether the model has been fitted.}
#'   \item{model_type}{Character string identifying the model type (\code{"pella_tomlinson"}).}
#'   \item{creation_date}{POSIXct timestamp of object creation.}
#' }
#'
#' @name classes
NULL

#' Validate a ProductionModel Object
#'
#' Checks that a ProductionModel list has the correct structure and that all
#' values satisfy the required constraints.
#'
#' @param x A ProductionModel object (list)
#' @return \code{TRUE} invisibly if valid; otherwise throws an error describing
#'   all validation failures.
#' @export
validate_ProductionModel <- function(x) {
  errors <- character(0)

  # ---- top-level structure ----
  required_fields <- c("parameters", "data", "results", "fitted", "model_type", "creation_date")
  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    errors <- c(errors, paste("Missing required fields:", paste(missing_fields, collapse = ", ")))
  }

  # ---- model_type ----
  if ("model_type" %in% names(x) && !identical(x$model_type, "pella_tomlinson")) {
    errors <- c(errors, "model_type must be 'pella_tomlinson'")
  }

  # ---- fitted ----
  if ("fitted" %in% names(x)) {
    if (!is.logical(x$fitted) || length(x$fitted) != 1) {
      errors <- c(errors, "fitted must be a single logical value")
    }
  }

  # ---- creation_date ----
  if ("creation_date" %in% names(x)) {
    if (!inherits(x$creation_date, "POSIXct") || length(x$creation_date) != 1) {
      errors <- c(errors, "creation_date must be a single POSIXct value")
    }
  }

  # ---- parameters ----
  if ("parameters" %in% names(x) && length(x$parameters) > 0) {
    core_required <- c("r", "K", "m", "sigma_obs")
    optional_params <- c("q", "sigma_proc", "sigma_process", "B0")
    param_names <- names(x$parameters)

    if (is.null(param_names)) {
      errors <- c(errors, "parameters must have names")
    } else {
      missing_core <- setdiff(core_required, param_names)
      if (length(missing_core) > 0) {
        errors <- c(errors, paste("Missing required parameters:", paste(missing_core, collapse = ", ")))
      }

      is_q_like <- grepl("^q(\\.|$)", param_names)
      is_b0_like <- grepl("^B0(\\.|$)", param_names)

      if (!any(is_q_like)) {
        errors <- c(errors, "At least one catchability parameter 'q' or 'q.<area>[.<label>]' must be provided")
      }

      allowed_core <- c(core_required, optional_params)
      extra_params <- param_names[!(param_names %in% allowed_core | is_q_like | is_b0_like)]
      if (length(extra_params) > 0) {
        errors <- c(errors, paste("Unknown parameters:", paste(extra_params, collapse = ", ")))
      }
    }

    # Individual parameter constraints
    if ("r" %in% names(x$parameters)) {
      if (!is.finite(x$parameters["r"]) || x$parameters["r"] <= 0) {
        errors <- c(errors, "Parameter 'r' must be positive and finite")
      }
    }
    if ("K" %in% names(x$parameters)) {
      if (!is.finite(x$parameters["K"]) || x$parameters["K"] <= 0) {
        errors <- c(errors, "Parameter 'K' must be positive and finite")
      }
    }
    if ("m" %in% names(x$parameters)) {
      if (!is.finite(x$parameters["m"]) || x$parameters["m"] <= 0) {
        errors <- c(errors, "Parameter 'm' must be positive and finite")
      }
    }

    q_like_idx <- grepl("^q(\\.|$)", names(x$parameters))
    if (any(q_like_idx)) {
      q_vals <- x$parameters[q_like_idx]
      if (any(!is.finite(q_vals) | q_vals <= 0)) {
        errors <- c(errors, "All catchability parameters 'q*' must be positive and finite")
      }
    }

    if ("sigma_proc" %in% names(x$parameters)) {
      if (!is.finite(x$parameters["sigma_proc"]) || x$parameters["sigma_proc"] <= 0) {
        errors <- c(errors, "Parameter 'sigma_proc' must be positive and finite")
      }
    }
    if ("sigma_process" %in% names(x$parameters)) {
      if (!is.finite(x$parameters["sigma_process"]) || x$parameters["sigma_process"] <= 0) {
        errors <- c(errors, "Parameter 'sigma_process' must be positive and finite")
      }
    }
    if ("sigma_obs" %in% names(x$parameters)) {
      if (!is.finite(x$parameters["sigma_obs"]) || x$parameters["sigma_obs"] <= 0) {
        errors <- c(errors, "Parameter 'sigma_obs' must be positive and finite")
      }
    }

    b0_like_idx <- grepl("^B0(\\.|$)", names(x$parameters))
    if (any(b0_like_idx)) {
      b0_vals <- x$parameters[b0_like_idx]
      if (any(!is.finite(b0_vals) | b0_vals <= 0)) {
        errors <- c(errors, "All initial biomass parameters 'B0*' must be positive and finite")
      }
    }
  }

  # ---- data ----
  if ("data" %in% names(x) && length(x$data) > 0) {
    required_data <- c("years", "catch", "cpue", "effort")
    data_names <- names(x$data)

    if (is.null(data_names)) {
      errors <- c(errors, "data must have names")
    } else {
      missing_data <- setdiff(required_data, data_names)
      if (length(missing_data) > 0) {
        errors <- c(errors, paste("Missing required data elements:", paste(missing_data, collapse = ", ")))
      }
    }

    if (all(required_data %in% names(x$data))) {
      n_years <- length(x$data$years)

      if (!is.integer(x$data$years) && !all(x$data$years == as.integer(x$data$years))) {
        errors <- c(errors, "years must be integers")
      }
      if (any(diff(x$data$years) <= 0)) {
        errors <- c(errors, "years must be in strictly increasing order")
      }

      first_dim_len <- function(v) {
        if (is.null(dim(v))) length(v) else dim(v)[1]
      }

      if (first_dim_len(x$data$catch) != n_years) {
        errors <- c(errors, "catch must be a vector of length years or have first dimension equal to years")
      }
      if (!is.numeric(x$data$catch) || any(!is.finite(c(x$data$catch))) || any(c(x$data$catch) < 0, na.rm = TRUE)) {
        errors <- c(errors, "catch data must be non-negative finite numbers")
      }

      if (first_dim_len(x$data$cpue) != n_years) {
        errors <- c(errors, "cpue must be a vector/matrix/array with first dimension equal to years")
      }
      if (!is.numeric(x$data$cpue)) {
        errors <- c(errors, "cpue must be numeric")
      } else {
        cpue_vals <- c(x$data$cpue)
        if (any(!is.na(cpue_vals) & (!is.finite(cpue_vals) | cpue_vals <= 0))) {
          errors <- c(errors, "cpue data must be positive finite numbers")
        }
      }

      if (first_dim_len(x$data$effort) != n_years) {
        errors <- c(errors, "effort must be a vector/matrix with first dimension equal to years")
      }
      if (!is.numeric(x$data$effort)) {
        errors <- c(errors, "effort must be numeric")
      } else {
        eff_vals <- c(x$data$effort)
        if (any(!is.na(eff_vals) & (!is.finite(eff_vals) | eff_vals <= 0))) {
          errors <- c(errors, "effort values must be positive and finite where present")
        }
      }
    }
  }

  # ---- results when fitted ----
  if ("fitted" %in% names(x) && isTRUE(x$fitted)) {
    if (length(x$results) == 0) {
      errors <- c(errors, "results must be present when fitted = TRUE")
    } else {
      required_results <- c("biomass", "harvest_rate", "msy", "bmsy", "convergence", "likelihood")
      missing_results <- setdiff(required_results, names(x$results))
      if (length(missing_results) > 0) {
        errors <- c(errors, paste("Missing required results when fitted:", paste(missing_results, collapse = ", ")))
      }
    }
  }

  if (length(errors) > 0) stop(paste(errors, collapse = "\n"))
  invisible(TRUE)
}
