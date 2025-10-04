#' S4 Class Definitions for Surplus Production Model
#'
#' This file contains the core S4 class definitions for the surplus production
#' model package, including the main ProductionModel class and associated
#' validation functions.
#'
#' @include generics.R
#' @name classes
NULL

#' Production Model S4 Class
#'
#' An S4 class to represent a Pella-Tomlinson surplus production model for
#' Antarctic toothfish stock assessment.
#'
#' @slot parameters A named numeric vector containing model parameters:
#'   \describe{
#'     \item{r}{Intrinsic growth rate (per year)}
#'     \item{K}{Carrying capacity (tonnes)}
#'     \item{m}{Pella-Tomlinson shape parameter (dimensionless)}
#'     \item{q}{Catchability coefficient (per unit effort)}
#'     \item{sigma_proc}{Process error standard deviation}
#'     \item{sigma_obs}{Observation error standard deviation}
#'   }
#' @slot data A list containing model data:
#'   \describe{
#'     \item{years}{Integer vector of years}
#'     \item{catch}{Numeric vector of catch data (tonnes)}
#'     \item{cpue}{Numeric vector of CPUE index data}
#'     \item{effort}{Numeric vector of fishing effort data}
#'   }
#' @slot results A list containing model fitting results (empty until fitted):
#'   \describe{
#'     \item{biomass}{Estimated biomass trajectory}
#'     \item{harvest_rate}{Estimated harvest rates}
#'     \item{msy}{Maximum sustainable yield estimate}
#'     \item{bmsy}{Biomass at MSY}
#'     \item{convergence}{Model convergence status}
#'     \item{likelihood}{Negative log-likelihood value}
#'   }
#' @slot fitted A logical indicating whether the model has been fitted
#' @slot model_type A character string identifying the model type ("pella_tomlinson")
#' @slot creation_date A POSIXct timestamp of when the object was created
#'
#' @examples
#' \dontrun{
#' # Create example data
#' years <- 2000:2020
#' catch <- rnorm(length(years), 1000, 100)
#' cpue <- exp(rnorm(length(years), log(1.5), 0.3))
#' effort <- catch / cpue
#'
#' # Create model object
#' model <- ProductionModel(
#'   years = years,
#'   catch = catch,
#'   cpue = cpue,
#'   effort = effort
#' )
#' }
#'
#' @export
setClass(
  "ProductionModel",
  slots = c(
    parameters = "numeric",
    data = "list",
    results = "list",
    fitted = "logical",
    model_type = "character",
    creation_date = "POSIXct"
  ),
  prototype = list(
    parameters = numeric(0),
    data = list(),
    results = list(),
    fitted = FALSE,
    model_type = "pella_tomlinson",
    creation_date = Sys.time()
  ),
  validity = function(object) {
    errors <- character(0)

    # Validate model_type
    if (!identical(object@model_type, "pella_tomlinson")) {
      errors <- c(errors, "model_type must be 'pella_tomlinson'")
    }

    # Validate fitted status
    if (!is.logical(object@fitted) || length(object@fitted) != 1) {
      errors <- c(errors, "fitted must be a single logical value")
    }

    # Validate creation_date
    if (!inherits(object@creation_date, "POSIXct") || length(object@creation_date) != 1) {
      errors <- c(errors, "creation_date must be a single POSIXct value")
    }

    # Validate parameters if present
    if (length(object@parameters) > 0) {
      required_params <- c("r", "K", "m", "q", "sigma_proc", "sigma_obs")
      param_names <- names(object@parameters)

      if (is.null(param_names)) {
        errors <- c(errors, "parameters must have names")
      } else {
        missing_params <- setdiff(required_params, param_names)
        if (length(missing_params) > 0) {
          errors <- c(errors, paste(
            "Missing required parameters:",
            paste(missing_params, collapse = ", ")
          ))
        }

        extra_params <- setdiff(param_names, required_params)
        if (length(extra_params) > 0) {
          errors <- c(errors, paste(
            "Unknown parameters:",
            paste(extra_params, collapse = ", ")
          ))
        }
      }

      # Validate parameter values
      if ("r" %in% names(object@parameters)) {
        if (!is.finite(object@parameters["r"]) || object@parameters["r"] <= 0) {
          errors <- c(errors, "Parameter 'r' must be positive and finite")
        }
      }

      if ("K" %in% names(object@parameters)) {
        if (!is.finite(object@parameters["K"]) || object@parameters["K"] <= 0) {
          errors <- c(errors, "Parameter 'K' must be positive and finite")
        }
      }

      if ("m" %in% names(object@parameters)) {
        if (!is.finite(object@parameters["m"]) || object@parameters["m"] <= 0) {
          errors <- c(errors, "Parameter 'm' must be positive and finite")
        }
      }

      if ("q" %in% names(object@parameters)) {
        if (!is.finite(object@parameters["q"]) || object@parameters["q"] <= 0) {
          errors <- c(errors, "Parameter 'q' must be positive and finite")
        }
      }

      if ("sigma_proc" %in% names(object@parameters)) {
        if (!is.finite(object@parameters["sigma_proc"]) || object@parameters["sigma_proc"] <= 0) {
          errors <- c(errors, "Parameter 'sigma_proc' must be positive and finite")
        }
      }

      if ("sigma_obs" %in% names(object@parameters)) {
        if (!is.finite(object@parameters["sigma_obs"]) || object@parameters["sigma_obs"] <= 0) {
          errors <- c(errors, "Parameter 'sigma_obs' must be positive and finite")
        }
      }
    }

    # Validate data if present
    if (length(object@data) > 0) {
      required_data <- c("years", "catch", "cpue", "effort")
      data_names <- names(object@data)

      if (is.null(data_names)) {
        errors <- c(errors, "data must have names")
      } else {
        missing_data <- setdiff(required_data, data_names)
        if (length(missing_data) > 0) {
          errors <- c(errors, paste(
            "Missing required data elements:",
            paste(missing_data, collapse = ", ")
          ))
        }
      }

      # Validate data consistency
      if (all(required_data %in% names(object@data))) {
        n_years <- length(object@data$years)

        if (length(object@data$catch) != n_years) {
          errors <- c(errors, "catch data length must match years length")
        }

        if (length(object@data$cpue) != n_years) {
          errors <- c(errors, "cpue data length must match years length")
        }

        if (length(object@data$effort) != n_years) {
          errors <- c(errors, "effort data length must match years length")
        }

        # Validate years
        if (!is.integer(object@data$years) && !all(object@data$years == as.integer(object@data$years))) {
          errors <- c(errors, "years must be integers")
        }

        if (any(diff(object@data$years) <= 0)) {
          errors <- c(errors, "years must be in strictly increasing order")
        }

        # Validate catch data
        if (!is.numeric(object@data$catch) || any(!is.finite(object@data$catch)) || any(object@data$catch < 0)) {
          errors <- c(errors, "catch data must be non-negative finite numbers")
        }

        # Validate CPUE data
        if (!is.numeric(object@data$cpue) || any(!is.finite(object@data$cpue)) || any(object@data$cpue <= 0)) {
          errors <- c(errors, "cpue data must be positive finite numbers")
        }

        # Validate effort data
        if (!is.numeric(object@data$effort) || any(!is.finite(object@data$effort)) || any(object@data$effort <= 0)) {
          errors <- c(errors, "effort data must be positive finite numbers")
        }
      }
    }

    # Validate results if fitted
    if (length(object@fitted) == 1 && object@fitted) {
      if (length(object@results) == 0) {
        errors <- c(errors, "results must be present when fitted = TRUE")
      } else {
        required_results <- c("biomass", "harvest_rate", "msy", "bmsy", "convergence", "likelihood")
        result_names <- names(object@results)

        missing_results <- setdiff(required_results, result_names)
        if (length(missing_results) > 0) {
          errors <- c(errors, paste(
            "Missing required results when fitted:",
            paste(missing_results, collapse = ", ")
          ))
        }
      }
    }

    if (length(errors) == 0) TRUE else errors
  }
)
