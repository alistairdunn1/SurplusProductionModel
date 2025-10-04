#' Constructor Functions for Production Model Objects
#'
#' This file contains constructor functions and methods for creating and
#' manipulating ProductionModel objects.
#'
#' @include generics.R classes.R
#' @name constructors
NULL

#' Create a Production Model Object
#'
#' Constructor function for creating a new ProductionModel object with
#' validated input data.
#'
#' @param years Integer vector of years for the assessment
#' @param catch Numeric vector of catch data (tonnes) corresponding to years
#' @param cpue Numeric vector of CPUE index data corresponding to years
#' @param effort Numeric vector of fishing effort data corresponding to years
#' @param parameters Optional named numeric vector of model parameters. If not
#'   provided, default initial values will be used.
#'
#' @return A ProductionModel object
#'
#' @details
#' The function validates all input data according to the data schema:
#' \itemize{
#'   \item Years must be integers in strictly increasing order
#'   \item Catch data must be non-negative finite numbers
#'   \item CPUE data must be positive finite numbers
#'   \item Effort data must be positive finite numbers
#'   \item All data vectors must have the same length
#' }
#'
#' Default parameter values are based on typical Antarctic toothfish
#' stock characteristics:
#' \itemize{
#'   \item r = 0.1 (intrinsic growth rate per year)
#'   \item K = 50000 (carrying capacity in tonnes)
#'   \item m = 2.0 (Pella-Tomlinson shape parameter, Schaefer model)
#'   \item q = 0.001 (catchability coefficient)
#'   \item sigma_proc = 0.1 (process error standard deviation)
#'   \item sigma_obs = 0.2 (observation error standard deviation)
#' }
#'
#' @examples
#' \dontrun{
#' # Create example data
#' years <- 2000:2020
#' catch <- rnorm(length(years), 1000, 100)
#' cpue <- exp(rnorm(length(years), log(1.5), 0.3))
#' effort <- catch / cpue
#'
#' # Create model with default parameters
#' model <- ProductionModel(
#'   years = years,
#'   catch = catch,
#'   cpue = cpue,
#'   effort = effort
#' )
#'
#' # Create model with custom parameters
#' custom_params <- c(
#'   r = 0.15,
#'   K = 60000,
#'   m = 1.5,
#'   q = 0.0015,
#'   sigma_proc = 0.12,
#'   sigma_obs = 0.18
#' )
#'
#' model_custom <- ProductionModel(
#'   years = years,
#'   catch = catch,
#'   cpue = cpue,
#'   effort = effort,
#'   parameters = custom_params
#' )
#' }
#'
#' @export
ProductionModel <- function(years, catch, cpue, effort, parameters = NULL) {
  # Validate inputs using checkmate
  checkmate::assert_integerish(years, any.missing = FALSE, min.len = 1)
  checkmate::assert_numeric(catch, any.missing = FALSE, lower = 0, finite = TRUE)
  checkmate::assert_numeric(cpue, any.missing = FALSE, lower = 0, finite = TRUE)
  checkmate::assert_numeric(effort, any.missing = FALSE, lower = 0, finite = TRUE)

  # Check that all data vectors have same length
  n_years <- length(years)
  if (length(catch) != n_years) {
    stop("catch data length (", length(catch), ") must match years length (", n_years, ")")
  }
  if (length(cpue) != n_years) {
    stop("cpue data length (", length(cpue), ") must match years length (", n_years, ")")
  }
  if (length(effort) != n_years) {
    stop("effort data length (", length(effort), ") must match years length (", n_years, ")")
  }

  # Check that years are in strictly increasing order
  if (any(diff(years) <= 0)) {
    stop("years must be in strictly increasing order")
  }

  # Convert years to integer if needed
  years <- as.integer(years)

  # Set default parameters if not provided
  if (is.null(parameters)) {
    parameters <- c(
      r = 0.1, # Intrinsic growth rate per year
      K = 50000, # Carrying capacity in tonnes
      m = 2.0, # Pella-Tomlinson shape (2 = Schaefer model)
      q = 0.001, # Catchability coefficient
      sigma_proc = 0.1, # Process error standard deviation
      sigma_obs = 0.2 # Observation error standard deviation
    )
  } else {
    # Validate provided parameters
    checkmate::assert_numeric(parameters, any.missing = FALSE, finite = TRUE)
    checkmate::assert_names(names(parameters),
      must.include = c("r", "K", "m", "q", "sigma_proc", "sigma_obs")
    )

    # Check parameter constraints
    if (parameters["r"] <= 0) stop("Parameter 'r' must be positive")
    if (parameters["K"] <= 0) stop("Parameter 'K' must be positive")
    if (parameters["m"] <= 0) stop("Parameter 'm' must be positive")
    if (parameters["q"] <= 0) stop("Parameter 'q' must be positive")
    if (parameters["sigma_proc"] <= 0) stop("Parameter 'sigma_proc' must be positive")
    if (parameters["sigma_obs"] <= 0) stop("Parameter 'sigma_obs' must be positive")
  }

  # Create data list
  data_list <- list(
    years = years,
    catch = catch,
    cpue = cpue,
    effort = effort
  )

  # Create the object
  new("ProductionModel",
    parameters = parameters,
    data = data_list,
    results = list(),
    fitted = FALSE,
    model_type = "pella_tomlinson",
    creation_date = Sys.time()
  )
}

#' Get Model Parameters
#'
#' Extract parameter values from a ProductionModel object.
#'
#' @param object A ProductionModel object
#' @return Named numeric vector of parameters
#' @export
setMethod("parameters", "ProductionModel", function(object) {
  object@parameters
})

#' Get Model Data
#'
#' Extract data from a ProductionModel object.
#'
#' @param object A ProductionModel object
#' @return List containing model data (years, catch, cpue, effort)
#' @export
setMethod("model_data", "ProductionModel", function(object) {
  object@data
})

#' Get Model Results
#'
#' Extract fitting results from a ProductionModel object.
#'
#' @param object A ProductionModel object
#' @return List containing model results (empty if not fitted)
#' @export
setMethod("results", "ProductionModel", function(object) {
  object@results
})

#' Check if Model is Fitted
#'
#' Check whether a ProductionModel object has been fitted.
#'
#' @param object A ProductionModel object
#' @return Logical indicating if model has been fitted
#' @export
setMethod("fitted", "ProductionModel", function(object) {
  object@fitted
})

#' Print Method for ProductionModel
#'
#' Print summary information about a ProductionModel object.
#'
#' @param x A ProductionModel object
#' @param ... Additional arguments (currently unused)
#' @return Invisibly returns the object
#' @export
setMethod("print", "ProductionModel", function(x, ...) {
  cat("Production Model Object\n")
  cat("======================\n\n")

  cat("Model Type:", x@model_type, "\n")
  cat("Created:", format(x@creation_date, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Fitted:", ifelse(x@fitted, "Yes", "No"), "\n\n")

  if (length(x@data) > 0) {
    cat("Data Summary:\n")
    cat(
      "  Years:", range(x@data$years)[1], "-", range(x@data$years)[2],
      "(", length(x@data$years), "observations)\n"
    )
    cat("  Catch range:", sprintf(
      "%.1f - %.1f tonnes",
      min(x@data$catch), max(x@data$catch)
    ), "\n")
    cat("  CPUE range:", sprintf(
      "%.3f - %.3f",
      min(x@data$cpue), max(x@data$cpue)
    ), "\n")
    cat("  Effort range:", sprintf(
      "%.1f - %.1f",
      min(x@data$effort), max(x@data$effort)
    ), "\n\n")
  }

  if (length(x@parameters) > 0) {
    cat("Parameters:\n")
    for (i in seq_along(x@parameters)) {
      cat("  ", names(x@parameters)[i], ":", sprintf("%.4f", x@parameters[i]), "\n")
    }
    cat("\n")
  }

  if (x@fitted && length(x@results) > 0) {
    cat("Results Summary:\n")
    if ("msy" %in% names(x@results)) {
      cat("  MSY:", sprintf("%.1f tonnes", x@results$msy), "\n")
    }
    if ("bmsy" %in% names(x@results)) {
      cat("  BMSY:", sprintf("%.1f tonnes", x@results$bmsy), "\n")
    }
    if ("likelihood" %in% names(x@results)) {
      cat("  Negative log-likelihood:", sprintf("%.2f", x@results$likelihood), "\n")
    }
    if ("convergence" %in% names(x@results)) {
      cat("  Convergence:", ifelse(x@results$convergence == 0, "Success", "Failed"), "\n")
    }
  }

  invisible(x)
})

#' Show Method for ProductionModel
#'
#' Show summary information about a ProductionModel object.
#'
#' @param object A ProductionModel object
#' @return Invisibly returns the object
#' @export
setMethod("show", "ProductionModel", function(object) {
  print(object)
})
