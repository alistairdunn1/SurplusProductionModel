#' Reference Points Calculation for Surplus Production Model
#'
#' This file contains functions for calculating biological reference points
#' from fitted Pella-Tomlinson surplus production models.
#'
#' @name reference-points
NULL

#' Calculate Reference Points from Fitted Model
#'
#' Calculates biological reference points (MSY, BMSY, FMSY) from a fitted
#' Pella-Tomlinson surplus production model.
#'
#' @param model_fit A fitted ProductionModel object (with fitted = TRUE)
#'
#' @return Named list containing reference points:
#'   \describe{
#'     \item{msy}{Maximum sustainable yield (tonnes)}
#'     \item{bmsy}{Biomass at MSY (tonnes)}
#'     \item{fmsy}{Fishing mortality at MSY (per year)}
#'     \item{current_status}{List with current biomass and status relative to reference points}
#'   }
#'
#' @details
#' Reference points are calculated using the Pella-Tomlinson formulation:
#'
#' MSY = r * K * (m-1)^((m-1)/m) / m
#'
#' BMSY = K * (m-1)^(1/m) / m
#'
#' FMSY = MSY / BMSY
#'
#' For special cases:
#' - Schaefer model (m=2): MSY = r*K/4, BMSY = K/2
#' - Fox model (m=1): Uses limiting case formulation
#'
#' Current status is calculated relative to reference points using the
#' final year biomass estimate from the fitted model.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a fitted model
#' ref_points <- calculate_reference_points(fitted_model)
#' print(ref_points)
#'
#' # Extract specific reference points
#' msy <- ref_points$msy
#' bmsy <- ref_points$bmsy
#' }
#'
#' @seealso \code{\link{fit_pella_tomlinson_model}} for model fitting
#'
#' @export
calculate_reference_points <- function(model_fit) {
  # Input validation
  if (!inherits(model_fit, "ProductionModel")) {
    stop("model_fit must be a ProductionModel object")
  }

  if (!model_fit@fitted) {
    stop("Model must be fitted before calculating reference points. Use fit_pella_tomlinson_model() first.")
  }

  if (length(model_fit@parameters) == 0) {
    stop("No fitted parameters found in model object")
  }

  # Extract fitted parameters
  parameters <- model_fit@parameters

  required_params <- c("r", "K", "m")
  if (!all(required_params %in% names(parameters))) {
    stop("Missing required parameters: ", paste(setdiff(required_params, names(parameters)), collapse = ", "))
  }

  r <- parameters[["r"]]
  K <- parameters[["K"]]
  m <- parameters[["m"]]

  # Validate parameter values
  if (any(c(r, K, m) <= 0)) {
    stop("All parameters must be positive")
  }

  # Calculate reference points based on shape parameter
  if (abs(m - 1) < 1e-6) {
    # Fox model (m ≈ 1) - use limiting case
    msy <- r * K / exp(1) # r * K / e
    bmsy <- K / exp(1) # K / e
  } else if (abs(m - 2) < 1e-6) {
    # Schaefer model (m ≈ 2) - analytical solution
    msy <- r * K / 4
    bmsy <- K / 2
  } else {
    # General Pella-Tomlinson case
    if (m <= 1) {
      # Handle edge case for very low m values
      warning("Shape parameter m <= 1 may give unrealistic reference points")
      msy <- r * K / 4 # Use Schaefer approximation
      bmsy <- K / 2
    } else {
      # Standard Pella-Tomlinson formulation
      msy <- r * K * (m - 1)^((m - 1) / m) / m
      bmsy <- K * (m - 1)^(1 / m) / m
    }
  }

  # Calculate FMSY
  if (bmsy > 0) {
    fmsy <- msy / bmsy
  } else {
    fmsy <- NA
    warning("BMSY is zero or negative, cannot calculate FMSY")
  }

  # Calculate current status if biomass results are available
  current_status <- NULL
  if ("biomass" %in% names(model_fit@results) && length(model_fit@results$biomass) > 0) {
    current_biomass <- tail(model_fit@results$biomass, 1) # Final year biomass

    # Calculate status ratios
    b_bmsy_ratio <- current_biomass / bmsy

    # Determine stock status
    if (b_bmsy_ratio > 1.0) {
      status <- "Above BMSY"
    } else if (b_bmsy_ratio > 0.5) {
      status <- "Below BMSY but above half BMSY"
    } else {
      status <- "Below half BMSY (overfished)"
    }

    current_status <- list(
      current_biomass = current_biomass,
      b_bmsy_ratio = b_bmsy_ratio,
      status = status
    )

    # Calculate current harvest rate if available
    if ("harvest_rate" %in% names(model_fit@results)) {
      current_harvest_rate <- tail(model_fit@results$harvest_rate, 1)
      f_fmsy_ratio <- current_harvest_rate / fmsy

      # Add harvest rate status
      if (f_fmsy_ratio > 1.0) {
        harvest_status <- "Overfishing occurring"
      } else {
        harvest_status <- "No overfishing"
      }

      current_status$current_harvest_rate <- current_harvest_rate
      current_status$f_fmsy_ratio <- f_fmsy_ratio
      current_status$harvest_status <- harvest_status
    }
  }

  # Package results
  reference_points <- list(
    msy = msy,
    bmsy = bmsy,
    fmsy = fmsy,
    current_status = current_status,
    model_type = "Pella-Tomlinson",
    shape_parameter = m,
    calculation_date = Sys.time()
  )

  # Add special case identification
  if (abs(m - 1) < 1e-6) {
    reference_points$special_case <- "Fox model (m ≈ 1)"
  } else if (abs(m - 2) < 1e-6) {
    reference_points$special_case <- "Schaefer model (m ≈ 2)"
  } else {
    reference_points$special_case <- "General Pella-Tomlinson model"
  }

  class(reference_points) <- c("pt_reference_points", "list")

  return(reference_points)
}

#' Print Method for Reference Points
#'
#' Print method for pt_reference_points objects.
#'
#' @param x A pt_reference_points object
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the object
#'
#' @export
print.pt_reference_points <- function(x, ...) {
  cat("Pella-Tomlinson Reference Points\n")
  cat("================================\n\n")

  cat("Model Type:", x$special_case, "\n")
  cat("Shape Parameter (m):", sprintf("%.3f", x$shape_parameter), "\n\n")

  cat("Reference Points:\n")
  cat("  MSY:", sprintf("%.1f tonnes", x$msy), "\n")
  cat("  BMSY:", sprintf("%.1f tonnes", x$bmsy), "\n")
  cat("  FMSY:", sprintf("%.3f per year", x$fmsy), "\n\n")

  if (!is.null(x$current_status)) {
    cat("Current Stock Status:\n")
    cat("  Current Biomass:", sprintf("%.1f tonnes", x$current_status$current_biomass), "\n")
    cat("  B/BMSY ratio:", sprintf("%.3f", x$current_status$b_bmsy_ratio), "\n")
    cat("  Status:", x$current_status$status, "\n")

    if (!is.null(x$current_status$current_harvest_rate)) {
      cat("  Current Harvest Rate:", sprintf("%.3f per year", x$current_status$current_harvest_rate), "\n")
      cat("  F/FMSY ratio:", sprintf("%.3f", x$current_status$f_fmsy_ratio), "\n")
      cat("  Harvest Status:", x$current_status$harvest_status, "\n")
    }
  }

  cat("\nCalculated:", format(x$calculation_date, "%Y-%m-%d %H:%M:%S"), "\n")

  invisible(x)
}

#' Biomass Estimation from Fitted Model
#'
#' Extract biomass estimates from a fitted ProductionModel object.
#'
#' @param model_fit A fitted ProductionModel object
#' @param years Optional vector of years to extract (defaults to all fitted years)
#'
#' @return Matrix with columns: year, biomass, harvest_rate (if available)
#'   or vector if single area
#'
#' @details
#' Extracts the estimated biomass trajectory from a fitted surplus production
#' model. If harvest rates were calculated during fitting, they are also included.
#'
#' @examples
#' \dontrun{
#' # Extract all biomass estimates
#' biomass_estimates <- estimate_biomass(fitted_model)
#'
#' # Extract specific years
#' recent_biomass <- estimate_biomass(fitted_model, years = 2015:2020)
#' }
#'
#' @export
estimate_biomass <- function(model_fit, years = NULL) {
  # Input validation
  if (!inherits(model_fit, "ProductionModel")) {
    stop("model_fit must be a ProductionModel object")
  }

  if (!model_fit@fitted) {
    stop("Model must be fitted before extracting biomass estimates")
  }

  if (!"biomass" %in% names(model_fit@results)) {
    stop("No biomass estimates found in fitted model")
  }

  # Extract data
  model_years <- model_fit@data$years
  biomass_values <- model_fit@results$biomass

  # Create results matrix
  results <- data.frame(
    year = model_years,
    biomass = biomass_values
  )

  # Add harvest rates if available
  if ("harvest_rate" %in% names(model_fit@results)) {
    results$harvest_rate <- model_fit@results$harvest_rate
  }

  # Filter to requested years if specified
  if (!is.null(years)) {
    years_available <- intersect(years, model_years)
    if (length(years_available) == 0) {
      stop("None of the requested years are available in the fitted model")
    }
    if (length(years_available) < length(years)) {
      warning("Some requested years are not available in the fitted model")
    }
    results <- results[results$year %in% years_available, ]
  }

  return(results)
}
