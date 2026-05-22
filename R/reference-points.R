#' Reference Points Calculation for Surplus Production Model
#'
#' This file contains functions for calculating biological reference points
#' from fitted Pella-Tomlinson surplus production models.
#'
#' @name reference-points
NULL

#' Set Package-Level Reference Point Defaults
#'
#' Sets package-wide defaults for depletion target reporting. These defaults
#' are used by functions such as \code{calculate_reference_points()},
#' \code{profile_likelihood()}, \code{bayesian_fit()},
#' \code{summary.ProductionModel()}, and \code{print.ProductionModel()}
#' when \code{biomass_target} and/or \code{baseline} are not supplied.
#'
#' @param biomass_target Optional numeric vector of default biomass fractions
#'   to use when a function call omits \code{biomass_target}. Pass
#'   \code{NULL} explicitly to clear the default.
#' @param baseline Optional default biomass baseline, one of
#'   \code{"auto"}, \code{"B0"}, or \code{"K"}.
#'
#' @return Invisibly returns the updated defaults as a named list.
#' @export
set_reference_point_defaults <- function(biomass_target = NULL,
                                         baseline = c("auto", "B0", "K")) {
  if (!missing(biomass_target)) {
    validate_biomass_target(biomass_target)
    do.call(options, stats::setNames(
      list(biomass_target),
      "SurplusProductionModel.biomass_target_default"
    ))
  }
  if (!missing(baseline)) {
    baseline <- match.arg(baseline)
    do.call(options, stats::setNames(
      list(baseline),
      "SurplusProductionModel.baseline_default"
    ))
  }

  invisible(get_reference_point_defaults())
}

#' Get Package-Level Reference Point Defaults
#'
#' Returns current package-level defaults for depletion target reporting.
#'
#' @return Named list with entries \code{biomass_target} and \code{baseline}.
#' @export
get_reference_point_defaults <- function() {
  list(
    biomass_target = getOption("SurplusProductionModel.biomass_target_default", NULL),
    baseline = getOption("SurplusProductionModel.baseline_default", "auto")
  )
}

#' Calculate Reference Points from Fitted Model
#'
#' Calculates biological reference points (MSY, BMSY, FMSY) from a fitted
#' Pella-Tomlinson surplus production model.
#'
#' @param model_fit A fitted ProductionModel object (with fitted = TRUE)
#' @param biomass_target Optional numeric vector of target biomass fractions.
#'   For example, `0.4` requests `B_40%B0` and `F_40%B0` style reference
#'   points. If omitted, uses package default from
#'   \code{get_reference_point_defaults()}.
#' @param baseline Character string indicating which biomass baseline to use
#'   for `biomass_target`: `"auto"` (default), `"B0"`, or `"K"`. In
#'   `"auto"` mode, the function uses fitted `B0` if available and otherwise
#'   falls back to `K`. If omitted, uses package default from
#'   \code{get_reference_point_defaults()}.
#'
#' @return Named list containing reference points:
#'   \describe{
#'     \item{msy}{Maximum sustainable yield (tonnes)}
#'     \item{bmsy}{Biomass at MSY (tonnes)}
#'     \item{fmsy}{Fishing mortality at MSY (per year)}
#'     \item{target_reference_points}{Optional data frame of user-defined
#'       biomass and fishing mortality targets for requested biomass fractions.}
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
#' For user-defined biomass fractions `x`, the function also calculates
#' analytical depletion-based targets:
#'
#' B_x = x * B_base
#'
#' F_x = r / m * (1 - (B_x / K)^(m - 1))
#'
#' where `B_base` is either fitted `B0` or `K`, depending on `baseline`.
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
#'
#' # Calculate B_40%B0 and F_40%B0 style targets
#' target_ref_points <- calculate_reference_points(
#'   fitted_model,
#'   biomass_target = 0.4,
#'   baseline = "B0"
#' )
#' }
#'
#' @seealso \code{\link{fit_pella_tomlinson_model}} for model fitting
#'
#' @export
calculate_reference_points <- function(model_fit, biomass_target = NULL, baseline = c("auto", "B0", "K")) {
  defaults <- resolve_reference_point_defaults(
    biomass_target = biomass_target,
    baseline = baseline,
    biomass_target_missing = missing(biomass_target),
    baseline_missing = missing(baseline)
  )
  biomass_target <- defaults$biomass_target
  baseline <- defaults$baseline

  # Input validation
  if (!inherits(model_fit, "ProductionModel")) {
    stop("model_fit must be a ProductionModel object")
  }

  if (!model_fit$fitted) {
    stop("Model must be fitted before calculating reference points. Use fit_pella_tomlinson_model() first.")
  }

  if (length(model_fit$parameters) == 0) {
    stop("No fitted parameters found in model object")
  }

  # Extract fitted parameters
  parameters <- model_fit$parameters

  required_params <- c("r", "K", "m")
  if (!all(required_params %in% names(parameters))) {
    stop("Missing required parameters: ", paste(setdiff(required_params, names(parameters)), collapse = ", "))
  }

  ref_core <- calculate_reference_points_from_parameters(
    parameters = parameters,
    biomass_target = biomass_target,
    baseline = baseline
  )

  # Calculate current status if biomass results are available
  current_status <- NULL
  if ("biomass" %in% names(model_fit$results) && length(model_fit$results$biomass) > 0) {
    current_biomass <- tail(model_fit$results$biomass, 1) # Final year biomass

    # Calculate status ratios
    b_bmsy_ratio <- current_biomass / ref_core$bmsy

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
    if ("harvest_rate" %in% names(model_fit$results)) {
      current_harvest_rate <- tail(model_fit$results$harvest_rate, 1)
      f_fmsy_ratio <- current_harvest_rate / ref_core$fmsy

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
    msy = ref_core$msy,
    bmsy = ref_core$bmsy,
    fmsy = ref_core$fmsy,
    target_reference_points = ref_core$target_reference_points,
    current_status = current_status,
    model_type = "Pella-Tomlinson",
    shape_parameter = parameters[["m"]],
    calculation_date = Sys.time()
  )

  # Add special case identification
  reference_points$special_case <- ref_core$special_case

  class(reference_points) <- c("pt_reference_points", "list")

  return(reference_points)
}

calculate_reference_points_from_parameters <- function(parameters,
                                                       biomass_target = NULL,
                                                       baseline = c("auto", "B0", "K"),
                                                       warn_on_invalid_m = TRUE) {
  required_params <- c("r", "K", "m")
  if (!all(required_params %in% names(parameters))) {
    stop("Missing required parameters: ", paste(setdiff(required_params, names(parameters)), collapse = ", "))
  }

  baseline <- match.arg(baseline)
  validate_biomass_target(biomass_target)

  r <- parameters[["r"]]
  K <- parameters[["K"]]
  m <- parameters[["m"]]

  if (any(c(r, K, m) <= 0)) {
    stop("All parameters must be positive")
  }

  if (abs(m - 1) < 1e-6) {
    msy <- r * K / exp(1)
    bmsy <- K / exp(1)
    special_case <- "Fox model (m ~ 1)"
  } else if (abs(m - 2) < 1e-6) {
    msy <- r * K / 4
    bmsy <- K / 2
    special_case <- "Schaefer model (m ~ 2)"
  } else if (m <= 1) {
    if (warn_on_invalid_m) {
      warning("Shape parameter m <= 1 may give unrealistic reference points")
    }
    msy <- r * K / 4
    bmsy <- K / 2
    special_case <- "General Pella-Tomlinson model"
  } else {
    msy <- r * K * (m - 1)^((m - 1) / m) / m
    bmsy <- K * (m - 1)^(1 / m) / m
    special_case <- "General Pella-Tomlinson model"
  }

  if (bmsy > 0) {
    fmsy <- msy / bmsy
  } else {
    fmsy <- NA_real_
    warning("BMSY is zero or negative, cannot calculate FMSY")
  }

  target_reference_points <- NULL
  if (!is.null(biomass_target)) {
    baseline_info <- resolve_reference_point_baseline(parameters, baseline)
    target_reference_points <- calculate_target_reference_points(
      parameters = parameters,
      biomass_target = biomass_target,
      baseline_name = baseline_info$name,
      baseline_biomass = baseline_info$value
    )
  }

  list(
    msy = msy,
    bmsy = bmsy,
    fmsy = fmsy,
    target_reference_points = target_reference_points,
    special_case = special_case
  )
}

validate_biomass_target <- function(biomass_target) {
  if (is.null(biomass_target)) return(invisible(NULL))

  if (!is.numeric(biomass_target) || any(!is.finite(biomass_target))) {
    stop("biomass_target must be a numeric vector of finite biomass fractions")
  }
  if (any(biomass_target <= 0 | biomass_target > 1)) {
    stop("All biomass_target values must be > 0 and <= 1")
  }

  invisible(NULL)
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
#' @method print pt_reference_points
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

  if (!is.null(x$target_reference_points) && nrow(x$target_reference_points) > 0) {
    cat("User-Defined Biomass Targets:\n")
    for (i in seq_len(nrow(x$target_reference_points))) {
      row <- x$target_reference_points[i, ]
      cat(" ", row$biomass_name, ":", sprintf("%.1f tonnes", row$biomass), "\n")
      cat(" ", row$f_name, ":", sprintf("%.3f per year", row$fishing_mortality), "\n")
    }
    cat("\n")
  }

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

resolve_reference_point_baseline <- function(parameters, baseline) {
  param_names <- names(parameters)
  b0_idx <- grepl("^B0(\\.|$)", param_names)
  has_b0 <- any(b0_idx)

  if (baseline == "auto") {
    baseline <- if (has_b0) "B0" else "K"
  }

  if (baseline == "B0") {
    if (!has_b0) {
      stop("baseline = 'B0' requested but no fitted B0 parameter was found")
    }

    if ("B0" %in% param_names) {
      baseline_biomass <- unname(parameters[["B0"]])
    } else {
      baseline_biomass <- sum(unname(parameters[b0_idx]))
    }
  } else {
    baseline_biomass <- unname(parameters[["K"]])
  }

  list(name = baseline, value = baseline_biomass)
}

calculate_target_reference_points <- function(parameters, biomass_target, baseline_name, baseline_biomass) {
  r <- parameters[["r"]]
  K <- parameters[["K"]]
  m <- parameters[["m"]]

  target_biomass <- biomass_target * baseline_biomass
  target_f <- calculate_target_fishing_mortality(r, K, m, target_biomass)
  target_f[!is.finite(target_f) | target_f < 0] <- NA_real_

  data.frame(
    fraction = biomass_target,
    baseline = baseline_name,
    biomass_name = vapply(biomass_target, format_reference_point_label, character(1), prefix = "B", baseline_name = baseline_name),
    biomass = target_biomass,
    f_name = vapply(biomass_target, format_reference_point_label, character(1), prefix = "F", baseline_name = baseline_name),
    fishing_mortality = target_f,
    stringsAsFactors = FALSE
  )
}

format_reference_point_label <- function(fraction, prefix, baseline_name) {
  pct <- sub("\\.?0+$", "", sprintf("%.4f", 100 * fraction))
  paste0(prefix, "_", pct, "%", baseline_name)
}

calculate_target_fishing_mortality <- function(r, K, m, target_biomass) {
  if (length(m) == 1L && abs(m - 1) < 1e-6) {
    return(r * log(K / target_biomass))
  }

  out <- r / m * (1 - (target_biomass / K)^(m - 1))
  fox_idx <- abs(m - 1) < 1e-6
  if (any(fox_idx)) {
    out[fox_idx] <- r[fox_idx] * log(K[fox_idx] / target_biomass[fox_idx])
  }
  out
}

extract_named_reference_value <- function(reference_points, name) {
  base_names <- c(MSY = "msy", BMSY = "bmsy", FMSY = "fmsy")
  if (name %in% names(base_names)) {
    return(reference_points[[base_names[[name]]]])
  }

  target_df <- reference_points$target_reference_points
  if (is.null(target_df) || nrow(target_df) == 0) return(NULL)

  b_idx <- match(name, target_df$biomass_name)
  if (!is.na(b_idx)) return(target_df$biomass[b_idx])

  f_idx <- match(name, target_df$f_name)
  if (!is.na(f_idx)) return(target_df$fishing_mortality[f_idx])

  NULL
}

resolve_reference_point_defaults <- function(biomass_target,
                                             baseline,
                                             biomass_target_missing = FALSE,
                                             baseline_missing = FALSE) {
  defaults <- get_reference_point_defaults()

  if (biomass_target_missing) {
    biomass_target <- defaults$biomass_target
  }
  if (baseline_missing) {
    baseline <- defaults$baseline
  }

  validate_biomass_target(biomass_target)
  baseline <- match.arg(baseline, c("auto", "B0", "K"))

  list(
    biomass_target = biomass_target,
    baseline = baseline
  )
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

  if (!model_fit$fitted) {
    stop("Model must be fitted before extracting biomass estimates")
  }

  if (!"biomass" %in% names(model_fit$results)) {
    stop("No biomass estimates found in fitted model")
  }

  # Extract data
  model_years <- model_fit$data$years
  biomass_values <- model_fit$results$biomass

  # Create results matrix
  results <- data.frame(
    year = model_years,
    biomass = biomass_values
  )

  # Add harvest rates if available
  if ("harvest_rate" %in% names(model_fit$results)) {
    results$harvest_rate <- model_fit$results$harvest_rate
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
