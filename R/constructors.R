#' Constructor and Methods for ProductionModel Objects
#'
#' This file contains the constructor function and S3 methods for creating
#' and manipulating ProductionModel objects.
#'
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
#' @return A ProductionModel object (S3 list)
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

  # Validate positive CPUE
  if (any(cpue <= 0)) {
    stop("cpue data must be positive finite numbers")
  }

  # Convert years to integer if needed
  years <- as.integer(years)

  # Set default parameters if not provided
  if (is.null(parameters)) {
    parameters <- c(
      r = 0.1,
      K = 50000,
      m = 2.0,
      q = 0.001,
      sigma_proc = 0.1,
      sigma_obs = 0.2
    )
  } else {
    checkmate::assert_numeric(parameters, any.missing = FALSE, finite = TRUE)
    checkmate::assert_names(names(parameters),
      must.include = c("r", "K", "m", "q", "sigma_proc", "sigma_obs")
    )

    if (parameters["r"] <= 0) stop("Parameter 'r' must be positive")
    if (parameters["K"] <= 0) stop("Parameter 'K' must be positive")
    if (parameters["m"] <= 0) stop("Parameter 'm' must be positive")
    if (parameters["q"] <= 0) stop("Parameter 'q' must be positive")
    if (parameters["sigma_proc"] <= 0) stop("Parameter 'sigma_proc' must be positive")
    if (parameters["sigma_obs"] <= 0) stop("Parameter 'sigma_obs' must be positive")
  }

  # Build the S3 object
  obj <- list(
    parameters = parameters,
    data = list(
      years = years,
      catch = catch,
      cpue = cpue,
      effort = effort
    ),
    results = list(),
    fitted = FALSE,
    model_type = "pella_tomlinson",
    creation_date = Sys.time()
  )
  class(obj) <- "ProductionModel"
  obj
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Get Model Parameters
#'
#' Extract parameter values from a ProductionModel object.
#'
#' @param object A ProductionModel object
#' @param ... Additional arguments (not used)
#' @return Named numeric vector of parameters
#' @export
parameters.ProductionModel <- function(object, ...) {
  object$parameters
}

#' Get Model Data
#'
#' Extract data from a ProductionModel object.
#'
#' @param object A ProductionModel object
#' @param ... Additional arguments (not used)
#' @return List containing model data (years, catch, cpue, effort)
#' @export
model_data.ProductionModel <- function(object, ...) {
  object$data
}

#' Get Model Results
#'
#' Extract fitting results from a ProductionModel object.
#'
#' @param object A ProductionModel object
#' @param ... Additional arguments (not used)
#' @return List containing model results (empty if not fitted)
#' @export
results.ProductionModel <- function(object, ...) {
  object$results
}

#' Check if Model is Fitted
#'
#' Check whether a ProductionModel object has been fitted.
#'
#' @param object A ProductionModel object
#' @param ... Additional arguments (not used)
#' @return Logical indicating if model has been fitted
#' @export
fitted.ProductionModel <- function(object, ...) {
  object$fitted
}

#' Print Method for ProductionModel
#'
#' Print summary information about a ProductionModel object.
#'
#' @param x A ProductionModel object
#' @param ... Additional arguments (currently unused)
#' @return Invisibly returns the object
#' @export
print.ProductionModel <- function(x, biomass_target = NULL, baseline = c("auto", "B0", "K"), ...) {
  defaults <- resolve_reference_point_defaults(
    biomass_target = biomass_target,
    baseline = baseline,
    biomass_target_missing = missing(biomass_target),
    baseline_missing = missing(baseline)
  )
  biomass_target <- defaults$biomass_target
  baseline <- defaults$baseline

  cat("Production Model Object\n")
  cat("======================\n\n")

  cat("Model Type:", x$model_type, "\n")
  cat("Created:", format(x$creation_date, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Fitted:", ifelse(x$fitted, "Yes", "No"), "\n\n")

  if (length(x$data) > 0) {
    cat("Data Summary:\n")
    yrs <- x$data$years
    cat("  Years:", range(yrs)[1], "-", range(yrs)[2],
        "(", length(yrs), "years)\n")
    # Area / label counts
    areas  <- x$data$areas
    n_areas <- if (!is.null(areas)) length(areas) else 1L
    cat("  Areas:", n_areas)
    if (n_areas > 1) cat(" (", paste(areas, collapse = ", "), ")")
    cat("\n")
    cpue <- x$data$cpue
    if (is.array(cpue) && length(dim(cpue)) == 3) {
      labs <- dimnames(cpue)[[3]]
      cat("  CPUE indices:", length(labs), "(", paste(labs, collapse = ", "), ")\n")
    }
    ctch <- x$data$catch
    if (is.matrix(ctch)) {
      cat("  Total catch range:", sprintf("%.0f - %.0f tonnes",
          min(rowSums(ctch)), max(rowSums(ctch))), "\n")
    } else {
      cat("  Catch range:", sprintf("%.0f - %.0f tonnes",
          min(ctch), max(ctch)), "\n")
    }
    cat("\n")
  }

  if (x$fitted && length(x$results) > 0) {
    cat("Convergence:",
        ifelse(x$results$convergence == 0, "Success", "FAILED"),
        "\n")
    cat("Neg. log-likelihood:", sprintf("%.2f", x$results$likelihood), "\n")
    if (!is.null(x$results$aic)) {
      cat("AIC:", sprintf("%.2f", x$results$aic), "\n")
    }
    cat("Fitting time:", sprintf("%.2f s", x$results$fitting_time), "\n\n")

    cat("Key Parameters:\n")
    p <- x$parameters
    for (nm in c("r", "K", "m", "sigma_obs")) {
      if (nm %in% names(p)) cat("  ", nm, "=", sprintf("%.4g", p[[nm]]), "\n")
    }
    cat("\n")

    if (!is.null(x$results$msy)) {
      cat("Reference Points:\n")
      cat("  MSY  =", sprintf("%.1f tonnes", x$results$msy), "\n")
      cat("  BMSY =", sprintf("%.1f tonnes", x$results$bmsy), "\n")
      cat("  FMSY =", sprintf("%.4f", x$results$fmsy), "\n")
      if (!is.null(biomass_target)) {
        ref_points <- calculate_reference_points(x, biomass_target = biomass_target, baseline = baseline)
        .print_target_reference_points_block(ref_points$target_reference_points, indent = "  ")
      }
    }
  } else if (length(x$parameters) > 0) {
    cat("Parameters:\n")
    for (i in seq_along(x$parameters)) {
      cat("  ", names(x$parameters)[i], ":", sprintf("%.4f", x$parameters[i]), "\n")
    }
  }
  cat("\n")

  invisible(x)
}


#' Summary Method for ProductionModel
#'
#' Produces a detailed summary of a fitted ProductionModel including
#' parameter estimates with standard errors, information criteria,
#' and per-area / per-index diagnostics.
#'
#' @param object A ProductionModel object
#' @param ... Additional arguments (currently unused)
#' @return Invisibly returns a list of summary components
#' @export
summary.ProductionModel <- function(object, biomass_target = NULL, baseline = c("auto", "B0", "K"), ...) {
  defaults <- resolve_reference_point_defaults(
    biomass_target = biomass_target,
    baseline = baseline,
    biomass_target_missing = missing(biomass_target),
    baseline_missing = missing(baseline)
  )
  biomass_target <- defaults$biomass_target
  baseline <- defaults$baseline

  if (!object$fitted) {
    cat("Unfitted ProductionModel -- no summary available.\n")
    return(invisible(NULL))
  }

  results <- object$results
  params  <- object$parameters
  se      <- results$std_errors

  cat("Summary: Pella-Tomlinson Surplus Production Model\n")
  cat("=================================================\n\n")

  # ---- Data overview ----
  yrs    <- object$data$years
  areas  <- object$data$areas
  n_areas <- if (!is.null(areas)) length(areas) else 1L
  cpue   <- object$data$cpue
  has_labels <- is.array(cpue) && length(dim(cpue)) == 3

  cat("Data:\n")
  cat("  Years:", min(yrs), "-", max(yrs), "(", length(yrs), ")\n")
  cat("  Areas:", n_areas, "\n")
  if (has_labels) {
    cat("  CPUE indices:", length(dimnames(cpue)[[3]]),
        "(", paste(dimnames(cpue)[[3]], collapse = ", "), ")\n")
  }
  n_obs <- results$n_observations
  n_par <- results$n_parameters
  cat("  Observations:", n_obs, "\n")
  cat("  Estimated parameters:", n_par, "\n\n")

  # ---- Parameters with SEs ----
  cat("Parameter Estimates:\n")
  # Build param table
  pnames <- names(params)
  ptable <- data.frame(
    Parameter = pnames,
    Estimate  = unname(params),
    stringsAsFactors = FALSE
  )
  if (!is.null(se)) {
    # Match SEs to natural-scale params via log-scale names
    log_names <- paste0("log_", gsub("\\.", "_", pnames))
    se_vals <- numeric(length(pnames))
    for (i in seq_along(pnames)) {
      # delta-method: SE(exp(x)) ~ exp(x) * SE(x)
      se_idx <- match(log_names[i], names(se))
      if (!is.na(se_idx)) {
        se_vals[i] <- params[i] * se[se_idx]
      } else {
        se_vals[i] <- NA
      }
    }
    ptable$SE <- se_vals
    ptable$CV <- ifelse(ptable$Estimate > 0, ptable$SE / ptable$Estimate, NA)
  }
  print(ptable, digits = 4, row.names = FALSE)
  cat("\n")

  # ---- Fit statistics ----
  cat("Fit Statistics:\n")
  cat("  Neg. log-likelihood:", sprintf("%.4f", results$likelihood), "\n")
  cat("  AIC:", sprintf("%.2f", results$aic), "\n")
  cat("  BIC:", sprintf("%.2f", results$bic), "\n")
  cat("  Hessian positive-definite:", results$hessian_valid, "\n")
  if (results$process_noise) cat("  Process noise: enabled (Laplace)\n")
  cat("\n")

  # ---- Reference points ----
  if (!is.null(results$msy)) {
    cat("Reference Points:\n")
    cat(sprintf("  MSY  = %.1f tonnes\n", results$msy))
    cat(sprintf("  BMSY = %.1f tonnes\n", results$bmsy))
    cat(sprintf("  FMSY = %.4f\n", results$fmsy))
    if (!is.null(biomass_target)) {
      ref_points <- calculate_reference_points(object, biomass_target = biomass_target, baseline = baseline)
      .print_target_reference_points_block(ref_points$target_reference_points, indent = "  ")
    }
    # Current status
    bio <- results$biomass
    terminal_bio <- if (is.matrix(bio)) sum(bio[nrow(bio), ]) else bio[length(bio)]
    cat(sprintf("  B_terminal / BMSY = %.2f\n", terminal_bio / results$bmsy))
    cat("\n")
  }

  # ---- Residual statistics ----
  res <- results$residuals
  if (!is.null(res)) {
    cat("Residual Diagnostics:\n")
    if (is.array(res) && length(dim(res)) == 3) {
      for (a in dimnames(res)[[2]]) {
        for (l in dimnames(res)[[3]]) {
          rv <- res[, a, l]
          rv <- rv[is.finite(rv)]
          if (length(rv) >= 3) {
            cat(sprintf("  [%s, %s] mean=%.3f  sd=%.3f  n=%d\n",
                        a, l, mean(rv), sd(rv), length(rv)))
          }
        }
      }
    } else if (is.matrix(res)) {
      for (a in colnames(res)) {
        rv <- res[, a]
        rv <- rv[is.finite(rv)]
        if (length(rv) >= 3) {
          cat(sprintf("  [%s] mean=%.3f  sd=%.3f  n=%d\n",
                      a, mean(rv), sd(rv), length(rv)))
        }
      }
    } else {
      rv <- res[is.finite(res)]
      cat(sprintf("  mean=%.3f  sd=%.3f  n=%d\n", mean(rv), sd(rv), length(rv)))
    }
  }

  invisible(list(parameters = ptable, results = results))
}

.print_target_reference_points_block <- function(target_reference_points, indent = "") {
  if (is.null(target_reference_points) || nrow(target_reference_points) == 0) {
    return(invisible(NULL))
  }

  cat(sprintf("%sUser-Defined Biomass Targets:\n", indent))
  for (i in seq_len(nrow(target_reference_points))) {
    row <- target_reference_points[i, ]
    cat(sprintf("%s  %s = %.1f tonnes\n", indent, row$biomass_name, row$biomass))
    cat(sprintf("%s  %s = %.4f\n", indent, row$f_name, row$fishing_mortality))
  }

  invisible(NULL)
}
