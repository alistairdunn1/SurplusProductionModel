#' Model Fitting Functions for Pella-Tomlinson Surplus Production Model
#'
#' This file contains the main model fitting infrastructure for the
#' Pella-Tomlinson surplus production model using RTMB for optimization.
#'
#' @name model-fitting
NULL

#' Fit Pella-Tomlinson Surplus Production Model
#'
#' Fits a Pella-Tomlinson surplus production model to catch and CPUE data
#' using maximum likelihood estimation via RTMB.
#'
#' @param data List containing model data with elements:
#'   \describe{
#'     \item{cpue_data}{Data frame with columns: year, cpue}
#'     \item{catch_data}{Data frame with columns: year, catch}
#'   }
#' @param params_init Named list of starting parameter values (optional).
#'   If NULL, starting values are generated automatically.
#' @param options List of optimization options (optional):
#'   \describe{
#'     \item{silent}{Logical, suppress RTMB output (default: TRUE)}
#'     \item{control}{List of control parameters for nlminb}
#'     \item{validate_data}{Logical, run data validation (default: TRUE)}
#'   }
#'
#' @return ProductionModel object with fitted results
#'
#' @details
#' This function implements the complete model fitting workflow:
#' 1. Data validation (if enabled)
#' 2. Data preprocessing and alignment
#' 3. Starting value generation (if not provided)
#' 4. RTMB objective function creation
#' 5. Maximum likelihood optimization using nlminb
#' 6. Standard error calculation via sdreport
#' 7. Convergence checking and diagnostics
#' 8. Results packaging in ProductionModel object
#'
#' The model uses a state-space formulation with:
#' - Process error in biomass dynamics
#' - Observation error in CPUE
#' - Log-transformed parameters for positivity constraints
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' data_list <- list(
#'   cpue_data = data.frame(year = 2000:2020, cpue = rnorm(21, 1.5, 0.3)),
#'   catch_data = data.frame(year = 2000:2020, catch = rnorm(21, 1000, 100))
#' )
#'
#' # Fit model
#' model_fit <- fit_pella_tomlinson_model(data_list)
#'
#' # View results
#' print(model_fit)
#' }
#' @include rtmb-objective.R
#' @export
fit_pella_tomlinson_model <- function(data, params_init = NULL, options = list()) {
    
  # Set default options
  default_options <- list(
    silent = TRUE,
    control = list(eval.max = 1000, iter.max = 500),
    validate_data = TRUE
  )
  options <- modifyList(default_options, options)

  # Input validation
  if (!is.list(data)) {
    stop("Argument 'data' must be a list")
  }

  if (!all(c("cpue_data", "catch_data") %in% names(data))) {
    stop("Data must contain 'cpue_data' and 'catch_data' elements")
  }

  # Validate data structure
  if (!is.data.frame(data$cpue_data) || !all(c("year", "cpue") %in% names(data$cpue_data))) {
    stop("cpue_data must be a data.frame with 'year' and 'cpue' columns")
  }

  if (!is.data.frame(data$catch_data) || !all(c("year", "catch") %in% names(data$catch_data))) {
    stop("catch_data must be a data.frame with 'year' and 'catch' columns")
  }

  # Data validation using validation functions
  if (options$validate_data) {
    cpue_validation <- validate_cpue_data(
      cpue = data$cpue_data$cpue,
      years = data$cpue_data$year,
      return_details = TRUE
    )
    catch_validation <- validate_catch_data(
      catch = data$catch_data$catch,
      years = data$catch_data$year,
      return_details = TRUE
    )

    if (!cpue_validation$valid) {
      stop("CPUE data validation failed: ", paste(cpue_validation$errors, collapse = "; "))
    }

    if (!catch_validation$valid) {
      stop("Catch data validation failed: ", paste(catch_validation$errors, collapse = "; "))
    }
  }

  # Data preprocessing: align years and merge datasets
  processed_data <- preprocess_model_data(data$cpue_data, data$catch_data)

  # Generate starting values if not provided
  if (is.null(params_init)) {
    params_init <- generate_starting_values(processed_data)
    message("Generated starting parameter values automatically")
  } else {
    # Validate provided starting values
    required_params <- c("log_r", "log_K", "log_m", "log_q", "log_sigma_proc", "log_sigma_obs", "log_B0")
    if (!all(required_params %in% names(params_init))) {
      stop(
        "params_init must contain all required parameters: ",
        paste(required_params, collapse = ", ")
      )
    }
  }

  # Create simple objective function for now (non-RTMB)
  # TODO: Implement full RTMB integration in future version
  #obj_fun <- create_simple_objective(processed_data, params_init)
  
  # create TMB objective function
  cmb <- function(f, d) function(p) f(p, d)
  obj_fun <- MakeADFun(cmb(rtmb_objective, processed_data), params_init)

  # Optimize using nlminb
  start_time <- Sys.time()

  opt_result <- nlminb(
    start     = obj_fun$par, 
    objective = obj_fun$fn, 
    gradient  = obj_fun$gr,
    control   = options$control
  )

  end_time <- Sys.time()
  fitting_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Check convergence
  if (opt_result$convergence != 0) {
    warning(
      "Model did not converge. Convergence code: ", opt_result$convergence,
      ". Message: ", opt_result$message
    )
  }

  # Calculate standard errors (simplified for MVP)
  std_errors <- NULL
  hessian_valid <- FALSE

  # For now, skip standard error calculation in MVP
  # Future versions will use RTMB::sdreport for proper uncertainty estimation

  # Transform parameters back to natural scale
  names(opt_result$par) <- names(params_init)
  fitted_params <- transform_parameters_to_natural(opt_result$par)

  # Calculate model results
  model_results <- calculate_model_results(fitted_params, processed_data)

  # Calculate reference points
  ref_points <- calculate_reference_points_internal(fitted_params)

  # Package results
  results_list <- list(
    parameters = fitted_params,
    std_errors = std_errors,
    biomass = model_results$biomass,
    harvest_rate = model_results$harvest_rate,
    fitted_cpue = model_results$fitted_cpue,
    residuals = model_results$residuals,
    msy = ref_points$msy,
    bmsy = ref_points$bmsy,
    fmsy = ref_points$fmsy,
    likelihood = opt_result$objective,
    convergence = opt_result$convergence,
    convergence_message = opt_result$message,
    fitting_time = fitting_time,
    hessian_valid = hessian_valid,
    n_parameters = length(opt_result$par),
    n_observations = length(processed_data$cpue),
    aic = 2 * length(opt_result$par) + 2 * opt_result$objective,
    bic = log(length(processed_data$cpue)) * length(opt_result$par) + 2 * opt_result$objective
  )

  # Create ProductionModel object
  fitted_model <- new("ProductionModel",
    parameters = fitted_params,
    data = list(
      years = processed_data$years,
      catch = processed_data$catch,
      cpue = processed_data$cpue,
      effort = processed_data$catch / processed_data$cpue
    ),
    results = list(),
    fitted = FALSE,
    model_type = "pella_tomlinson",
    creation_date = Sys.time()
  )

  # Update results slot
  fitted_model@results <- results_list
  fitted_model@fitted <- TRUE

  return(fitted_model)
}

#' Preprocess Model Data
#'
#' Aligns and merges CPUE and catch data for model fitting.
#'
#' @param cpue_data Data frame with year and cpue columns
#' @param catch_data Data frame with year and catch columns
#'
#' @return List with aligned data vectors
#'
#' @keywords internal
preprocess_model_data <- function(cpue_data, catch_data) {
  # Ensure data is sorted by year
  cpue_data <- cpue_data[order(cpue_data$year), ]
  catch_data <- catch_data[order(catch_data$year), ]

  # Find common years
  common_years <- intersect(cpue_data$year, catch_data$year)

  if (length(common_years) < 3) {
    stop("Insufficient overlapping years between CPUE and catch data (minimum 3 required)")
  }

  # Filter to common years
  cpue_subset <- cpue_data[cpue_data$year %in% common_years, ]
  catch_subset <- catch_data[catch_data$year %in% common_years, ]

  # Merge datasets
  merged_data <- merge(cpue_subset, catch_subset, by = "year", all = FALSE)
  merged_data <- merged_data[order(merged_data$year), ]

  # Check for missing values
  if (any(is.na(merged_data$cpue))) {
    warning("Missing CPUE values detected and will be handled in likelihood")
  }

  if (any(is.na(merged_data$catch))) {
    stop("Missing catch values are not allowed")
  }

  # Return processed data
  list(
    years = merged_data$year,
    cpue = merged_data$cpue,
    catch = merged_data$catch
  )
}

#' Transform Parameters to Natural Scale
#'
#' Transforms log-transformed parameters back to natural scale.
#'
#' @param log_params Named vector of log-transformed parameters
#'
#' @return Named vector of parameters on natural scale
#'
#' @keywords internal
transform_parameters_to_natural <- function(log_params) {
  natural_params <- exp(log_params)

  # Rename to natural scale parameter names
  names(natural_params) <- gsub("^log_", "", names(natural_params))

  return(natural_params)
}

#' Calculate Model Results
#'
#' Calculate biomass trajectory, harvest rates, and fitted values.
#'
#' @param parameters Named vector of fitted parameters on natural scale
#' @param data List of model data
#'
#' @return List with model results
#'
#' @keywords internal
calculate_model_results <- function(parameters, data) {
  n_years <- length(data$years)
  biomass <- numeric(n_years)
  harvest_rate <- numeric(n_years)
  fitted_cpue <- numeric(n_years)

  # Extract parameters
  r <- parameters[["r"]]
  K <- parameters[["K"]]
  m <- parameters[["m"]]
  q <- parameters[["q"]]
  B0 <- parameters[["B0"]]

  # Calculate biomass trajectory
  biomass[1] <- B0

  for (t in 1:(n_years - 1)) {
    # Pella-Tomlinson production
    if (biomass[t] > 0 && K > 0 && m > 0) {
      production <- r * biomass[t] * (1 - (biomass[t] / K)^(m - 1)) / m
    } else {
      production <- 0
    }

    # Next year's biomass
    biomass[t + 1] <- biomass[t] + production - data$catch[t]
    biomass[t + 1] <- max(biomass[t + 1], 0.01) # Prevent negative biomass
  }

  # Calculate harvest rates and fitted CPUE
  for (t in 1:n_years) {
    harvest_rate[t] <- data$catch[t] / biomass[t]
    fitted_cpue[t] <- q * biomass[t]
  }

  # Calculate residuals
  residuals <- log(data$cpue) - log(fitted_cpue)
  residuals[is.infinite(residuals)] <- NA # Handle missing CPUE

  list(
    biomass = biomass,
    harvest_rate = harvest_rate,
    fitted_cpue = fitted_cpue,
    residuals = residuals
  )
}

#' Calculate Reference Points (Internal)
#'
#' Calculate biological reference points from fitted parameters.
#'
#' @param parameters Named vector of fitted parameters
#'
#' @return List with reference points
#'
#' @keywords internal
calculate_reference_points_internal <- function(parameters) {
  r <- parameters[["r"]]
  K <- parameters[["K"]]
  m <- parameters[["m"]]

  # Calculate reference points using Pella-Tomlinson formulas
  if (m <= 1) {
    # Handle edge case
    msy <- r * K / 4 # Approximate for m close to 1
    bmsy <- K / 2
  } else {
    msy <- r * K * (m - 1)^((m - 1) / m) / m
    bmsy <- K * (m - 1)^(1 / m) / m
  }

  fmsy <- msy / bmsy

  list(
    msy = msy,
    bmsy = bmsy,
    fmsy = fmsy
  )
}
