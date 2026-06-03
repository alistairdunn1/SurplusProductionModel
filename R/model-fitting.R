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
#'     \item{cpue_data}{Data frame with columns: year, cpue, optional: area}
#'     \item{catch_data}{Data frame with columns: year, catch, optional: area}
#'     \item{movement}{optional list with movement inputs: distance_matrix (area by area matrix), attractiveness (named numeric by area)}
#'   }
#' @param params_init Named list of starting parameter values (optional).
#'   If NULL, starting values are generated automatically.
#' @param options List of optimization options (optional):
#'   \describe{
#'     \item{silent}{Logical, suppress RTMB output (default: TRUE)}
#'     \item{show_starting_values_message}{Logical, print a message when
#'       starting values are generated automatically (default: TRUE).}
#'     \item{control}{List of control parameters for nlminb}
#'     \item{validate_data}{Logical, run data validation (default: TRUE)}
#'     \item{process_noise}{Logical, enable state-space process deviations as
#'       RTMB random effects integrated with the Laplace approximation
#'       (default: FALSE).}
#'     \item{process_error_structure}{Character, one of \code{"iid"}
#'       (default) or \code{"ar1"}. Applies only when
#'       \code{process_noise = TRUE}.}
#'     \item{env_covariates}{Optional character vector of environmental
#'       covariate column names expected in \code{data$env_data}. When
#'       supplied, linear covariate effects are included in biomass dynamics.}
#'     \item{env_lag}{Non-negative integer lag (years) applied to
#'       environmental covariates (default: 0).}
#'     \item{env_scale}{Logical, center and scale covariates before fitting
#'       (default: TRUE).}
#'     \item{n_starts}{Integer, number of random restarts (default: 1).
#'       When > 1, the optimizer is run from \code{n_starts} different
#'       starting vectors (the original plus jittered versions) and the
#'       run with the lowest objective is retained.}
#'     \item{jitter_sd}{Numeric, standard deviation of log-normal jitter
#'       applied to starting values for multi-start (default: 0.2).}
#'     \item{priors}{Optional named list of priors on model parameters.
#'       Each element is a list with \code{dist} and distribution-specific
#'       fields. Supported distributions are:\cr
#'       \code{normal}: \code{list(dist = "normal", mean = ..., sd = ...)}\cr
#'       \code{lognormal}: \code{list(dist = "lognormal", meanlog = ..., sdlog = ...)}\cr
#'       Parameter names can be supplied in natural form (e.g. \code{"r"},
#'       \code{"K"}) or log form (e.g. \code{"log_r"}, \code{"log_K"}).}
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
#' By default the model uses deterministic biomass dynamics with observation
#' error in CPUE. When \code{options$process_noise = TRUE}, the fitter adds a
#' state-space formulation with process deviations in biomass dynamics as RTMB
#' random effects.
#'
#' In both modes, parameters are log-transformed where needed to enforce
#' positivity constraints.
#'
#' @references Pella, J. J.; Tomlinson, P. K. (1969). A generalised stock production model. Inter-American Tropical Tuna Commission Bulletin 13, 419-496.
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
#' # Fit model with informative priors on r and K
#' fit_pr <- fit_pella_tomlinson_model(
#'   data_list,
#'   options = list(
#'     priors = list(
#'       r = list(dist = "lognormal", meanlog = log(0.2), sdlog = 0.4),
#'       K = list(dist = "lognormal", meanlog = log(6000), sdlog = 0.5)
#'     )
#'   )
#' )
#'
#' # View results
#' print(model_fit)
#' }
#'
#' @export
fit_pella_tomlinson_model <- function(data, params_init = NULL, options = list()) {
  # Set default options
  default_options <- list(
    silent = TRUE,
    show_starting_values_message = TRUE,
    control = list(eval.max = 1000, iter.max = 500),
    validate_data = TRUE,
    process_noise = FALSE,
    process_error_structure = "iid",
    env_covariates = NULL,
    env_lag = 0,
    env_scale = TRUE,
    priors = NULL
  )
  options <- modifyList(default_options, options)

  use_process_noise <- isTRUE(options$process_noise)
  process_error_structure <- match.arg(
    as.character(options$process_error_structure %||% "iid"),
    c("iid", "ar1")
  )

  env_covariates <- options$env_covariates
  if (is.null(env_covariates)) {
    env_covariates <- character(0)
  }
  if (!is.character(env_covariates)) {
    stop("options$env_covariates must be NULL or a character vector")
  }
  env_covariates <- unique(env_covariates[nzchar(env_covariates)])

  env_lag <- as.integer(options$env_lag %||% 0)
  if (!is.finite(env_lag) || length(env_lag) != 1 || env_lag < 0) {
    stop("options$env_lag must be a single non-negative integer")
  }
  env_scale <- isTRUE(options$env_scale)

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

  # Attach movement inputs (if provided) to processed data
  if (is.list(data$movement)) {
    # Accept distance_matrix, attractiveness, movement_rate, decay
    if (!is.null(data$movement$distance_matrix)) processed_data$distance_matrix <- data$movement$distance_matrix
    if (!is.null(data$movement$attractiveness)) processed_data$attractiveness <- data$movement$attractiveness
    if (!is.null(data$movement$movement_rate)) processed_data$movement_rate <- data$movement$movement_rate
    if (!is.null(data$movement$decay)) processed_data$decay <- data$movement$decay
  }

  # Optional environmental covariates used in biomass transition equation
  env_info <- .prepare_environmental_covariates(
    env_data = data$env_data,
    years = processed_data$years,
    areas = processed_data$areas %||% "A1",
    covariate_names = env_covariates,
    lag = env_lag,
    scale_covariates = env_scale
  )
  if (!is.null(env_info)) {
    processed_data$env_array <- env_info$env_array
    processed_data$env_covariate_names <- env_info$covariate_names
    processed_data$env_scaling <- env_info$scaling
  }

  # Generate starting values if not provided
  if (is.null(params_init)) {
    params_init <- generate_starting_values(processed_data)
    if (isTRUE(options$show_starting_values_message)) {
      message("Generated starting parameter values automatically")
    }
  } else {
    # Validate provided starting values
    # Required globals
    required_params <- c("log_r", "log_K", "log_m", "log_sigma_proc", "log_sigma_obs")
    # Per-area and per-index params when multi-area / multi-index
    if (!is.null(processed_data$areas)) {
      area_suffix <- paste0(".", processed_data$areas)
      if (!is.null(processed_data$labels)) {
        # q per area and label
        label_suffix <- paste0(".", processed_data$labels)
        q_names <- as.vector(outer(paste0("log_q", area_suffix), label_suffix, paste0))
        required_params <- c(required_params, q_names, paste0("log_B0", area_suffix))
      } else {
        required_params <- c(required_params, paste0("log_q", area_suffix), paste0("log_B0", area_suffix))
      }
    } else {
      required_params <- c(required_params, "log_q", "log_B0")
    }
    if (!all(required_params %in% names(params_init))) {
      stop(
        "params_init must contain all required parameters: ",
        paste(required_params, collapse = ", ")
      )
    }
  }

  # ---- Build RTMB data environment and parameter list ------------------

  # Normalise parameter names: RTMB getAll() uses R variable names, so we use

  # underscores (log_q_A1) rather than dots (log_q.A1).
  rename_dots <- function(x) {
    nms <- names(x)
    # Replace dots in suffixes like log_q.A1 -> log_q_A1
    nms <- gsub("\\.", "_", nms)
    names(x) <- nms
    x
  }

  params_init <- rename_dots(params_init)

  # Build data environment for the RTMB objective
  n_years <- length(processed_data$years)
  areas <- if (!is.null(processed_data$areas)) processed_data$areas else "A1"
  n_areas <- length(areas)
  has_labels <- !is.null(processed_data$labels)

  # Normalise single-area parameter names:
  # log_q -> log_q_A1, log_B0 -> log_B0_A1  (objective always uses area suffix)
  if (n_areas == 1 && !has_labels) {
    a <- areas[1]
    q_key <- paste0("log_q_", a)
    b0_key2 <- paste0("log_B0_", a)
    if ("log_q" %in% names(params_init) && !q_key %in% names(params_init)) {
      params_init[[q_key]] <- params_init[["log_q"]]
      params_init[["log_q"]] <- NULL
    }
    if ("log_B0" %in% names(params_init) && !b0_key2 %in% names(params_init)) {
      params_init[[b0_key2]] <- params_init[["log_B0"]]
      params_init[["log_B0"]] <- NULL
    }
  }

  # Add default starting values for optional environmental effects and AR1.
  if (!is.null(processed_data$env_covariate_names)) {
    for (cov_name in processed_data$env_covariate_names) {
      beta_name <- paste0("beta_", cov_name)
      if (!beta_name %in% names(params_init)) {
        params_init[[beta_name]] <- 0
      }
    }
  }
  if (use_process_noise && identical(process_error_structure, "ar1") && !"theta_rho" %in% names(params_init)) {
    params_init[["theta_rho"]] <- atanh(0.2)
  }

  # catch as matrix [n_years x n_areas] always
  catch_mat <- if (is.matrix(processed_data$catch)) {
    processed_data$catch
  } else {
    matrix(processed_data$catch, ncol = 1, dimnames = list(NULL, areas))
  }

  # cpue as matrix or 3-d array
  if (has_labels) {
    cpue_obs <- processed_data$cpue # [year x area x label]
  } else {
    cpue_obs <- if (is.matrix(processed_data$cpue)) {
      processed_data$cpue
    } else {
      matrix(processed_data$cpue, ncol = 1, dimnames = list(NULL, areas))
    }
  }

  rtmb_data <- list(
    n_years = n_years,
    n_areas = n_areas,
    areas = areas,
    catch_mat = catch_mat,
    cpue_obs = cpue_obs,
    labels = if (has_labels) processed_data$labels else NULL,
    process_error_structure = process_error_structure
  )

  if (!is.null(processed_data$env_array)) {
    rtmb_data$env_array <- processed_data$env_array
    rtmb_data$env_covariate_names <- processed_data$env_covariate_names
  }

  # Attach movement data
  if (!is.null(processed_data$movement_rate)) {
    rtmb_data$movement_rate <- processed_data$movement_rate
    rtmb_data$decay <- processed_data$decay
    rtmb_data$distance_matrix <- processed_data$distance_matrix
    rtmb_data$attractiveness <- processed_data$attractiveness
  }

  # ---- Build parameter list for MakeADFun ---------------------------
  # RTMB parameters must be a named list (not a flat vector).
  # Each element is a scalar, vector, or matrix.
  # Strip element names first to prevent unlist from concatenating
  # outer and inner names (e.g. log_q_A1 = c(A1=val) -> log_q_A1.A1).
  params_init <- lapply(params_init, unname)
  rtmb_parms <- as.list(unlist(params_init))

  # Handle process noise / random effects
  if (use_process_noise) {
    # Add proc_dev random effects [n_years-1 x n_areas]
    rtmb_parms$proc_dev <- matrix(0, n_years - 1, n_areas)
    # Ensure log_sigma_proc is present
    if (is.null(rtmb_parms$log_sigma_proc)) {
      rtmb_parms$log_sigma_proc <- log(0.05)
    }
  }

  # ---- Build map for fixed parameters --------------------------------
  # options$fixed_params is a named list of log-scale values.
  # In RTMB, we set the parameter to the desired value in rtmb_parms,
  # and use map = list(param_name = factor(NA)) to fix it.
  fixed <- options$fixed_params
  map_list <- list()
  if (!is.null(fixed) && length(fixed) > 0) {
    fix_names <- gsub("\\.", "_", names(fixed))
    for (fn in fix_names) {
      if (!fn %in% names(rtmb_parms)) {
        stop("fixed_params name '", fn, "' not found in parameter list")
      }
      rtmb_parms[[fn]] <- unname(fixed[[gsub("_", ".", fn)]])
      # If parameter is NA, try the underscore version
      if (is.null(rtmb_parms[[fn]])) rtmb_parms[[fn]] <- unname(fixed[[fn]])
      map_list[[fn]] <- factor(NA)
    }
  }

  priors <- .prepare_model_priors(
    priors = options$priors,
    rtmb_parms = rtmb_parms,
    areas = areas
  )
  if (!is.null(priors)) {
    rtmb_data$priors <- priors
  }

  # ---- Create RTMB AD object ----------------------------------------
  rtmb_obj_fn <- create_rtmb_objective(rtmb_data)

  random_effects <- if (use_process_noise) "proc_dev" else NULL

  obj <- RTMB::MakeADFun(
    func       = rtmb_obj_fn,
    parameters = rtmb_parms,
    random     = random_effects,
    map        = map_list,
    silent     = options$silent
  )

  # ---- Optimise (with optional multi-start) ----------------------------
  n_starts <- max(1L, as.integer(options$n_starts %||% 1L))
  jitter_sd <- as.numeric(options$jitter_sd %||% 0.2)

  start_time <- Sys.time()

  best_opt <- NULL
  best_nll <- Inf

  for (istart in seq_len(n_starts)) {
    if (istart == 1L) {
      start_par <- obj$par
    } else {
      # Jitter the starting parameters on the log scale
      start_par <- obj$par + rnorm(length(obj$par), 0, jitter_sd)
    }
    this_opt <- tryCatch(
      nlminb(
        start     = start_par,
        objective = obj$fn,
        gradient  = obj$gr,
        control   = options$control
      ),
      error = function(e) NULL
    )
    if (!is.null(this_opt) && is.finite(this_opt$objective) &&
      this_opt$objective < best_nll) {
      best_nll <- this_opt$objective
      best_opt <- this_opt
    }
  }

  # Fall back to the first run if all starts failed
  if (is.null(best_opt)) {
    best_opt <- nlminb(
      start     = obj$par,
      objective = obj$fn,
      gradient  = obj$gr,
      control   = options$control
    )
  }
  opt_result <- best_opt

  end_time <- Sys.time()
  fitting_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Check convergence
  if (opt_result$convergence != 0) {
    warning(
      "Model did not converge. Convergence code: ", opt_result$convergence,
      ". Message: ", opt_result$message
    )
  }

  # ---- sdreport for standard errors -----------------------------------
  std_errors <- NULL
  hessian_valid <- FALSE
  sdr <- NULL

  tryCatch(
    {
      sdr <- RTMB::sdreport(obj)
      hessian_valid <- sdr$pdHess
      # Fixed parameter SEs
      summ_fixed <- summary(sdr, "fixed")
      std_errors <- setNames(summ_fixed[, "Std. Error"], rownames(summ_fixed))
    },
    error = function(e) {
      warning("sdreport failed: ", e$message)
    }
  )

  # ---- Extract results ------------------------------------------------
  # First, transform parameters back to natural scale (needed for
  # calculate_model_results and reference points).
  # Reconstruct full log-parameter vector (including fixed)
  full_log_par <- unlist(rtmb_parms)
  # Overwrite estimated values from optimised result
  est_names <- names(opt_result$par)
  for (nm in est_names) {
    full_log_par[nm] <- opt_result$par[nm]
  }
  # Keep only scalar parameters (exclude proc_dev, etc.)
  scalar_par_names <- grep("^(log_r|log_K|log_m|log_sigma|log_q|log_B0)", names(full_log_par), value = TRUE)
  scalar_par_names <- unique(c(
    scalar_par_names,
    grep("^(beta_|theta_rho)$", names(full_log_par), value = TRUE)
  ))
  log_par_scalar <- full_log_par[scalar_par_names]
  # Rename back to dot-separated for backward compatibility.
  # Only convert area/label suffixes, NOT core parameter names like log_K.
  for (a in areas) {
    names(log_par_scalar) <- gsub(paste0("_", a), paste0(".", a),
      names(log_par_scalar),
      fixed = TRUE
    )
  }
  if (has_labels) {
    for (l in processed_data$labels) {
      names(log_par_scalar) <- gsub(paste0("_", l), paste0(".", l),
        names(log_par_scalar),
        fixed = TRUE
      )
    }
  }
  fitted_params <- transform_parameters_to_natural(log_par_scalar)

  # Recompute derived quantities in plain R from fitted parameters.
  # This avoids RTMB REPORT/advector complications (NaN promotion,
  # dimension loss) while still using ADREPORT for sdreport SEs.
  model_results <- calculate_model_results(fitted_params, processed_data)
  B_est <- model_results$biomass
  hr_est <- model_results$harvest_rate
  fc_est <- model_results$fitted_cpue
  residuals <- model_results$residuals

  # SE for derived quantities from ADREPORT (biomass SEs)
  biomass_se <- NULL
  if (!is.null(sdr) && hessian_valid) {
    tryCatch(
      {
        summ_report <- summary(sdr, "report")
        # Extract B_mat rows
        b_rows <- grepl("^B_mat$", rownames(summ_report))
        if (any(b_rows)) {
          biomass_se <- summ_report[b_rows, "Std. Error"]
          if (n_areas == 1 && !has_labels) {
            biomass_se <- as.numeric(biomass_se)
          } else {
            biomass_se <- matrix(biomass_se,
              nrow = n_years, ncol = n_areas,
              dimnames = list(as.character(processed_data$years), areas)
            )
          }
        }
      },
      error = function(e) NULL
    )
  }

  # Reference points
  ref_points <- calculate_reference_points_internal(fitted_params)

  # Count observations and effective parameters
  n_obs <- if (has_labels) sum(!is.na(cpue_obs)) else sum(!is.na(cpue_obs))
  n_est_pars <- length(opt_result$par)

  # Package results
  results_list <- list(
    parameters = fitted_params,
    std_errors = std_errors,
    biomass = B_est,
    biomass_se = biomass_se,
    harvest_rate = hr_est,
    fitted_cpue = fc_est,
    residuals = residuals,
    msy = ref_points$msy,
    bmsy = ref_points$bmsy,
    fmsy = ref_points$fmsy,
    opt_par = opt_result$par,
    likelihood = opt_result$objective,
    convergence = opt_result$convergence,
    convergence_message = opt_result$message,
    fitting_time = fitting_time,
    hessian_valid = hessian_valid,
    n_parameters = n_est_pars,
    n_observations = n_obs,
    aic = 2 * n_est_pars + 2 * opt_result$objective,
    bic = log(n_obs) * n_est_pars + 2 * opt_result$objective,
    process_noise = use_process_noise,
    process_error_structure = if (use_process_noise) process_error_structure else "none",
    rho = if ("rho" %in% names(fitted_params)) fitted_params[["rho"]] else NA_real_,
    env_effects = {
      env_idx <- grepl("^beta[._]", names(fitted_params))
      if (any(env_idx)) fitted_params[env_idx] else NULL
    },
    env_scaling = processed_data$env_scaling %||% NULL,
    priors = priors,
    rtmb_obj = obj,
    sdreport = sdr
  )

  # Create ProductionModel S3 object
  effort <- {
    if (is.array(processed_data$cpue) && length(dim(processed_data$cpue)) == 3) {
      cpue_mean <- apply(processed_data$cpue, c(1, 2), function(x) mean(x, na.rm = TRUE))
      processed_data$catch / cpue_mean
    } else if (is.matrix(processed_data$cpue)) {
      processed_data$catch / processed_data$cpue
    } else {
      processed_data$catch / processed_data$cpue
    }
  }

  fitted_model <- structure(
    list(
      parameters = fitted_params,
      data = list(
        years = processed_data$years,
        areas = processed_data$areas,
        catch = processed_data$catch,
        cpue = processed_data$cpue,
        effort = effort,
        env_scaling = processed_data$env_scaling %||% NULL
      ),
      results = results_list,
      fitted = TRUE,
      model_type = "pella_tomlinson",
      creation_date = Sys.time()
    ),
    class = "ProductionModel"
  )

  return(fitted_model)
}

.prepare_model_priors <- function(priors, rtmb_parms, areas) {
  if (is.null(priors)) {
    return(NULL)
  }

  if (!is.list(priors) || is.null(names(priors)) || any(names(priors) == "")) {
    stop("options$priors must be a named list of prior specifications")
  }

  normalized <- vector("list", length(priors))
  idx <- 0L

  for (param_name in names(priors)) {
    spec <- priors[[param_name]]
    if (!is.list(spec)) {
      stop("Prior for '", param_name, "' must be supplied as a list")
    }

    fit_param <- .normalize_prior_parameter_name(param_name, rtmb_parms, areas)
    dist <- spec$dist %||% spec$distribution
    if (is.null(dist) || !is.character(dist) || length(dist) != 1) {
      stop("Prior for '", param_name, "' must supply 'dist' or 'distribution'")
    }
    dist <- tolower(dist)

    idx <- idx + 1L
    if (dist %in% c("normal", "gaussian")) {
      mean <- spec$mean
      sd <- spec$sd
      if (!is.numeric(mean) || length(mean) != 1 || !is.finite(mean)) {
        stop("Normal prior for '", param_name, "' must supply a finite scalar 'mean'")
      }
      if (!is.numeric(sd) || length(sd) != 1 || !is.finite(sd) || sd <= 0) {
        stop("Normal prior for '", param_name, "' must supply a positive finite scalar 'sd'")
      }
      normalized[[idx]] <- list(param = fit_param, dist = "normal", mean = mean, sd = sd)
    } else if (dist == "lognormal") {
      meanlog <- spec$meanlog
      sdlog <- spec$sdlog
      if (!grepl("^log_", fit_param)) {
        stop("Lognormal priors are only supported for log-scale parameters; use a normal prior for '", param_name, "'")
      }
      if (!is.numeric(meanlog) || length(meanlog) != 1 || !is.finite(meanlog)) {
        stop("Lognormal prior for '", param_name, "' must supply a finite scalar 'meanlog'")
      }
      if (!is.numeric(sdlog) || length(sdlog) != 1 || !is.finite(sdlog) || sdlog <= 0) {
        stop("Lognormal prior for '", param_name, "' must supply a positive finite scalar 'sdlog'")
      }
      normalized[[idx]] <- list(param = fit_param, dist = "lognormal", meanlog = meanlog, sdlog = sdlog)
    } else {
      stop(
        "Unsupported prior distribution '", dist, "' for '", param_name,
        "'. Supported distributions are 'normal' and 'lognormal'"
      )
    }
  }

  normalized
}

.normalize_prior_parameter_name <- function(param_name, rtmb_parms, areas) {
  raw_name <- gsub("\\.", "_", param_name)
  raw_name <- sub("^log_B_initial$", "log_B0", raw_name)
  raw_name <- sub("^B_initial$", "B0", raw_name)
  raw_name <- sub("^log_B_initial_", "log_B0_", raw_name)
  raw_name <- sub("^B_initial_", "B0_", raw_name)
  candidates <- raw_name

  if (!grepl("^log_", raw_name)) {
    candidates <- c(candidates, paste0("log_", raw_name))
  }

  if (length(areas) == 1L) {
    area_name <- areas[1]
    if (raw_name %in% c("q", "B0", "log_q", "log_B0")) {
      prefixed <- if (grepl("^log_", raw_name)) raw_name else paste0("log_", raw_name)
      candidates <- c(candidates, paste0(prefixed, "_", area_name))
    }
  }

  matches <- unique(candidates[candidates %in% names(rtmb_parms)])
  if (length(matches) == 0) {
    stop("Prior parameter '", param_name, "' does not match any fitted parameter")
  }

  matches[1]
}

.prepare_environmental_covariates <- function(env_data,
                                              years,
                                              areas,
                                              covariate_names,
                                              lag = 0,
                                              scale_covariates = TRUE) {
  if (length(covariate_names) == 0) {
    return(NULL)
  }
  if (is.null(env_data)) {
    stop("data$env_data is required when options$env_covariates is supplied")
  }

  env_data <- as.data.frame(env_data, stringsAsFactors = FALSE)
  if (!all(c("year", covariate_names) %in% names(env_data))) {
    missing_cols <- setdiff(c("year", covariate_names), names(env_data))
    stop("env_data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!"area" %in% names(env_data)) {
    env_data[["area"]] <- areas[1]
  }

  env_data[["year"]] <- as.integer(env_data[["year"]])
  env_data[["area"]] <- as.character(env_data[["area"]])
  env_data <- env_data[env_data[["year"]] %in% years & env_data[["area"]] %in% areas, , drop = FALSE]
  if (nrow(env_data) == 0) {
    stop("No environmental covariate records overlap model years/areas")
  }

  clean_names <- gsub("[^A-Za-z0-9]", "_", covariate_names)
  clean_names <- make.names(clean_names, unique = TRUE)

  n_years <- length(years)
  n_areas <- length(areas)
  n_cov <- length(covariate_names)
  env_full <- array(NA_real_, dim = c(n_years, n_areas, n_cov), dimnames = list(year = years, area = areas, covariate = clean_names))

  for (ic in seq_len(n_cov)) {
    cov <- covariate_names[ic]
    for (ia in seq_len(n_areas)) {
      area_i <- areas[ia]
      rows <- env_data[["area"]] == area_i
      if (!any(rows)) {
        rows <- rep(TRUE, nrow(env_data))
      }
      idx <- match(years, env_data[["year"]][rows])
      vals <- env_data[[cov]][rows][idx]
      vals <- as.numeric(vals)
      if (all(!is.finite(vals))) {
        stop("Environmental covariate '", cov, "' has no finite values for area '", area_i, "'")
      }
      if (any(!is.finite(vals))) {
        fill <- mean(vals[is.finite(vals)], na.rm = TRUE)
        vals[!is.finite(vals)] <- fill
      }
      env_full[, ia, ic] <- vals
    }
  }

  n_trans <- max(n_years - 1L, 1L)
  env_array <- array(0, dim = c(n_trans, n_areas, n_cov), dimnames = list(transition = seq_len(n_trans), area = areas, covariate = clean_names))
  for (t in seq_len(n_trans)) {
    src_idx <- max(1L, t - lag)
    env_array[t, , ] <- env_full[src_idx, , ]
  }

  scaling <- NULL
  if (scale_covariates) {
    centers <- numeric(n_cov)
    scales <- numeric(n_cov)
    for (ic in seq_len(n_cov)) {
      vals <- as.numeric(env_array[, , ic])
      centers[ic] <- mean(vals)
      scales[ic] <- stats::sd(vals)
      if (!is.finite(scales[ic]) || scales[ic] <= 0) {
        scales[ic] <- 1
      }
      env_array[, , ic] <- (env_array[, , ic] - centers[ic]) / scales[ic]
    }
    scaling <- list(
      names = clean_names,
      center = centers,
      scale = scales,
      lag = lag
    )
  } else {
    scaling <- list(names = clean_names, center = rep(0, n_cov), scale = rep(1, n_cov), lag = lag)
  }

  list(
    env_array = env_array,
    covariate_names = clean_names,
    scaling = scaling
  )
}

#' Preprocess Model Data
#'
#' Aligns and merges CPUE and catch data for model fitting.
#'
#' @param cpue_data Data frame with year, cpue, optional area column
#' @param catch_data Data frame with year, catch, optional area column
#'
#' @return List with aligned data vectors
#'
#' @keywords internal
preprocess_model_data <- function(cpue_data, catch_data) {
  # Coerce to plain data.frame (strip extra attributes/classes that interfere

  # with column access inside the package namespace)
  cpue_data <- as.data.frame(cpue_data, stringsAsFactors = FALSE)
  catch_data <- as.data.frame(catch_data, stringsAsFactors = FALSE)

  # Normalize area column to factor if present
  if (!"area" %in% names(cpue_data)) cpue_data[["area"]] <- factor("A1")
  if (!"area" %in% names(catch_data)) catch_data[["area"]] <- factor("A1")

  cpue_data[["area"]] <- as.factor(cpue_data[["area"]])
  catch_data[["area"]] <- as.factor(catch_data[["area"]])

  # If label column present, treat as multiple indices per area
  has_label <- "label" %in% names(cpue_data)
  if (has_label) cpue_data[["label"]] <- as.factor(cpue_data[["label"]])

  # Ensure data is sorted
  if (has_label) {
    cpue_data <- cpue_data[order(cpue_data[["year"]], cpue_data[["area"]], cpue_data[["label"]]), ]
  } else {
    cpue_data <- cpue_data[order(cpue_data[["year"]], cpue_data[["area"]]), ]
  }
  catch_data <- catch_data[order(catch_data[["year"]], catch_data[["area"]]), ]

  # Model years = complete annual sequence from earliest to latest year in

  # either dataset.  The model needs every year for correct biomass dynamics;
  # CPUE observations may start later or have gaps (stored as NA) and catch

  # for years without data defaults to 0.
  data_years <- sort(union(unique(catch_data[["year"]]), unique(cpue_data[["year"]])))
  all_years <- seq(min(data_years), max(data_years))
  common_areas <- sort(intersect(levels(cpue_data[["area"]]), levels(catch_data[["area"]])))
  labels <- if (has_label) sort(levels(cpue_data[["label"]])) else NULL

  overlap_years <- sort(intersect(unique(cpue_data[["year"]]), unique(catch_data[["year"]])))
  if (length(overlap_years) < 3) stop("Insufficient overlapping years between CPUE and catch data (minimum 3 required)")
  if (length(common_areas) < 1) stop("No overlapping areas between CPUE and catch data")

  # Defensive check: years must be sequential with increment 1
  stopifnot(
    "Model years must be sequential with increment 1" =
      length(all_years) >= 2 && all(diff(all_years) == 1L)
  )

  # Filter to model years and areas
  cpue_rows <- cpue_data[["year"]] %in% all_years & cpue_data[["area"]] %in% common_areas
  if (has_label) cpue_rows <- cpue_rows & cpue_data[["label"]] %in% labels
  cpue_cols <- c("year", "area", "cpue", if (has_label) "label" else NULL)
  cpue_subset <- cpue_data[cpue_rows, cpue_cols]
  catch_subset <- catch_data[catch_data[["year"]] %in% all_years & catch_data[["area"]] %in% common_areas, c("year", "area", "catch")]

  # Create matrices/arrays  (using all_years as index)
  year_index <- match(cpue_subset[["year"]], all_years)
  area_index <- match(cpue_subset[["area"]], common_areas)
  if (has_label) {
    label_index <- match(cpue_subset[["label"]], labels)
    cpue_arr <- array(NA_real_,
      dim = c(length(all_years), length(common_areas), length(labels)),
      dimnames = list(year = all_years, area = common_areas, label = labels)
    )
    for (i in seq_len(nrow(cpue_subset))) {
      cpue_arr[year_index[i], area_index[i], label_index[i]] <- cpue_subset[["cpue"]][i]
    }
  } else {
    cpue_mat <- matrix(NA_real_, nrow = length(all_years), ncol = length(common_areas), dimnames = list(all_years, common_areas))
    cpue_mat[cbind(year_index, area_index)] <- cpue_subset[["cpue"]]
  }

  year_index_c <- match(catch_subset[["year"]], all_years)
  area_index_c <- match(catch_subset[["area"]], common_areas)
  catch_mat <- matrix(0, nrow = length(all_years), ncol = length(common_areas), dimnames = list(all_years, common_areas))
  catch_mat[cbind(year_index_c, area_index_c)] <- catch_subset[["catch"]]

  # Check for missing catch — years with CPUE but no catch get zero catch (reasonable default)
  # Only error if catch is explicitly NA (not just missing from the data)
  if (any(is.na(catch_mat))) {
    warning("Missing catch values set to 0 for years without catch data")
    catch_mat[is.na(catch_mat)] <- 0
  }
  # Allow missing CPUE, likelihood can handle
  if (has_label) {
    if (any(is.na(cpue_arr))) warning("Missing CPUE values detected and will be handled in likelihood")
  } else {
    if (any(is.na(cpue_mat))) warning("Missing CPUE values detected and will be handled in likelihood")
  }

  out <- list(
    years = as.integer(all_years),
    areas = as.character(common_areas),
    catch = if (ncol(catch_mat) == 1) as.numeric(catch_mat[, 1]) else catch_mat
  )
  if (has_label) {
    out$labels <- labels
    out$cpue <- cpue_arr
  } else {
    out$cpue <- if (ncol(cpue_mat) == 1) as.numeric(cpue_mat[, 1]) else cpue_mat
  }
  out
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
  natural_params <- log_params
  nms <- names(log_params)

  if (length(log_params) > 0) {
    log_idx <- grepl("^log_", nms)
    if (any(log_idx)) {
      natural_params[log_idx] <- exp(log_params[log_idx])
      nms[log_idx] <- sub("^log_", "", nms[log_idx])
    }
    if ("theta_rho" %in% nms) {
      rho_idx <- which(nms == "theta_rho")[1]
      natural_params[rho_idx] <- tanh(log_params[rho_idx])
      nms[rho_idx] <- "rho"
    }
  }
  names(natural_params) <- nms

  # Backward-compatibility: if only one area and parameters are suffixed, also provide unsuffixed aliases
  if (!is.null(names(natural_params))) {
    pnames <- names(natural_params)
    # detect single-area suffixed q.* and B0.*
    q_like <- grep("^q\\.[A-Za-z0-9_]+$", pnames, value = TRUE)
    b0_like <- grep("^B0\\.[A-Za-z0-9_]+$", pnames, value = TRUE)
    if (length(q_like) == 1 && !("q" %in% pnames)) {
      natural_params <- c(natural_params, q = unname(natural_params[q_like]))
    }
    if (length(b0_like) == 1 && !("B0" %in% pnames)) {
      natural_params <- c(natural_params, B0 = unname(natural_params[b0_like]))
    }

    # Clarity alias: expose B_initial alongside B0 without breaking backward compatibility
    pnames <- names(natural_params)
    if ("B0" %in% pnames && !("B_initial" %in% pnames)) {
      natural_params <- c(natural_params, B_initial = unname(natural_params["B0"]))
    }
    b0_named <- grep("^B0\\.[A-Za-z0-9_]+$", pnames, value = TRUE)
    if (length(b0_named) > 0) {
      for (nm in b0_named) {
        alias_nm <- sub("^B0\\.", "B_initial.", nm)
        if (!(alias_nm %in% names(natural_params))) {
          natural_params <- c(natural_params, setNames(unname(natural_params[nm]), alias_nm))
        }
      }
    }
  }

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
  multi_area <- is.matrix(data$catch) || is.matrix(data$cpue) || !is.null(data$areas)
  areas <- if (!is.null(data$areas)) data$areas else if (is.matrix(data$cpue)) colnames(data$cpue) else "A1"
  n_areas <- if (multi_area) length(areas) else 1L
  has_labels <- !is.null(data$labels)
  has_env <- !is.null(data$env_array)

  # Extract global parameters
  r <- parameters[["r"]]
  K <- parameters[["K"]]
  m <- parameters[["m"]]

  # Build per-area parameters (q and B0)
  get_param_area <- function(params, base, areas, fallback = NULL) {
    out <- numeric(length(areas))
    names(out) <- areas
    pnames <- names(params)
    for (a in areas) {
      key <- paste0(base, ".", a)
      if (!is.null(pnames) && key %in% pnames) {
        out[[a]] <- unname(params[key])
      } else if (!is.null(fallback) && !is.null(pnames) && fallback %in% pnames) {
        out[[a]] <- unname(params[fallback])
      } else {
        stop("Missing parameter ", key, " and no fallback provided")
      }
    }
    out
  }
  if (multi_area) {
    # Per-area q or per-index q
    if (!has_labels) {
      # Try per-area q first
      q_vec <- try(get_param_area(parameters, "q", areas, fallback = "q"), silent = TRUE)
      if (inherits(q_vec, "try-error") || any(is.na(q_vec))) {
        # fallback to global q
        if (!is.null(names(parameters)) && "q" %in% names(parameters)) {
          q_vec <- rep(parameters[["q"]], length(areas))
          names(q_vec) <- areas
        } else {
          stop("Missing catchability q parameters")
        }
      }
    } else {
      # Build per-index q array [area x label]
      pnames <- names(parameters)
      labels <- data$labels
      q_arr <- array(NA_real_, dim = c(length(areas), length(labels)), dimnames = list(area = areas, label = labels))
      for (a in areas) {
        for (l in labels) {
          key <- paste0("q.", a, ".", l)
          if (!is.null(pnames) && key %in% pnames) {
            q_arr[a, l] <- parameters[[key]]
          } else if (!is.null(pnames) && paste0("q.", a) %in% pnames) {
            q_arr[a, l] <- parameters[[paste0("q.", a)]]
          } else if (!is.null(pnames) && "q" %in% pnames) {
            q_arr[a, l] <- parameters[["q"]]
          } else {
            stop("Missing catchability parameter for area ", a, ", label ", l)
          }
        }
      }
    }
    B0_vec <- get_param_area(parameters, "B0", areas, fallback = "B0")
  } else {
    area1 <- areas[[1]]
    q_key <- paste0("q.", area1)
    b0_key <- paste0("B0.", area1)
    pnames <- names(parameters)
    q_val <- if (!is.null(pnames) && q_key %in% pnames) unname(parameters[q_key]) else unname(parameters["q"])
    b0_val <- if (!is.null(pnames) && b0_key %in% pnames) unname(parameters[b0_key]) else unname(parameters["B0"])
    if (is.null(q_val) || is.null(b0_val)) stop("Missing q or B0 parameter for single-area model")
    q_vec <- c(q_val)
    names(q_vec) <- area1
    B0_vec <- c(b0_val)
    names(B0_vec) <- area1
  }

  # Initialize outputs
  biomass <- matrix(NA_real_, nrow = n_years, ncol = n_areas, dimnames = list(as.character(data$years), areas))
  harvest_rate <- biomass
  fitted_cpue <- biomass

  # Coerce data to matrices (always — even single area returns vectors from preprocess)
  catch_mat <- if (is.matrix(data$catch)) data$catch else matrix(data$catch, ncol = 1, dimnames = list(NULL, areas))
  if (!has_labels) {
    cpue_mat <- if (is.matrix(data$cpue)) data$cpue else matrix(data$cpue, ncol = 1, dimnames = list(NULL, areas))
  } else {
    cpue_arr <- data$cpue
  }

  # Biomass recursion per area
  biomass[1, ] <- as.numeric(B0_vec)

  env_term <- matrix(0, nrow = max(n_years - 1, 1), ncol = n_areas)
  if (has_env) {
    env_array <- data$env_array
    cov_names <- dimnames(env_array)[[3]]
    beta_vec <- numeric(length(cov_names))
    for (ic in seq_along(cov_names)) {
      key_u <- paste0("beta_", cov_names[ic])
      key_d <- paste0("beta.", cov_names[ic])
      if (!is.null(names(parameters)) && key_u %in% names(parameters)) {
        beta_vec[ic] <- parameters[[key_u]]
      } else if (!is.null(names(parameters)) && key_d %in% names(parameters)) {
        beta_vec[ic] <- parameters[[key_d]]
      } else {
        beta_vec[ic] <- 0
      }
    }
    for (t in seq_len(n_years - 1)) {
      for (ia in seq_len(n_areas)) {
        env_term[t, ia] <- sum(env_array[t, ia, ] * beta_vec)
      }
    }
  }

  for (t in 1:(n_years - 1)) {
    Bt <- biomass[t, ]
    production <- ifelse(Bt > 0 & K > 0 & m > 0, r * Bt * (1 - (Bt / K)^(m - 1)) / m, 0)
    production[!is.finite(production)] <- 0
    B_det <- pmax(Bt + production - catch_mat[t, ], 0.01)
    if (has_env) {
      biomass[t + 1, ] <- pmax(exp(log(B_det) + env_term[t, ]), 0.01)
    } else {
      biomass[t + 1, ] <- B_det
    }
  }

  # Apply gravity movement redistribution (if present in data)
  if (!is.null(data$movement_rate) && data$movement_rate > 0) {
    move_rate <- data$movement_rate
    decay_val <- data$decay
    dist_mat <- data$distance_matrix
    attract <- data$attractiveness
    nA <- n_areas
    # Build movement kernel (same as in objective function)
    W <- matrix(0, nA, nA)
    for (a in seq_len(nA)) {
      for (b in seq_len(nA)) {
        W[a, b] <- attract[b] * exp(-decay_val * dist_mat[a, b])
      }
    }
    Kmat <- W / rowSums(W)
    for (t in 2:n_years) {
      Bt <- biomass[t, ]
      biomass[t, ] <- pmax((1 - move_rate) * Bt + move_rate * as.numeric(Kmat %*% Bt), 0.01)
    }
  }

  # Derived quantities
  for (t in 1:n_years) {
    harvest_rate[t, ] <- catch_mat[t, ] / biomass[t, ]
    if (!has_labels) {
      fitted_cpue[t, ] <- as.numeric(q_vec) * biomass[t, ]
    }
  }

  if (has_labels) {
    # Build fitted cpue array [year x area x label]
    labels <- data$labels
    fitted_cpue_arr <- array(NA_real_, dim = c(n_years, n_areas, length(labels)), dimnames = list(year = as.character(data$years), area = areas, label = labels))
    for (t in 1:n_years) {
      for (a in seq_along(areas)) {
        for (l in seq_along(labels)) {
          fitted_cpue_arr[t, a, l] <- q_arr[a, l] * biomass[t, a]
        }
      }
    }
  }

  # Residuals (log-scale), ignore NA cpue
  if (!has_labels) {
    withCallingHandlers(
      {
        residuals <- log(cpue_mat) - log(fitted_cpue)
      },
      warning = function(w) {}
    )
    residuals[!is.finite(residuals)] <- NA
  } else {
    residuals <- log(cpue_arr) - log(fitted_cpue_arr)
    residuals[!is.finite(residuals)] <- NA
  }

  # If single area, drop to vectors for backward compatibility
  if (!multi_area) {
    biomass <- as.numeric(biomass[, 1])
    harvest_rate <- as.numeric(harvest_rate[, 1])
    if (!has_labels) {
      fitted_cpue <- as.numeric(fitted_cpue[, 1])
      residuals <- as.numeric(residuals[, 1])
    } else {
      fitted_cpue <- as.numeric(fitted_cpue_arr[, 1, 1])
      residuals <- as.numeric(residuals[, 1, 1])
    }
  }

  list(
    biomass = biomass,
    harvest_rate = harvest_rate,
    fitted_cpue = if (!has_labels) fitted_cpue else fitted_cpue_arr,
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
  ref <- calculate_reference_points_from_parameters(
    parameters = parameters,
    warn_on_invalid_m = FALSE
  )

  list(
    msy = ref$msy,
    bmsy = ref$bmsy,
    fmsy = ref$fmsy
  )
}
