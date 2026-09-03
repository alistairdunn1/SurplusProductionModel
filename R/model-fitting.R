#' Model Fitting Functions for Pella-Tomlinson Surplus Production Model
#'
#' This file contains the main model fitting infrastructure for the
#' Pella-Tomlinson surplus production model using RTMB for optimization.
#'
#' @name model-fitting
NULL

#' Prepare Model Data
#'
#' Public wrapper around internal model-data preprocessing used by
#' \code{fit_pella_tomlinson_model}. This helper returns the aligned data
#' object expected by low-level fitting utilities.
#'
#' @param cpue_data Data frame with year, cpue, optional area/label columns,
#'   and optional observation-uncertainty inputs such as cv, se, or obs_sd_log.
#' @param catch_data Data frame with year, catch, and optional area column.
#'
#' @return List with aligned model inputs for fitting internals.
#' @export
prepare_model_data <- function(cpue_data, catch_data) {
  preprocess_model_data(cpue_data, catch_data)
}

#' Prepare Starting Values
#'
#' Public wrapper around internal starting-value generation used by
#' \code{fit_pella_tomlinson_model}. Optionally sets all K starts to a
#' common value.
#'
#' @param processed_data List returned by \code{prepare_model_data}.
#' @param k_start Optional positive numeric scalar. If supplied, all
#'   \code{log_K} starts are set to \code{log(k_start)}.
#'
#' @return Named list of starting parameter values on the working scale.
#' @export
prepare_starting_values <- function(processed_data, k_start = NULL) {
  starts <- generate_starting_values(processed_data)

  if (!is.null(k_start)) {
    checkmate::assert_number(
      k_start,
      lower = .Machine$double.eps,
      .var.name = "k_start"
    )
    k_keys <- grep("^log_K(\\.|_|$)", names(starts), value = TRUE)
    if (length(k_keys) == 0) {
      stop("No log_K parameters found in generated starting values.", call. = FALSE)
    }
    starts[k_keys] <- log(k_start)
  }

  starts
}

#' Fit Pella-Tomlinson Surplus Production Model
#'
#' Fits a Pella-Tomlinson surplus production model to catch and CPUE data
#' using maximum likelihood estimation via RTMB.
#'
#' @param data List containing model data with elements:
#'   \describe{
#'     \item{cpue_data}{Data frame with columns: year, cpue, optional: area, label, and observation-uncertainty inputs such as cv, se, or obs_sd_log}
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
#'     \item{estimate_movement_rate}{Logical, estimate \code{movement_rate}
#'       as a model parameter when movement inputs are supplied
#'       (default: FALSE).}
#'     \item{movement_rate_start}{Optional positive numeric starting value
#'       for \code{movement_rate} when \code{estimate_movement_rate = TRUE}.
#'       If NULL, the fitter uses \code{data$movement$movement_rate} when
#'       provided, otherwise 0.12. Must lie in \code{(0, 1)}.}
#'     \item{spinup_years}{Non-negative integer number of deterministic
#'       pre-data years (zero catch) used to resolve the unfished equilibrium
#'       before the first observation year (default: 0, or 50 when movement
#'       inputs are supplied). The spin-up starts from the per-area carrying
#'       capacities and, when movement is present, converges to the joint
#'       (redistributed) unfished spatial equilibrium, which differs from the
#'       per-area \code{K}. Initial depletion \code{d0} is applied as a
#'       multiplier to that equilibrium, so the spin-up resolves the
#'       movement/unfished level while \code{d0} sets depletion; the two are
#'       complementary. Without movement the spin-up is a no-op (each area sits
#'       at its own \code{K}) and \code{d0} is applied directly to \code{K}.}
#'     \item{n_starts}{Integer, number of random restarts (default: 1).
#'       When > 1, the optimizer is run from \code{n_starts} different
#'       starting vectors (the original plus jittered versions) and the
#'       run with the lowest objective is retained.}
#'     \item{jitter_sd}{Numeric, standard deviation of log-normal jitter
#'       applied to starting values for multi-start (default: 0.2).}
#'     \item{area_k_shares}{Deprecated and ignored. Carrying capacity is now
#'       estimated independently for each area in multi-area fits.}
#'     \item{priors}{Optional named list of priors on model parameters.
#'       Each element is a list with \code{dist} and distribution-specific
#'       fields. Supported distributions are:\cr
#'       \code{normal}: \code{list(dist = "normal", mean = ..., sd = ...)}\cr
#'       \code{lognormal}: \code{list(dist = "lognormal", meanlog = ..., sdlog = ...)}\cr
#'       \code{exponential}: \code{list(dist = "exponential", rate = ...)}
#'       (applied on natural scale for log-transformed positive parameters)
#'       Parameter names can be supplied in natural form (e.g. \code{"r"},
#'       \code{"K"}) or log form (e.g. \code{"log_r"}, \code{"log_K"}).}
#'   }
#'
#' @return ProductionModel object with fitted results
#'
#' @details
#' This function implements the complete model fitting workflow:
#' To fit the Fox model, fix \code{m = 1} (or \code{log_m = 0}) through
#' \code{fixed_params}. The RTMB objective then uses the analytic Fox limit.
#'
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
#' The reported \code{aic} and \code{bic} use the optimised objective and the
#' number of estimated fixed parameters. When \code{process_noise = TRUE} the
#' objective is the Laplace-marginal negative log-likelihood, so these are
#' marginal information criteria and carry the usual conditional-AIC caveats
#' for mixed models (\code{results$ic_type} records which form applies).
#'
#' Carrying capacity \code{K}, catchability \code{q}, and initial depletion
#' \code{d0} are jointly only weakly identified from a single relative CPUE
#' index without a strong one-way-trip contrast. Supplying informative priors
#' (see the \code{priors} option) on one or more of these, or fixing \code{d0},
#' is recommended when the index is uninformative about absolute scale.
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
    estimate_movement_rate = FALSE,
    movement_rate_start = NULL,
    spinup_years = 0L,
    area_k_shares = NULL,
    priors = NULL,
    calculate_se = TRUE
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
  estimate_movement_rate <- isTRUE(options$estimate_movement_rate)
  movement_rate_start <- options$movement_rate_start
  spinup_years <- as.integer(options$spinup_years %||% 0L)
  if (!is.finite(spinup_years) || length(spinup_years) != 1 || spinup_years < 0) {
    stop("options$spinup_years must be a single non-negative integer")
  }

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
    if (!is.null(data$movement$transition_matrix)) processed_data$transition_matrix <- data$movement$transition_matrix
    if (!is.null(data$movement$distance_matrix)) processed_data$distance_matrix <- data$movement$distance_matrix
    if (!is.null(data$movement$attractiveness)) processed_data$attractiveness <- data$movement$attractiveness
    if (!is.null(data$movement$movement_rate)) processed_data$movement_rate <- data$movement$movement_rate
    if (!is.null(data$movement$decay)) processed_data$decay <- data$movement$decay
  }

  if (!is.null(processed_data$areas)) {
    area_levels <- as.character(processed_data$areas)

    if (!is.null(processed_data$transition_matrix)) {
      tm <- processed_data$transition_matrix
      if (!is.matrix(tm) || nrow(tm) != length(area_levels) || ncol(tm) != length(area_levels)) {
        stop("movement transition_matrix dimensions must match the number of model areas")
      }
      if (!is.null(rownames(tm)) && !is.null(colnames(tm))) {
        if (!all(area_levels %in% rownames(tm)) || !all(area_levels %in% colnames(tm))) {
          stop("movement transition_matrix row/column names must include all model areas")
        }
        tm <- tm[area_levels, area_levels, drop = FALSE]
      }
      if (any(!is.finite(tm)) || any(tm < 0) || any(abs(rowSums(tm) - 1) > 1e-10)) {
        stop("movement transition_matrix must be finite, non-negative, and have rows that sum to one")
      }
      processed_data$transition_matrix <- tm
    }

    if (!is.null(processed_data$distance_matrix)) {
      dm <- processed_data$distance_matrix
      if (!is.matrix(dm)) {
        stop("movement distance_matrix must be a matrix")
      }
      if (nrow(dm) != length(area_levels) || ncol(dm) != length(area_levels)) {
        stop("movement distance_matrix dimensions must match the number of model areas")
      }
      dm_rn <- rownames(dm)
      dm_cn <- colnames(dm)
      if (!is.null(dm_rn) && !is.null(dm_cn)) {
        if (!all(area_levels %in% dm_rn) || !all(area_levels %in% dm_cn)) {
          stop("movement distance_matrix row/column names must include all model areas")
        }
        processed_data$distance_matrix <- dm[area_levels, area_levels, drop = FALSE]
      }
    }

    if (!is.null(processed_data$attractiveness)) {
      att <- processed_data$attractiveness
      if (!is.numeric(att) || length(att) != length(area_levels)) {
        stop("movement attractiveness must be a numeric vector with one value per model area")
      }
      if (!is.null(names(att))) {
        if (!all(area_levels %in% names(att))) {
          stop("movement attractiveness names must include all model areas")
        }
        processed_data$attractiveness <- as.numeric(att[area_levels])
      } else {
        processed_data$attractiveness <- as.numeric(att)
      }
    }
  }

  has_transition_matrix <- !is.null(processed_data$transition_matrix)
  has_movement_inputs <- has_transition_matrix || (!is.null(processed_data$distance_matrix) &&
    !is.null(processed_data$attractiveness) &&
    !is.null(processed_data$decay))

  # The spin-up resolves the unfished spatial equilibrium, which under movement
  # differs from the per-area carrying capacities (biomass is redistributed
  # between areas). Initial depletion d0 is applied to that equilibrium, so the
  # spin-up and d0 are complementary rather than in conflict. When movement is
  # present but no spin-up was requested, enable a default spin-up so the
  # equilibrium (and hence status baseline) is well defined.
  if (has_movement_inputs && spinup_years == 0L) {
    spinup_years <- 50L
    if (isTRUE(options$show_starting_values_message)) {
      message(
        "Movement inputs supplied: using a 50-year zero-catch spin-up to ",
        "resolve the unfished spatial equilibrium (initial depletion d0 is ",
        "applied to that equilibrium)."
      )
    }
  }

  # Retain the resolved option for the plain-R reconstruction performed after
  # optimisation. The RTMB objective and reported biomass must use the same
  # pre-fishing spatial equilibrium.
  processed_data$spinup_years <- spinup_years

  if (estimate_movement_rate && !has_movement_inputs) {
    stop(
      "options$estimate_movement_rate = TRUE requires data$movement with distance_matrix, attractiveness, and decay"
    )
  }
  if (estimate_movement_rate && has_transition_matrix) {
    stop("estimate_movement_rate cannot be used with a complete transition_matrix")
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
    required_params <- c("log_r", "log_m", "log_sigma_proc", "log_sigma_obs")
    # Per-area and per-index params when multi-area / multi-index
    if (!is.null(processed_data$areas)) {
      area_suffix <- paste0(".", processed_data$areas)
      k_names <- paste0("log_K", area_suffix)
      if (!is.null(processed_data$labels)) {
        # q per area and label
        label_suffix <- paste0(".", processed_data$labels)
        q_names <- as.vector(outer(paste0("log_q", area_suffix), label_suffix, paste0))
        required_params <- c(required_params, k_names, q_names, "log_d0")
      } else {
        required_params <- c(required_params, k_names, paste0("log_q", area_suffix), "log_d0")
      }
    } else {
      required_params <- c(required_params, "log_K", "log_q", "log_d0")
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
  # log_q/log_K -> log_q_<area>/log_K_<area> (objective always uses area suffix when areas are present)
  if (n_areas == 1) {
    a <- areas[1]
    q_key <- paste0("log_q_", a)
    if ("log_q" %in% names(params_init) && !q_key %in% names(params_init)) {
      params_init[[q_key]] <- params_init[["log_q"]]
      params_init[["log_q"]] <- NULL
    }
    k_key <- paste0("log_K_", a)
    if ("log_K" %in% names(params_init) && !k_key %in% names(params_init)) {
      params_init[[k_key]] <- params_init[["log_K"]]
      params_init[["log_K"]] <- NULL
    }
  }

  if (!"log_d0" %in% names(params_init)) {
    params_init[["log_d0"]] <- log(0.8)
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
  if (estimate_movement_rate && !"log_movement_rate" %in% names(params_init)) {
    move_start <- movement_rate_start %||% processed_data$movement_rate %||% 0.12
    if (!is.numeric(move_start) || length(move_start) != 1 || !is.finite(move_start) || move_start <= 0 || move_start >= 1) {
      stop("movement_rate_start must be a single finite numeric value in (0, 1)")
    }
    params_init[["log_movement_rate"]] <- stats::qlogis(as.numeric(move_start))
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

  cpue_sd <- processed_data$cpue_sd %||% NULL
  if (!has_labels && !is.null(cpue_sd) && !is.matrix(cpue_sd)) {
    cpue_sd <- matrix(cpue_sd, ncol = 1, dimnames = list(NULL, areas))
  }

  rtmb_data <- list(
    n_years = n_years,
    n_areas = n_areas,
    areas = areas,
    catch_mat = catch_mat,
    cpue_obs = cpue_obs,
    cpue_sd = cpue_sd,
    labels = if (has_labels) processed_data$labels else NULL,
    process_error_structure = process_error_structure,
    spinup_years = spinup_years
  )

  if (!is.null(processed_data$env_array)) {
    rtmb_data$env_array <- processed_data$env_array
    rtmb_data$env_covariate_names <- processed_data$env_covariate_names
  }

  # Attach movement data
  if (has_movement_inputs) {
    if (has_transition_matrix) {
      rtmb_data$transition_matrix <- processed_data$transition_matrix
    }
    if (!is.null(processed_data$movement_rate)) {
      rtmb_data$movement_rate <- as.numeric(processed_data$movement_rate)
    }
    if (!has_transition_matrix) {
      rtmb_data$decay <- processed_data$decay
      rtmb_data$distance_matrix <- processed_data$distance_matrix
      rtmb_data$attractiveness <- processed_data$attractiveness
    }
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

  # Process-error variance is absent from the deterministic objective. Keep
  # its supplied value for result compatibility, but do not optimise an
  # unidentifiable parameter.
  if (!use_process_noise && "log_sigma_proc" %in% names(rtmb_parms)) {
    map_list$log_sigma_proc <- factor(NA)
  }

  if (!is.null(fixed) && length(fixed) > 0) {
    for (orig in names(fixed)) {
      # RTMB parameter names use underscores (log_K_A1); fixed_params may be
      # supplied with dot-separated area/label suffixes (log_K.A1).
      fn <- gsub("\\.", "_", orig)
      # Single-area fits store per-area names, so a bare "log_K"/"log_q" maps
      # to its area-suffixed counterpart when there is exactly one area.
      if (!fn %in% names(rtmb_parms) && n_areas == 1L) {
        suffixed <- paste0(fn, "_", areas[1])
        if (suffixed %in% names(rtmb_parms)) fn <- suffixed
      }
      if (!fn %in% names(rtmb_parms)) {
        stop("fixed_params name '", orig, "' not found in parameter list")
      }
      rtmb_parms[[fn]] <- unname(fixed[[orig]])
      map_list[[fn]] <- factor(NA)
    }
  }
  # RTMB cannot branch on the AD-valued shape parameter.  If the Fox shape is
  # fixed exactly (m = 1; log_m = 0), pass a data-level flag so the objective
  # evaluates the analytic Fox production equation.
  rtmb_data$fox_mode <- "log_m" %in% names(map_list) &&
    isTRUE(all.equal(as.numeric(rtmb_parms$log_m), 0))

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

  # Retain a data-derived alternative to a supplied sequential warm start.
  # A management trajectory can move into a different likelihood basin, so
  # jittering only the preceding fit does not provide independent starts.
  automatic_start <- rename_dots(generate_starting_values(processed_data))
  if (n_areas == 1L) {
    area <- areas[[1]]
    q_name <- paste0("log_q_", area)
    k_name <- paste0("log_K_", area)
    if ("log_q" %in% names(automatic_start)) {
      automatic_start[[q_name]] <- automatic_start[["log_q"]]
    }
    if ("log_K" %in% names(automatic_start)) {
      automatic_start[[k_name]] <- automatic_start[["log_K"]]
    }
  }
  automatic_start_par <- obj$par
  common_start_names <- intersect(
    names(automatic_start_par),
    names(automatic_start)
  )
  automatic_start_par[common_start_names] <- unlist(
    automatic_start[common_start_names],
    use.names = FALSE
  )

  lower_bounds <- rep(-Inf, length(obj$par))
  upper_bounds <- rep(Inf, length(obj$par))
  names(lower_bounds) <- names(upper_bounds) <- names(obj$par)
  if ("log_d0" %in% names(obj$par)) {
    # Initial depletion is a proportion of unfished biomass.
    lower_bounds[["log_d0"]] <- log(1e-4)
    upper_bounds[["log_d0"]] <- 0
  }
  clamp_to_bounds <- function(par) {
    pmin(pmax(par, lower_bounds), upper_bounds)
  }
  automatic_start_par <- clamp_to_bounds(automatic_start_par)

  # ---- Optimise (with optional multi-start) ----------------------------
  n_starts <- max(1L, as.integer(options$n_starts %||% 1L))
  jitter_sd <- as.numeric(options$jitter_sd %||% 0.2)

  start_time <- Sys.time()

  best_opt <- NULL
  best_nll <- Inf
  best_converged <- FALSE

  for (istart in seq_len(n_starts)) {
    if (istart == 1L) {
      start_par <- obj$par
    } else if (istart == 2L && !is.null(params_init)) {
      start_par <- automatic_start_par
    } else {
      # Use deterministic dispersed offsets on the log scale. Optimisation
      # must not consume the simulation random-number stream because that
      # changes future observations and parallel reproducibility.
      start_centre <- if (istart %% 2L == 0L) {
        automatic_start_par
      } else {
        obj$par
      }
      coordinate <- seq_along(obj$par)
      phase <- (
        coordinate * 0.618033988749895 +
          (istart - 2L) * 0.414213562373095
      ) %% 1
      phase <- pmin(pmax(phase, 1e-6), 1 - 1e-6)
      start_par <- start_centre + stats::qnorm(phase) * jitter_sd
    }
    this_opt <- tryCatch(
      nlminb(
        start     = clamp_to_bounds(start_par),
        objective = obj$fn,
        gradient  = obj$gr,
        lower     = lower_bounds,
        upper     = upper_bounds,
        control   = options$control
      ),
      error = function(e) NULL
    )
    if (!is.null(this_opt) && is.finite(this_opt$objective)) {
      this_converged <- isTRUE(this_opt$convergence == 0)
      replace_best <- is.null(best_opt) ||
        (this_converged && !best_converged) ||
        (this_converged == best_converged &&
          this_opt$objective < best_nll)

      if (replace_best) {
        best_nll <- this_opt$objective
        best_opt <- this_opt
        best_converged <- this_converged
      }
    }
  }

  # nlminb can report false convergence or an iteration limit after reaching
  # a useful parameter region. Restart once from the best finite candidate so
  # that convergence is assessed from that region rather than from another
  # jittered initial value. A non-converged restart remains a failed fit.
  if (!is.null(best_opt) && !best_converged) {
    restart_opt <- tryCatch(
      nlminb(
        start = clamp_to_bounds(best_opt$par),
        objective = obj$fn,
        gradient = obj$gr,
        lower = lower_bounds,
        upper = upper_bounds,
        control = options$control
      ),
      error = function(e) NULL
    )
    if (!is.null(restart_opt) && is.finite(restart_opt$objective)) {
      restart_converged <- isTRUE(restart_opt$convergence == 0)
      if (restart_converged || restart_opt$objective < best_nll) {
        best_opt <- restart_opt
        best_nll <- restart_opt$objective
        best_converged <- restart_converged
      }
    }
  }

  # Use an independent quasi-Newton convergence check when PORT cannot certify
  # the best finite solution. This addresses false-convergence diagnostics
  # without accepting a non-converged estimate or changing the objective.
  if (!is.null(best_opt) && !best_converged) {
    bfgs_control <- list(
      maxit = as.integer(options$control$iter.max %||% 500L),
      factr = 1e7,
      pgtol = as.numeric(options$control$rel.tol %||% 1e-8)
    )
    bfgs_opt <- tryCatch(
      stats::optim(
        par = best_opt$par,
        fn = obj$fn,
        gr = obj$gr,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        control = bfgs_control
      ),
      error = function(e) NULL
    )
    if (!is.null(bfgs_opt) &&
      is.finite(bfgs_opt$value) &&
      isTRUE(bfgs_opt$convergence == 0)) {
      best_opt <- list(
        par = bfgs_opt$par,
        objective = bfgs_opt$value,
        convergence = 0L,
        iterations = unname(bfgs_opt$counts[["function"]]),
        evaluations = unname(bfgs_opt$counts),
        message = "relative convergence (L-BFGS-B verification)"
      )
      best_nll <- bfgs_opt$value
      best_converged <- TRUE
    }
  }

  # Fall back to the first run if all starts failed
  if (is.null(best_opt)) {
    best_opt <- nlminb(
      start     = clamp_to_bounds(obj$par),
      objective = obj$fn,
      gradient  = obj$gr,
      lower     = lower_bounds,
      upper     = upper_bounds,
      control   = options$control
    )
  }
  opt_result <- best_opt

  end_time <- Sys.time()
  fitting_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Check convergence
  if (opt_result$convergence != 0) {
    gradient_at_solution <- tryCatch(
      obj$gr(opt_result$par),
      error = function(e) rep(NA_real_, length(opt_result$par))
    )
    warning(
      "Model did not converge. Convergence code: ", opt_result$convergence,
      ". Message: ", opt_result$message,
      ". Objective: ", signif(opt_result$objective, 8),
      ". Maximum absolute gradient: ",
      signif(max(abs(gradient_at_solution), na.rm = TRUE), 8),
      ". Parameters: ", paste(
        paste0(names(opt_result$par), "=", signif(opt_result$par, 6)),
        collapse = ", "
      )
    )
  }

  # ---- sdreport for standard errors -----------------------------------
  std_errors <- NULL
  hessian_valid <- FALSE
  sdr <- NULL

  if (isTRUE(options$calculate_se)) {
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
  }

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
  scalar_par_names <- grep("^(log_r|log_K|log_m|log_sigma|log_q|log_d0)", names(full_log_par), value = TRUE)
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

  # Recover the conditional modes of the fitted process deviations. When
  # standard errors were calculated, sdreport provides the modes used for its
  # ADREPORT trajectory. Otherwise, update the inner Laplace solution at the
  # optimum and restore the original parameter dimensions with parList().
  process_deviations <- NULL
  if (use_process_noise) {
    expected_dimensions <- c(n_years - 1L, n_areas)
    if (!is.null(sdr) &&
        length(sdr$par.random) == prod(expected_dimensions)) {
      process_deviations <- matrix(
        as.numeric(sdr$par.random),
        nrow = expected_dimensions[[1]],
        ncol = expected_dimensions[[2]]
      )
    } else {
      invisible(obj$fn(opt_result$par))
      conditional_parameters <- obj$env$parList()
      process_deviations <- conditional_parameters$proc_dev
    }
    if (!is.matrix(process_deviations) ||
        !identical(dim(process_deviations), expected_dimensions) ||
        any(!is.finite(process_deviations))) {
      stop(
        "RTMB did not return a finite process-deviation matrix with dimensions ",
        paste(expected_dimensions, collapse = " x "),
        call. = FALSE
      )
    }
    dimnames(process_deviations) <- list(
      year = as.character(processed_data$years[-n_years]),
      area = areas
    )
  }

  # Recompute derived quantities in plain R from fitted parameters.
  # This avoids RTMB REPORT/advector complications (NaN promotion,
  # dimension loss) while still using ADREPORT for sdreport SEs.
  model_results <- calculate_model_results(
    fitted_params,
    processed_data,
    process_deviations = process_deviations
  )
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
    # Retain the complete optimised parameter vector on the fitting scale.
    # This includes mapped parameters and permits sequential assessments to
    # start from the preceding converged solution.
    log_parameters = log_par_scalar,
    std_errors = std_errors,
    biomass = B_est,
    biomass_se = biomass_se,
    # Unfished spatial equilibrium and total unfished biomass (B0). Under
    # movement B0 differs from the sum of per-area carrying capacities and is
    # the appropriate baseline for depletion/status reporting.
    b_unfished = model_results$b_unfished,
    b0_total = model_results$b0_total,
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
    # With process deviations the objective is the Laplace-marginal negative
    # log-likelihood and n_est_pars counts only the fixed effects, so these are
    # marginal AIC/BIC. They are not directly comparable with the conditional
    # AIC of a mixed model and understate the flexibility contributed by the
    # random effects; interpret with the usual conditional-AIC caveats.
    aic = 2 * n_est_pars + 2 * opt_result$objective,
    bic = log(n_obs) * n_est_pars + 2 * opt_result$objective,
    ic_type = if (use_process_noise) "marginal (Laplace); see conditional-AIC caveats" else "standard",
    process_noise = use_process_noise,
    process_error_structure = if (use_process_noise) process_error_structure else "none",
    process_deviations = process_deviations,
    rho = if ("rho" %in% names(fitted_params)) fitted_params[["rho"]] else NA_real_,
    movement_rate = if ("movement_rate" %in% names(fitted_params)) fitted_params[["movement_rate"]] else (processed_data$movement_rate %||% NA_real_),
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
        cpue_sd = processed_data$cpue_sd %||% NULL,
        effort = effort,
        spinup_years = spinup_years,
        movement = if (has_movement_inputs) {
          if (has_transition_matrix) list(
            transition_matrix = processed_data$transition_matrix
          ) else list(
            distance_matrix = processed_data$distance_matrix,
            attractiveness = processed_data$attractiveness,
            decay = processed_data$decay,
            movement_rate = if ("movement_rate" %in% names(fitted_params)) {
              as.numeric(fitted_params[["movement_rate"]])
            } else {
              as.numeric(processed_data$movement_rate %||% NA_real_)
            }
          )
        } else {
          NULL
        },
        env_scaling = processed_data$env_scaling %||% NULL
      ),
      results = results_list,
      fitted = isTRUE(opt_result$convergence == 0),
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

    fit_params <- .normalize_prior_parameter_name(param_name, rtmb_parms, areas)
    dist <- spec$dist %||% spec$distribution
    if (is.null(dist) || !is.character(dist) || length(dist) != 1) {
      stop("Prior for '", param_name, "' must supply 'dist' or 'distribution'")
    }
    dist <- tolower(dist)

    for (fit_param in fit_params) {
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
      } else if (dist == "beta") {
        shape1 <- spec$shape1 %||% spec$alpha
        shape2 <- spec$shape2 %||% spec$beta
        if (!grepl("^log_", fit_param)) {
          stop("Beta priors are only supported for log-scale parameters constrained to (0, 1); use parameter '", param_name, "' on log scale")
        }
        if (!is.numeric(shape1) || length(shape1) != 1 || !is.finite(shape1) || shape1 <= 0) {
          stop("Beta prior for '", param_name, "' must supply a positive finite scalar 'shape1' (or 'alpha')")
        }
        if (!is.numeric(shape2) || length(shape2) != 1 || !is.finite(shape2) || shape2 <= 0) {
          stop("Beta prior for '", param_name, "' must supply a positive finite scalar 'shape2' (or 'beta')")
        }
        normalized[[idx]] <- list(param = fit_param, dist = "beta", shape1 = shape1, shape2 = shape2)
      } else if (dist %in% c("exponential", "exp")) {
        rate <- spec$rate %||% spec$lambda
        if (!grepl("^log_", fit_param)) {
          stop("Exponential priors are only supported for log-scale positive parameters; use parameter '", param_name, "' on log scale")
        }
        if (!is.numeric(rate) || length(rate) != 1 || !is.finite(rate) || rate <= 0) {
          stop("Exponential prior for '", param_name, "' must supply a positive finite scalar 'rate' (or 'lambda')")
        }
        normalized[[idx]] <- list(param = fit_param, dist = "exponential", rate = rate)
      } else {
        stop(
          "Unsupported prior distribution '", dist, "' for '", param_name,
          "'. Supported distributions are 'normal', 'lognormal', 'beta', and 'exponential'"
        )
      }
    }
  }

  normalized
}

.normalize_prior_parameter_name <- function(param_name, rtmb_parms, areas) {
  raw_name <- gsub("\\.", "_", param_name)
  candidates <- raw_name

  if (!grepl("^log_", raw_name)) {
    candidates <- c(candidates, paste0("log_", raw_name))
  }

  if (raw_name %in% c("K", "log_K") && length(areas) >= 1L) {
    prefixed <- if (grepl("^log_", raw_name)) raw_name else paste0("log_", raw_name)
    candidates <- c(candidates, paste0(prefixed, "_", areas))
  }

  if (length(areas) == 1L) {
    area_name <- areas[1]
    if (raw_name %in% c("q", "d0", "log_q", "log_d0")) {
      prefixed <- if (grepl("^log_", raw_name)) raw_name else paste0("log_", raw_name)
      candidates <- c(candidates, paste0(prefixed, "_", area_name))
    }
  }

  matches <- unique(candidates[candidates %in% names(rtmb_parms)])
  if (length(matches) == 0) {
    stop("Prior parameter '", param_name, "' does not match any fitted parameter")
  }

  matches
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
#' @param cpue_data Data frame with year, cpue, optional area/label columns, and optional observation-uncertainty inputs such as cv, se, or obs_sd_log
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

  derive_cpue_sd_log <- function(df) {
    n <- nrow(df)
    out <- rep(NA_real_, n)

    direct_cols <- c("obs_sd_log", "cpue_sd_log", "se_log", "sigma_log")
    direct_col <- direct_cols[direct_cols %in% names(df)][1]
    if (!is.na(direct_col)) {
      out <- as.numeric(df[[direct_col]])
      return(out)
    }

    cv_col <- c("cv", "CV")
    cv_col <- cv_col[cv_col %in% names(df)][1]
    if (!is.na(cv_col)) {
      cv <- as.numeric(df[[cv_col]])
      valid <- is.finite(cv) & cv >= 0
      out[valid] <- sqrt(log1p(cv[valid]^2))
      return(out)
    }

    se_col <- c("se", "SE", "obs_se", "cpue_se")
    se_col <- se_col[se_col %in% names(df)][1]
    if (!is.na(se_col)) {
      se <- as.numeric(df[[se_col]])
      cpue <- as.numeric(df[["cpue"]])
      valid <- is.finite(se) & se >= 0 & is.finite(cpue) & cpue > 0
      out[valid] <- sqrt(log1p((se[valid] / cpue[valid])^2))
    }

    out
  }

  cpue_data[["obs_sd_log"]] <- derive_cpue_sd_log(cpue_data)

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
  cpue_cols <- c("year", "area", "cpue", if (has_label) "label" else NULL, "obs_sd_log")
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
    cpue_sd_arr <- array(NA_real_,
      dim = c(length(all_years), length(common_areas), length(labels)),
      dimnames = list(year = all_years, area = common_areas, label = labels)
    )
    for (i in seq_len(nrow(cpue_subset))) {
      cpue_arr[year_index[i], area_index[i], label_index[i]] <- cpue_subset[["cpue"]][i]
      cpue_sd_arr[year_index[i], area_index[i], label_index[i]] <- cpue_subset[["obs_sd_log"]][i]
    }
  } else {
    cpue_mat <- matrix(NA_real_, nrow = length(all_years), ncol = length(common_areas), dimnames = list(all_years, common_areas))
    cpue_sd_mat <- matrix(NA_real_, nrow = length(all_years), ncol = length(common_areas), dimnames = list(all_years, common_areas))
    cpue_mat[cbind(year_index, area_index)] <- cpue_subset[["cpue"]]
    cpue_sd_mat[cbind(year_index, area_index)] <- cpue_subset[["obs_sd_log"]]
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
    if (any(is.finite(cpue_sd_arr))) out$cpue_sd <- cpue_sd_arr
  } else {
    out$cpue <- if (ncol(cpue_mat) == 1) as.numeric(cpue_mat[, 1]) else cpue_mat
    if (any(is.finite(cpue_sd_mat))) {
      out$cpue_sd <- if (ncol(cpue_sd_mat) == 1) as.numeric(cpue_sd_mat[, 1]) else cpue_sd_mat
    }
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
    move_idx <- nms == "log_movement_rate"
    if (any(move_idx)) {
      natural_params[move_idx] <- stats::plogis(log_params[move_idx])
      nms[move_idx] <- "movement_rate"
    }
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

  return(natural_params)
}

#' Calculate Model Results
#'
#' Calculate biomass trajectory, harvest rates, and fitted values.
#'
#' @param parameters Named vector of fitted parameters on natural scale
#' @param data List of model data
#' @param process_deviations Optional matrix of fitted log-scale process
#'   deviations, with transition years in rows and areas in columns
#'
#' @return List with model results
#'
#' @keywords internal
calculate_model_results <- function(parameters, data, process_deviations = NULL) {
  n_years <- length(data$years)
  spinup_years <- as.integer(data$spinup_years %||% 0L)
  multi_area <- is.matrix(data$catch) || is.matrix(data$cpue) || !is.null(data$areas)
  areas <- if (!is.null(data$areas)) data$areas else if (is.matrix(data$cpue)) colnames(data$cpue) else "A1"
  n_areas <- if (multi_area) length(areas) else 1L
  has_labels <- !is.null(data$labels)
  has_env <- !is.null(data$env_array)
  has_process_deviations <- !is.null(process_deviations)

  if (has_process_deviations) {
    expected_dimensions <- c(n_years - 1L, n_areas)
    if (!is.matrix(process_deviations) ||
        !identical(dim(process_deviations), expected_dimensions) ||
        any(!is.finite(process_deviations))) {
      stop(
        "process_deviations must be a finite matrix with dimensions ",
        paste(expected_dimensions, collapse = " x "),
        call. = FALSE
      )
    }
  }

  # Extract global parameters
  r <- parameters[["r"]]
  m <- parameters[["m"]]

  # Build per-area parameters (q)
  get_param_area <- function(params, base, areas) {
    out <- numeric(length(areas))
    names(out) <- areas
    pnames <- names(params)
    for (a in areas) {
      key <- paste0(base, ".", a)
      if (!is.null(pnames) && key %in% pnames) {
        out[[a]] <- unname(params[key])
      } else {
        stop("Missing required area-specific parameter ", key)
      }
    }
    out
  }
  if (multi_area) {
    # Per-area q or per-index q
    if (!has_labels) {
      q_vec <- get_param_area(parameters, "q", areas)
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
          } else {
            stop("Missing catchability parameter for area ", a, ", label ", l)
          }
        }
      }
    }
  } else {
    area1 <- areas[[1]]
    q_key <- paste0("q.", area1)
    pnames <- names(parameters)
    q_val <- if (!is.null(pnames) && q_key %in% pnames) unname(parameters[q_key]) else unname(parameters["q"])
    if (is.null(q_val)) stop("Missing q parameter for single-area model")
    q_vec <- c(q_val)
    names(q_vec) <- area1
  }

  # Initial depletion and area-specific carrying capacities
  pnames <- names(parameters)
  if (!is.null(pnames) && "d0" %in% pnames) {
    d0 <- as.numeric(parameters[["d0"]])
  } else {
    stop("Missing d0 parameter")
  }
  if (!is.finite(d0) || d0 <= 0) {
    stop("Parameter d0 must be positive and finite")
  }

  K_vec <- if (multi_area) {
    get_param_area(parameters, "K", areas)
  } else {
    area1 <- areas[[1]]
    pnames <- names(parameters)
    if (!is.null(pnames) && "K" %in% pnames) {
      out <- c(unname(parameters[["K"]]))
      names(out) <- area1
      out
    } else {
      stop("Missing K parameter for single-area model")
    }
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

  has_transition_matrix <- !is.null(data$transition_matrix)
  move_rate <- if (has_transition_matrix) {
    1
  } else if (!is.null(names(parameters)) && "movement_rate" %in% names(parameters)) {
    as.numeric(parameters[["movement_rate"]])
  } else {
    as.numeric(data$movement_rate %||% NA_real_)
  }
  has_movement <- has_transition_matrix || (is.finite(move_rate) && move_rate > 0 &&
    !is.null(data$distance_matrix) && !is.null(data$attractiveness) && !is.null(data$decay))
  if (has_transition_matrix) {
    Kmat <- data$transition_matrix
  } else if (has_movement) {
    decay_val <- data$decay
    dist_mat <- data$distance_matrix
    attract <- data$attractiveness
    if (!is.null(rownames(dist_mat)) && !is.null(colnames(dist_mat))) {
      if (!all(areas %in% rownames(dist_mat)) || !all(areas %in% colnames(dist_mat))) {
        stop("movement distance_matrix row/column names must include all model areas")
      }
      dist_mat <- dist_mat[areas, areas, drop = FALSE]
    }
    if (!is.null(names(attract))) {
      if (!all(areas %in% names(attract))) {
        stop("movement attractiveness names must include all model areas")
      }
      attract <- as.numeric(attract[areas])
    } else {
      attract <- as.numeric(attract)
    }
    names(attract) <- areas
    nA <- n_areas
    W <- matrix(0, nA, nA)
    rownames(W) <- areas
    colnames(W) <- areas
    for (a in seq_len(nA)) {
      for (b in seq_len(nA)) {
        W[a, b] <- attract[b] * exp(-decay_val * dist_mat[a, b])
      }
    }
    Kmat <- W / rowSums(W)
  }

  # Biomass recursion per area.
  # Resolve the unfished spatial equilibrium, then apply initial depletion.
  # The spin-up starts from the per-area carrying capacities and, when movement
  # is present, converges to the joint (redistributed) unfished equilibrium
  # rather than the per-area K. Initial depletion d0 scales that equilibrium, so
  # the spin-up resolves the movement/unfished level while d0 sets depletion.
  b_unfished <- as.numeric(K_vec)
  if (spinup_years > 0) {
    b_spin <- as.numeric(K_vec)
    for (s in seq_len(spinup_years)) {
      production <- .pt_production(b_spin, r, K_vec, m)
      production[!is.finite(production)] <- 0
      b_next <- pmax(b_spin + production, 0.01)
      if (has_movement) {
        b_next <- pmax((1 - move_rate) * b_next + move_rate * as.numeric(t(Kmat) %*% b_next), 0.01)
      }
      b_spin <- b_next
    }
    b_unfished <- b_spin
  }
  names(b_unfished) <- areas
  biomass[1, ] <- d0 * b_unfished

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
    production <- .pt_production(Bt, r, K_vec, m)
    production[!is.finite(production)] <- 0
    B_det <- pmax(Bt + production - catch_mat[t, ], 0.01)
    b_next <- if (has_env || has_process_deviations) {
      log_increment <- env_term[t, ]
      if (has_process_deviations) {
        log_increment <- log_increment + process_deviations[t, ]
      }
      pmax(exp(log(B_det) + log_increment), 0.01)
    } else {
      B_det
    }
    # Apply gravity movement within the annual step so the redistributed
    # biomass carries forward as the next year's starting state. This couples
    # areas over time and matches the RTMB objective (a genuine movement
    # process rather than a post-hoc reshaping of the reported trajectory).
    if (has_movement) {
      b_next <- pmax(
        (1 - move_rate) * b_next + move_rate * as.numeric(t(Kmat) %*% b_next),
        0.01
      )
    }
    biomass[t + 1, ] <- b_next
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
    residuals = residuals,
    # Unfished spatial equilibrium (per area) and its total. Under movement this
    # is the joint redistributed equilibrium and is the appropriate unfished
    # baseline (B0) for depletion and status, distinct from the sum of per-area
    # carrying capacities.
    b_unfished = b_unfished,
    b0_total = sum(b_unfished)
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
