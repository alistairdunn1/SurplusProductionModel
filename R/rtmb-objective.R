#' RTMB Objective Function for Pella-Tomlinson Model
#'
#' This file contains the RTMB objective function implementation for the
#' Pella-Tomlinson surplus production model with state-space formulation.
#'
#' @name rtmb-objective
NULL

#' Create RTMB Objective Function for Pella-Tomlinson Model
#'
#' Creates an RTMB objective function for parameter estimation using maximum
#' likelihood estimation with automatic differentiation.
#'
#' @param data List containing model data with elements:
#'   \describe{
#'     \item{years}{Integer vector of years}
#'     \item{catch}{Numeric vector of catch data (tonnes)}
#'     \item{cpue}{Numeric vector of CPUE index data}
#'   }
#' @param parameters Named list of starting parameter values:
#'   \describe{
#'     \item{log_r}{Log intrinsic growth rate}
#'     \item{log_K}{Log carrying capacity}
#'     \item{log_m}{Log Pella-Tomlinson shape parameter}
#'     \item{log_q}{Log catchability coefficient}
#'     \item{log_sigma_proc}{Log process error standard deviation}
#'     \item{log_sigma_obs}{Log observation error standard deviation}
#'     \item{log_B0}{Log initial biomass}
#'   }
#'
#' @return RTMB objective function object
#'
#' @details
#' The objective function implements the negative log-likelihood for the
#' Pella-Tomlinson surplus production model in state-space form:
#'
#' Production function: \deqn{P(B) = r .B . (1 - (B/K)^(m-1)) / m}
#'
#' State equation: \deqn{B[t+1] = B[t] + P(B[t]) - C[t] + ε[t]}
#' where ε[t] ~ N(0, σ²_process)
#'
#' Observation equation: \deqn{CPUE[t] = q . B[t] . exp(η[t])}
#' where η[t] ~ N(0, σ²_obs)
#'
#' Parameters are log-transformed to ensure positivity constraints.
#'
#' @keywords internal
#' Create Simple Objective Function for Pella-Tomlinson Model
#'
#' Creates a simple objective function for parameter estimation using maximum
#' likelihood estimation without automatic differentiation (MVP version).
#'
#' @param data List containing model data with elements 'years', 'catch', 'cpue'
#' @param initial_params Named list of starting parameter values (log-transformed)
#'
#' @return Function that takes parameter vector and returns negative log-likelihood
#'
#' @import RTMB
#' @keywords internal
create_simple_objective <- function(data, initial_params) {
  # Store parameter names for reference
  param_names <- names(initial_params)

  # Return objective function
  function(par) {
    # Convert parameter vector back to named list
    if (is.null(names(par))) {
      names(par) <- param_names
    }

    # Transform parameters from log scale
    r <- exp(par[["log_r"]])
    K <- exp(par[["log_K"]])
    m <- exp(par[["log_m"]])
    q <- exp(par[["log_q"]])
    sigma_obs <- exp(par[["log_sigma_obs"]])
    B0 <- exp(par[["log_B0"]])

    # Data dimensions
    n_years <- length(data$years)

    # Initialize biomass trajectory
    B <- numeric(n_years)
    B[1] <- B0

    # Initialize negative log-likelihood
    nll <- 0

    # Calculate biomass dynamics
    for (t in 1:(n_years - 1)) {
      # Production using Pella-Tomlinson
      if (B[t] > 0 && K > 0 && m > 0) {
        production <- r * B[t] * (1 - (B[t] / K)^(m - 1)) / m
      } else {
        production <- 0
      }

      # Next year's biomass
      B_next <- B[t] + production - data$catch[t]
      B[t + 1] <- max(B_next, 0.01) # Prevent negative biomass
    }

    # Calculate observation likelihood (CPUE)
    for (t in 1:n_years) {
      if (!is.na(data$cpue[t]) && data$cpue[t] > 0 && B[t] > 0) {
        # Expected CPUE
        cpue_pred <- q * B[t]

        if (cpue_pred > 0) {
          # Log-normal likelihood
          log_cpue_obs <- log(data$cpue[t])
          log_cpue_pred <- log(cpue_pred)

          nll <- nll + 0.5 * log(2 * pi * sigma_obs^2) +
            0.5 * (log_cpue_obs - log_cpue_pred)^2 / sigma_obs^2
        } else {
          nll <- nll + 1000 # Large penalty
        }
      }
    }

    # Add penalties for boundary conditions
    if (any(B <= 0.01)) {
      nll <- nll + 1000 * sum(B <= 0.01)
    }

    # Parameter bounds penalties
    if (r <= 0 || r > 2) nll <- nll + 1000
    if (K <= 0 || K > 1e6) nll <- nll + 1000
    if (m <= 0 || m > 10) nll <- nll + 1000
    if (q <= 0 || q > 1) nll <- nll + 1000
    if (sigma_obs <= 0 || sigma_obs > 2) nll <- nll + 1000
    if (B0 <= 0 || B0 > K * 2) nll <- nll + 1000

    return(nll)
  }
}
#' @keywords internal
rtmb_objective <- function(pars, data) {
    
    getAll(data, pars, warn = FALSE)
    
    message("ADFun 1")
    # Transform parameters from log scale
    r         <- exp(log_r)
    K         <- exp(log_K)
    m         <- exp(log_m)
    q         <- exp(log_q)
    sigma_obs <- exp(log_sigma_obs)
    B0         <- exp(log_B0)
    
    # observed values
    cpue <- OBS(cpue)
    message("ADFun 1b")
    # Data dimensions
    n_years <- length(years)
    
    # Initialize biomass trajectory
    B <- AD(numeric(n_years))
    #ADREPORT(B)
    #
    B[1] <- B0

    # Calculate biomass dynamics
    for (t in 1:(n_years - 1)) {
        
        # Production using Pella-Tomlinson
        #production <- r * B[t] * (1 - (B[t] / K))
        message("ADFun 1d")
        # Next year's biomass
        B[t + 1] <- B[t] + r * B[t] * (1 - B[t] / K) - catch[t]
    }
    message("ADFun 2")
    #REPORT(B)
    # Calculate observation likelihood (CPUE)
    cpue_pred <- AD(numeric(n_years))
    #cpue_loc  <- logical(n_years)
    for (t in 1:n_years) {
        #if (!is.na(cpue[t]) && cpue[t] > 0 && B[t] > 0) {
            
            # Expected CPUE
            cpue_pred[t] <- q * B[t]
            message("ADFun 3a ")
            # location of data in vector
            #cpue_loc[t] <- 1 > 0
        #}
    }
    message("ADFun 3a")
    REPORT(cpue_pred)
    #cpue[cpue_loc] %~% dlnorm(log(cpue_pred[cpue_loc]) - sigma_obs^2 / 2, sigma_obs) 
    #cpue %~% dlnorm(log(cpue_pred), sdlog = 0.1) 
    
    #nll_cpue <- -dlnorm(x = cpue, meanlog = log(cpue_pred), sdlog = sigma_obs, log = TRUE)
    nll_cpue <- -dlnorm(x = cpue, mean = log(cpue_pred), sd = 0.1, log = TRUE)
    
    message("ADFun 3b")
    #REPORT(cpue_pred)
    return(sum(nll_cpue))
}

#' Generate Starting Parameter Values
#'
#' Generate reasonable starting values for Pella-Tomlinson model parameters
#' based on data characteristics.
#'
#' @param data List containing model data with 'years', 'catch', and 'cpue'
#'
#' @return Named list of log-transformed starting parameter values
#'
#' @details
#' Starting values are estimated using simple heuristics:
#' - r: Estimated from CPUE trend (rough productivity estimate)
#' - K: Estimated as multiple of maximum observed biomass proxy
#' - m: Default to 2 (Schaefer model)
#' - q: Estimated from CPUE and biomass proxy relationship
#' - sigma_proc: Default moderate process error
#' - sigma_obs: Estimated from CPUE variability
#' - B0: Estimated from initial CPUE and catchability proxy
#'
#' @keywords internal
generate_starting_values <- function(data) {
  n_years <- length(data$years)

  # Estimate intrinsic growth rate from CPUE trend
  cpue_trend <- lm(log(data$cpue) ~ data$years)
  r_est <- abs(coef(cpue_trend)[2]) * 4 # Scale up for reasonable r
  r_est <- max(0.05, min(r_est, 1.0)) # Constrain to reasonable range

  # Estimate carrying capacity as multiple of biomass proxy
  # Use average catch / average harvest rate proxy
  avg_catch <- mean(data$catch, na.rm = TRUE)
  K_est <- avg_catch * 10 # Rough multiple, assuming moderate exploitation

  # Shape parameter (default to Schaefer)
  m_est <- 2.0

  # Estimate catchability from CPUE and biomass relationship
  # Assume CPUE = q * B, and B ~ K/2 on average
  q_est <- mean(data$cpue, na.rm = TRUE) / (K_est / 2)
  q_est <- max(1e-6, min(q_est, 1.0)) # Reasonable bounds

  # Process error - moderate default
  sigma_proc_est <- 0.2

  # Observation error from CPUE variability
  sigma_obs_est <- sd(log(data$cpue), na.rm = TRUE)
  sigma_obs_est <- max(0.1, min(sigma_obs_est, 1.0)) # Reasonable bounds

  # Initial biomass estimate
  B0_est <- data$cpue[1] / q_est # From first CPUE observation
  B0_est <- max(K_est * 0.3, min(B0_est, K_est * 1.2)) # Reasonable bounds

  # Return log-transformed parameters
  starting_values <- list(
    log_r = log(r_est),
    log_K = log(K_est),
    log_m = log(m_est),
    log_q = log(q_est),
    log_sigma_proc = log(sigma_proc_est),
    log_sigma_obs = log(sigma_obs_est),
    log_B0 = log(B0_est)
  )

  return(starting_values)
}

#' Calculate Pella-Tomlinson Production Function
#'
#' Calculate production at given biomass level using Pella-Tomlinson formulation.
#'
#' @param B Biomass level (tonnes)
#' @param r Intrinsic growth rate (per year)
#' @param K Carrying capacity (tonnes)
#' @param m Shape parameter (dimensionless)
#'
#' @return Production value (tonnes per year)
#'
#' @details
#' Implements the Pella-Tomlinson production function:
#' P(B) = r * B * (1 - (B/K)^(m-1)) / m
#'
#' Special cases:
#' - m = 1: Fox model
#' - m = 2: Schaefer model
#'
#' @examples
#' \dontrun{
#' # Schaefer model example
#' production <- pella_tomlinson_production(1000, r = 0.3, K = 5000, m = 2)
#' }
#'
#' @export
pella_tomlinson_production <- function(B, r, K, m) {
  # Input validation
  if (any(c(B, r, K, m) < 0)) {
    stop("All parameters must be non-negative")
  }

  if (K == 0) {
    stop("Carrying capacity K cannot be zero")
  }

  if (m == 0) {
    stop("Shape parameter m cannot be zero")
  }

  # Calculate production
  if (B == 0) {
    return(0)
  }

  # Handle potential numerical issues
  if (B > K) {
    return(0) # No production above carrying capacity
  }

  # Pella-Tomlinson production function
  production <- r * B * (1 - (B / K)^(m - 1)) / m

  return(production)
}
