#' RTMB Objective Function for Pella-Tomlinson Model
#'
#' This file contains the RTMB objective function implementation for the
#' Pella-Tomlinson surplus production model with state-space formulation.
#'
#' @name rtmb-objective
NULL

#' Create RTMB Objective Function for Pella-Tomlinson Model
#'
#' Creates an RTMB-compatible objective function for the Pella-Tomlinson
#' surplus production model.  When \code{proc_dev} random effects are supplied
#' in the parameter list, the function implements a state-space formulation
#' where log-scale process deviations are integrated out via the Laplace
#' approximation.
#'
#' @param data_env Environment (or list coerced to environment) holding the
#'   model data visible to the objective function.  Required elements:
#'   \code{n_years}, \code{n_areas}, \code{catch_mat}, \code{cpue_obs},
#'   \code{areas}, and optionally movement matrices.
#'
#' @return A function of a single argument \code{parms} (a named list)
#'   suitable for passing to \code{RTMB::MakeADFun}.
#'
#' @details
#' The function uses \code{RTMB::getAll()} to unpack parameters and
#' \code{RTMB::ADREPORT()} to report derived quantities (biomass trajectory,
#' harvest rates, fitted CPUE, natural-scale parameters).
#'
#' Production function:
#' \deqn{P(B) = r \cdot B \cdot \bigl(1 - (B/K)^{m-1}\bigr) / m}
#'
#' State equation (deterministic):
#' \deqn{B_{t+1} = B_t + P(B_t) - C_t}
#'
#' State equation (state-space, when \code{proc_dev} random effects present):
#' \deqn{\log B_{t+1} = \log\bigl(B_t + P(B_t) - C_t\bigr) + \epsilon_t}
#' where \eqn{\epsilon_t \sim N(0, \sigma_{\text{proc}}^2)}.
#'
#' Observation equation (lognormal):
#' \deqn{\log \text{CPUE}_t = \log(q \cdot B_t) + \eta_t, \quad
#'        \eta_t \sim N(0, \sigma_{\text{obs}}^2)}
#'
#' @keywords internal
create_rtmb_objective <- function(data_env) {
  # Coerce list to environment so closure can capture it

  if (is.list(data_env)) data_env <- list2env(data_env, parent = emptyenv())

  function(parms) {
    # Unpack ALL parameters into local scope (RTMB-tracked AD variables)
    RTMB::getAll(parms, warn = FALSE)

    # ---- Transform global parameters from log scale ----
    r         <- exp(log_r)
    K         <- exp(log_K)
    m         <- exp(log_m)
    sigma_obs <- exp(log_sigma_obs)

    n_years <- data_env$n_years
    n_areas <- data_env$n_areas
    areas   <- data_env$areas

    has_proc <- exists("log_sigma_proc", inherits = FALSE) &&
                exists("proc_dev", inherits = FALSE)
    if (has_proc) {
      sigma_proc <- exp(log_sigma_proc)
    }

    # ---- Extract per-area B0 and q ----
    # IMPORTANT: In closures passed to MakeADFun, element assignment into
    # plain numeric() vectors silently loses the advector class.  We must
    # pre-allocate with RTMB::advector() so the container is AD-aware.
    B0_vec <- RTMB::advector(numeric(n_areas))
    for (ia in seq_len(n_areas)) {
      key <- paste0("log_B0_", areas[ia])
      B0_vec[ia] <- exp(parms[[key]])
    }

    # q – may be per-area or per-area-label
    has_labels <- !is.null(data_env$labels)
    if (has_labels) {
      labels <- data_env$labels
      n_labels <- length(labels)
      q_arr <- RTMB::advector(matrix(0, n_areas, n_labels))
      for (ia in seq_len(n_areas)) {
        for (il in seq_len(n_labels)) {
          key <- paste0("log_q_", areas[ia], "_", labels[il])
          q_arr[ia, il] <- exp(parms[[key]])
        }
      }
    } else {
      q_vec <- RTMB::advector(numeric(n_areas))
      for (ia in seq_len(n_areas)) {
        key <- paste0("log_q_", areas[ia])
        q_vec[ia] <- exp(parms[[key]])
      }
    }

    # ---- Retrieve data ----
    catch_mat <- data_env$catch_mat   # [n_years x n_areas]
    cpue_obs  <- data_env$cpue_obs    # matrix or 3-d array

    # ---- Biomass dynamics ----
    # Use a list-of-vectors to store biomass (avoids matrix ops that strip
    # RTMB AD class attributes).  B[[t]][ia] = biomass in year t, area ia.
    B <- vector("list", n_years)
    B[[1]] <- B0_vec

    nll <- 0

    for (t in seq_len(n_years - 1)) {
      # Element-wise Pella-Tomlinson production per area
      B_next <- RTMB::advector(numeric(n_areas))
      for (ia in seq_len(n_areas)) {
        Bt_ia <- B[[t]][ia]
        prod_ia <- r * Bt_ia * (1 - (Bt_ia / K)^(m - 1)) / m
        B_det_ia <- Bt_ia + prod_ia - catch_mat[t, ia]
        # Soft lower bound: RTMB cannot branch on AD types, so use
        # a differentiable approximation to max(x, eps).
        # sqrt(x^2 + eps^2) ≈ |x| for |x| >> eps.
        B_det_ia <- 0.5 * (B_det_ia + sqrt(B_det_ia * B_det_ia + 4e-8)) + 1e-8

        if (has_proc) {
          log_B_next_ia <- log(B_det_ia) + proc_dev[t, ia]
          nll <- nll - RTMB::dnorm(proc_dev[t, ia], 0, sigma_proc, log = TRUE)
          B_next[ia] <- exp(log_B_next_ia)
        } else {
          B_next[ia] <- B_det_ia
        }
      }
      B[[t + 1]] <- B_next
    }

    # Apply gravity movement (if present) --- deterministic redistribution
    if (!is.null(data_env$movement_rate) && data_env$movement_rate > 0) {
      move_rate <- data_env$movement_rate
      decay_val <- data_env$decay
      dist_mat  <- data_env$distance_matrix
      attract   <- data_env$attractiveness
      nA <- n_areas
      # Build movement kernel (constant, plain R – no AD)
      W <- matrix(0, nA, nA)
      for (a in seq_len(nA)) {
        for (b in seq_len(nA)) {
          W[a, b] <- attract[b] * exp(-decay_val * dist_mat[a, b])
        }
      }
      rs <- rowSums(W)
      Kmat <- W / rs
      # Redistribute biomass for years 2..n_years
      for (t in 2:n_years) {
        Bt_old <- B[[t]]
        Bt_new <- RTMB::advector(numeric(nA))
        for (ia in seq_len(nA)) {
          moved <- 0
          for (ib in seq_len(nA)) {
            moved <- moved + Kmat[ia, ib] * Bt_old[ib]
          }
          Bt_new[ia] <- (1 - move_rate) * Bt_old[ia] + move_rate * moved
          # Differentiable lower bound (no if-branch on AD types)
          Bt_new[ia] <- 0.5 * (Bt_new[ia] + sqrt(Bt_new[ia] * Bt_new[ia] + 4e-8)) + 1e-8
        }
        B[[t]] <- Bt_new
      }
    }

    # ---- Observation likelihood (lognormal CPUE) ----
    if (has_labels) {
      labels <- data_env$labels
      n_labels <- length(labels)
      for (t in seq_len(n_years)) {
        for (ia in seq_len(n_areas)) {
          for (il in seq_len(n_labels)) {
            cpue_val <- cpue_obs[t, ia, il]
            if (!is.na(cpue_val) && cpue_val > 0) {
              cpue_pred <- q_arr[ia, il] * B[[t]][ia]
              nll <- nll - RTMB::dnorm(log(cpue_val), log(cpue_pred), sigma_obs, log = TRUE)
            }
          }
        }
      }
    } else {
      for (t in seq_len(n_years)) {
        for (ia in seq_len(n_areas)) {
          cpue_val <- cpue_obs[t, ia]
          if (!is.na(cpue_val) && cpue_val > 0) {
            cpue_pred <- q_vec[ia] * B[[t]][ia]
            nll <- nll - RTMB::dnorm(log(cpue_val), log(cpue_pred), sigma_obs, log = TRUE)
          }
        }
      }
    }

    # ---- Convert B list to matrix for ADREPORT ------------------------
    B_mat <- RTMB::advector(matrix(0, n_years, n_areas))
    for (t in seq_len(n_years)) {
      for (ia in seq_len(n_areas)) {
        B_mat[t, ia] <- B[[t]][ia]
      }
    }

    # ---- ADREPORT derived quantities (for sdreport SEs) ----
    RTMB::ADREPORT(B_mat)
    RTMB::ADREPORT(r)
    RTMB::ADREPORT(K)
    RTMB::ADREPORT(m)
    RTMB::ADREPORT(sigma_obs)
    if (has_proc) RTMB::ADREPORT(sigma_proc)

    # Harvest rate (advector for ADREPORT)
    harvest_rate <- RTMB::advector(matrix(0, n_years, n_areas))
    for (t in seq_len(n_years)) {
      for (ia in seq_len(n_areas)) {
        harvest_rate[t, ia] <- catch_mat[t, ia] / B_mat[t, ia]
      }
    }
    RTMB::ADREPORT(harvest_rate)

    # Fitted CPUE (advector for ADREPORT)
    if (has_labels) {
      fitted_cpue <- RTMB::advector(array(0, dim = c(n_years, n_areas, n_labels)))
      for (t in seq_len(n_years)) {
        for (ia in seq_len(n_areas)) {
          for (il in seq_len(n_labels)) {
            fitted_cpue[t, ia, il] <- q_arr[ia, il] * B_mat[t, ia]
          }
        }
      }
      RTMB::ADREPORT(fitted_cpue)
    } else {
      fitted_cpue <- RTMB::advector(matrix(0, n_years, n_areas))
      for (t in seq_len(n_years)) {
        for (ia in seq_len(n_areas)) {
          fitted_cpue[t, ia] <- q_vec[ia] * B_mat[t, ia]
        }
      }
      RTMB::ADREPORT(fitted_cpue)
    }

    nll
  }
}

#' Create Simple Objective Function (non-RTMB fallback)
#'
#' Creates a plain R objective function for parameter estimation using maximum
#' likelihood estimation without automatic differentiation.  This is used as a
#' fallback when RTMB is not available.
#'
#' @param data List containing model data (from \code{preprocess_model_data}).
#' @param initial_params Named list of starting parameter values (log-transformed)
#'
#' @return Function that takes a parameter vector and returns the negative
#'   log-likelihood scalar.
#'
#' @keywords internal
create_simple_objective <- function(data, initial_params) {
  # Store parameter names for reference
  param_names <- names(initial_params)

  # Safe check: does a named key exist in par?
  has_key <- function(par, key) key %in% names(par)

  # Helper to extract per-area/label parameters (supports single-area fallback)
  get_q_param <- function(par, areas, labels = NULL) {
    if (is.null(labels)) {
      vals <- numeric(length(areas))
      names(vals) <- areas
      for (a in areas) {
        key <- paste0("log_q.", a)
        if (has_key(par, key)) {
          vals[[a]] <- exp(par[[key]])
        } else if (has_key(par, "log_q")) {
          vals[[a]] <- exp(par[["log_q"]])
        } else {
          # fallback: try any log_q.*
          q_keys <- grep("^log_q", names(par), value = TRUE)
          if (length(q_keys) > 0) {
            vals[[a]] <- exp(par[[q_keys[1]]])
          } else {
            # Instead of error, set to default (1) and warn
            warning(paste("Missing parameter:", key, "- using default q=1"))
            vals[[a]] <- 1
          }
        }
      }
      return(vals)
    } else {
      vals <- array(NA_real_, dim = c(length(areas), length(labels)), dimnames = list(area = areas, label = labels))
      for (a in areas) {
        for (l in labels) {
          key <- paste0("log_q.", a, ".", l)
          area_key <- paste0("log_q.", a)
          if (has_key(par, key)) {
            vals[a, l] <- exp(par[[key]])
          } else if (has_key(par, area_key)) {
            vals[a, l] <- exp(par[[area_key]])
          } else if (has_key(par, "log_q")) {
            vals[a, l] <- exp(par[["log_q"]])
          } else {
            # fallback: try any log_q.*
            q_keys <- grep("^log_q", names(par), value = TRUE)
            if (length(q_keys) > 0) {
              vals[a, l] <- exp(par[[q_keys[1]]])
            } else {
              warning(paste("Missing parameter:", key, "- using default q=1"))
              vals[a, l] <- 1
            }
          }
        }
      }
      return(vals)
    }
  }
  get_B0_param <- function(par, areas) {
    vals <- numeric(length(areas))
    names(vals) <- areas
    for (a in areas) {
      key <- paste0("log_B0.", a)
      if (has_key(par, key)) {
        vals[[a]] <- exp(par[[key]])
      } else if (has_key(par, "log_B0")) {
        vals[[a]] <- exp(par[["log_B0"]])
      } else {
        # fallback: try any B0
        b0_keys <- grep("^log_B0", names(par), value = TRUE)
        if (length(b0_keys) > 0) {
          vals[[a]] <- exp(par[[b0_keys[1]]])
        } else {
          stop("Missing parameter: ", key, " and no default provided")
        }
      }
    }
    vals
  }

  # Return objective function
  function(par) {
    if (is.null(names(par))) names(par) <- param_names

    # Transform global parameters from log scale
    r <- exp(par[["log_r"]])
    K <- exp(par[["log_K"]])
    m <- exp(par[["log_m"]])
    sigma_obs <- exp(par[["log_sigma_obs"]])
    sigma_proc <- if (has_key(par, "log_sigma_proc")) exp(par[["log_sigma_proc"]]) else NA_real_

    n_years <- length(data$years)
    areas <- data$areas
    n_areas <- length(areas)
    has_label <- !is.null(data$labels)
    labels <- if (has_label) data$labels else NULL

    # Extract per-area/label parameters
    q_param <- get_q_param(par, areas, labels)
    B0_vec <- get_B0_param(par, areas)

    # Initialize biomass trajectory matrix [year x area]
    B <- matrix(NA_real_,
      nrow = n_years, ncol = n_areas,
      dimnames = list(as.character(data$years), areas)
    )
    B[1, ] <- as.numeric(B0_vec)

    # Initialize negative log-likelihood
    nll <- 0

    # Catch as matrix
    catch_mat <- if (is.matrix(data$catch)) data$catch else matrix(data$catch, ncol = 1, dimnames = list(NULL, areas))

    # Calculate biomass dynamics per area, with optional gravity movement
    for (t in 1:(n_years - 1)) {
      Bt <- B[t, ]
      production <- ifelse(Bt > 0 & K > 0 & m > 0,
        r * Bt * (1 - (Bt / K)^(m - 1)) / m,
        0
      )
      # Guard against NaN/Inf from extreme parameter combinations
      production[!is.finite(production)] <- 0
      B_next <- Bt + production - catch_mat[t, ]

      # Gravity movement kernel (if enabled)
      move_rate <- data$movement_rate %||% 0
      decay <- data$decay %||% 1
      dist_mat <- data$distance_matrix
      attract <- data$attractiveness
      if (!is.null(move_rate) && move_rate > 0 && !is.null(dist_mat) && !is.null(attract)) {
        # Build movement weights W_{ab} = A_b * exp(-decay * D_{ab})
        nA <- length(B_next)
        W <- matrix(0, nA, nA)
        for (a in seq_len(nA)) {
          for (b in seq_len(nA)) {
            W[a, b] <- attract[b] * exp(-decay * dist_mat[a, b])
          }
        }
        # Row-normalize weights
        Kmat <- W / rowSums(W)
        # Redistribute biomass: B_a'[t+1] = (1-move_rate)*B_a[t+1] + move_rate*sum_b B_b[t+1]*K_{ba}
        B_next <- (1 - move_rate) * B_next + move_rate * as.numeric(Kmat %*% B_next)
      }
      B[t + 1, ] <- pmax(B_next, 0.01)
    }

    # Calculate observation likelihood (CPUE) per area/label
    if (has_label) {
      cpue_arr <- data$cpue
      for (t in 1:n_years) {
        Bt <- B[t, ]
        for (a in seq_along(areas)) {
          for (l in seq_along(labels)) {
            cpue_obs <- cpue_arr[t, a, l]
            q_val <- q_param[a, l]
            cpue_pred <- q_val * Bt[a]
            valid <- !is.na(cpue_obs) && cpue_obs > 0 && Bt[a] > 0 && cpue_pred > 0
            if (valid) {
              log_obs <- log(cpue_obs)
              log_pred <- log(cpue_pred)
              res2 <- (log_obs - log_pred)^2
              nll <- nll + 0.5 * log(2 * pi * sigma_obs^2) + 0.5 * res2 / sigma_obs^2
            }
          }
        }
      }
    } else {
      cpue_mat <- if (is.matrix(data$cpue)) data$cpue else matrix(data$cpue, ncol = 1, dimnames = list(NULL, areas))
      for (t in 1:n_years) {
        Bt <- B[t, ]
        cpue_obs <- cpue_mat[t, ]
        cpue_pred <- as.numeric(q_param) * Bt
        valid <- !is.na(cpue_obs) & cpue_obs > 0 & Bt > 0 & cpue_pred > 0
        if (any(valid)) {
          log_obs <- log(cpue_obs[valid])
          log_pred <- log(cpue_pred[valid])
          res2 <- (log_obs - log_pred)^2
          nll <- nll + sum(0.5 * log(2 * pi * sigma_obs^2) + 0.5 * res2 / sigma_obs^2)
        }
      }
    }

    # Add penalties for boundary conditions (NA-safe)
    bad_B <- !is.finite(B) | B <= 0.01
    if (any(bad_B)) {
      nll <- nll + 1000 * sum(bad_B)
    }

    # Parameter bounds penalties (globals)
    if (!is.finite(r) || r <= 0 || r > 2) nll <- nll + 1000
    if (!is.finite(K) || K <= 0 || K > 1e9) nll <- nll + 1000
    if (!is.finite(m) || m <= 0 || m > 10) nll <- nll + 1000
    if (!is.finite(sigma_obs) || sigma_obs <= 0 || sigma_obs > 2) nll <- nll + 1000
    # Per-area/label bounds (NA-safe)
    bad_q <- !is.finite(q_param) | q_param <= 0 | q_param > 1
    if (any(bad_q)) nll <- nll + 1000 * sum(bad_q)
    bad_B0 <- !is.finite(B0_vec) | B0_vec <= 0 | B0_vec > K * 2
    if (any(bad_B0)) nll <- nll + 500 * sum(bad_B0)

    # Final guard: if nll became NaN/Inf return large penalty
    if (!is.finite(nll)) nll <- 1e8

    return(nll)
  }
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
  multi_area <- is.matrix(data$cpue) || is.matrix(data$catch)
  has_label <- !is.null(data$labels)
  if (has_label) {
    areas <- data$areas
  } else {
    areas <- if (!is.null(data$areas)) data$areas else if (multi_area) colnames(data$cpue) else "A1"
  }
  n_areas <- length(areas)

  # Coerce to matrices for simple handling
  catch_mat <- if (is.matrix(data$catch)) data$catch else matrix(data$catch, ncol = 1)
  if (has_label) {
    # collapse across labels for trend/variance heuristics
    cpue_mat <- apply(data$cpue, c(1, 2), function(x) mean(x, na.rm = TRUE))
  } else {
    cpue_mat <- if (multi_area) data$cpue else matrix(data$cpue, ncol = 1)
  }

  # Estimate intrinsic growth rate from CPUE trend (use total index if multi-area)
  cpue_total <- rowSums(cpue_mat, na.rm = TRUE)
  cpue_trend <- try(lm(log(pmax(cpue_total, 1e-8)) ~ data$years), silent = TRUE)
  r_est <- if (inherits(cpue_trend, "try-error")) 0.2 else abs(coef(cpue_trend)[2]) * 4
  r_est <- max(0.05, min(r_est, 1.0)) # Constrain to reasonable range

  # Estimate carrying capacity using MSY heuristic (Schaefer: MSY = rK/4)
  # max_catch is a lower bound on MSY, so K >= 4*max_catch/r
  max_catch <- max(catch_mat, na.rm = TRUE)
  avg_catch <- mean(catch_mat, na.rm = TRUE)
  K_from_r <- 4 * max_catch / r_est
  K_from_catch <- avg_catch * 30 # broader multiplier for slow-growing species
  K_est <- max(100, sqrt(K_from_r * K_from_catch)) # geometric mean

  # Shape parameter (default to Schaefer)
  m_est <- 2.0

  # Process error - moderate default
  sigma_proc_est <- 0.2

  # Observation error from CPUE variability (use all cpue entries)
  sigma_obs_est <- sd(log(as.numeric(cpue_mat)), na.rm = TRUE)
  sigma_obs_est <- max(0.1, min(ifelse(is.na(sigma_obs_est), 0.3, sigma_obs_est), 1.0))

  # Per-area catchability and initial biomass
  starting_values <- list(
    log_r = log(r_est),
    log_K = log(K_est),
    log_m = log(m_est),
    log_sigma_proc = log(sigma_proc_est),
    log_sigma_obs = log(sigma_obs_est)
  )

  if (n_areas == 1 && !has_label) {
    cpue_a <- cpue_mat[, 1]
    mean_cpue_a <- mean(cpue_a, na.rm = TRUE)
    q_est_a <- mean_cpue_a / (K_est / 2)
    q_est_a <- max(1e-8, min(ifelse(is.na(q_est_a) || !is.finite(q_est_a), 1e-3, q_est_a), 1.0))
    first_non_na <- which(!is.na(cpue_a))[1]
    cpue0 <- if (!is.na(first_non_na)) cpue_a[first_non_na] else mean_cpue_a
    if (is.na(cpue0) || !is.finite(cpue0)) cpue0 <- max(mean_cpue_a, 1e-6)
    B0_est_a <- pmax(K_est * 0.1, pmin(cpue0 / q_est_a, K_est * 1.2))
    starting_values[["log_q"]] <- log(q_est_a)
    starting_values[["log_B0"]] <- log(B0_est_a)
  } else if (!has_label) {
    for (j in seq_along(areas)) {
      a <- areas[[j]]
      cpue_a <- cpue_mat[, j]
      mean_cpue_a <- mean(cpue_a, na.rm = TRUE)
      q_est_a <- mean_cpue_a / (K_est / 2)
      q_est_a <- max(1e-8, min(ifelse(is.na(q_est_a) || !is.finite(q_est_a), 1e-3, q_est_a), 1.0))
      first_non_na <- which(!is.na(cpue_a))[1]
      cpue0 <- if (!is.na(first_non_na)) cpue_a[first_non_na] else mean_cpue_a
      if (is.na(cpue0) || !is.finite(cpue0)) cpue0 <- max(mean_cpue_a, 1e-6)
      B0_est_a <- pmax(K_est * 0.1, pmin(cpue0 / q_est_a, K_est * 1.2))
      starting_values[[paste0("log_q.", a)]] <- log(q_est_a)
      starting_values[[paste0("log_B0.", a)]] <- log(B0_est_a)
    }
  } else {
    # per area and label q, and per area B0
    labels <- data$labels
    # B0 by area using mean across labels
    for (j in seq_along(areas)) {
      a <- areas[[j]]
      cpue_a_all <- as.numeric(data$cpue[, j, ])
      mean_cpue_a <- mean(cpue_a_all, na.rm = TRUE)
      # derive a generic B0 for area a
      q_tmp <- max(1e-6, min(mean_cpue_a / (K_est / 2), 1.0))
      cpue0 <- mean_cpue_a
      B0_est_a <- pmax(K_est * 0.1, pmin(cpue0 / q_tmp, K_est * 1.2))
      starting_values[[paste0("log_B0.", a)]] <- log(B0_est_a)
      # now per-label q
      for (l in seq_along(labels)) {
        mean_cpue_al <- mean(data$cpue[, j, l], na.rm = TRUE)
        q_est_al <- mean_cpue_al / (K_est / 2)
        q_est_al <- max(1e-8, min(ifelse(is.na(q_est_al) || !is.finite(q_est_al), 1e-3, q_est_al), 1.0))
        starting_values[[paste0("log_q.", a, ".", labels[[l]])]] <- log(q_est_al)
      }
    }
  }

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
