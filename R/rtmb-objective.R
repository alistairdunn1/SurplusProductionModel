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
#' Production function (standard Pella-Tomlinson):
#' \deqn{P(B) = \frac{r}{m-1} \cdot B \cdot \bigl(1 - (B/K)^{m-1}\bigr)}
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
    r <- exp(log_r)
    m <- exp(log_m)
    sigma_obs <- exp(log_sigma_obs)

    n_years <- data_env$n_years
    n_areas <- data_env$n_areas
    areas <- data_env$areas

    has_proc <- exists("log_sigma_proc", inherits = FALSE) &&
      exists("proc_dev", inherits = FALSE)
    if (has_proc) {
      sigma_proc <- exp(log_sigma_proc)
    }

    has_movement_inputs <- !is.null(data_env$distance_matrix) &&
      !is.null(data_env$attractiveness) &&
      !is.null(data_env$decay)
    has_movement_param <- "log_movement_rate" %in% names(parms)
    if (has_movement_param) {
      move_rate <- 1 / (1 + exp(-parms[["log_movement_rate"]]))
    } else {
      move_rate <- data_env$movement_rate %||% 0
    }

    process_error_structure <- data_env$process_error_structure %||% "iid"
    use_ar1 <- has_proc && identical(process_error_structure, "ar1")
    if (use_ar1) {
      rho <- tanh(parms[["theta_rho"]])
    }

    has_env <- !is.null(data_env$env_array) && !is.null(data_env$env_covariate_names)
    if (has_env) {
      env_array <- data_env$env_array # [n_years-1 x n_areas x n_cov]
      cov_names <- data_env$env_covariate_names
      n_cov <- length(cov_names)
      beta_cov <- RTMB::advector(numeric(n_cov))
      for (ic in seq_len(n_cov)) {
        beta_name <- paste0("beta_", cov_names[ic])
        beta_cov[ic] <- parms[[beta_name]]
      }
    }

    # ---- Extract initial depletion and per-area carrying capacity ----
    d0 <- exp(parms[["log_d0"]])
    K_vec <- RTMB::advector(numeric(n_areas))
    B_initial_vec <- RTMB::advector(numeric(n_areas))
    for (ia in seq_len(n_areas)) {
      k_key <- paste0("log_K_", areas[ia])
      if (k_key %in% names(parms)) {
        K_vec[ia] <- exp(parms[[k_key]])
      } else if (n_areas == 1L && "log_K" %in% names(parms)) {
        K_vec[ia] <- exp(parms[["log_K"]])
      } else {
        stop("Missing carrying-capacity parameter for area ", areas[ia])
      }
      # Spin-up starts from the undepleted per-area carrying capacity; initial
      # depletion d0 is applied after the unfished (movement) equilibrium is
      # resolved, so d0 is not erased by the spin-up.
      B_initial_vec[ia] <- K_vec[ia]
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
    catch_mat <- data_env$catch_mat # [n_years x n_areas]
    cpue_obs <- data_env$cpue_obs # matrix or 3-d array
    cpue_sd <- data_env$cpue_sd
    spinup_years <- as.integer(data_env$spinup_years %||% 0L)

    if (!has_labels && !is.matrix(cpue_obs)) {
      cpue_obs <- matrix(cpue_obs, ncol = n_areas, dimnames = list(NULL, areas))
    }
    if (!has_labels && !is.null(cpue_sd) && !is.matrix(cpue_sd)) {
      cpue_sd <- matrix(cpue_sd, ncol = n_areas, dimnames = list(NULL, areas))
    }

    # Pre-build movement kernel when movement inputs are available.
    if (has_movement_inputs) {
      decay_val <- data_env$decay
      dist_mat <- data_env$distance_matrix
      attract <- data_env$attractiveness
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
      rs <- rowSums(W)
      Kmat <- W / rs
    }

    # ---- Biomass dynamics ----
    # Use a list-of-vectors to store biomass (avoids matrix ops that strip
    # RTMB AD class attributes).  B[[t]][ia] = biomass in year t, area ia.
    B <- vector("list", n_years)

    # Resolve the unfished spatial equilibrium, then apply initial depletion.
    # With movement the joint zero-catch equilibrium differs from the per-area
    # carrying capacities because biomass is redistributed between areas; the
    # deterministic spin-up (from K, zero catch) converges to that equilibrium.
    # Without movement each area sits at its own K and the spin-up is a no-op.
    # Initial depletion d0 is applied as a multiplier afterwards, so the spin-up
    # resolves the movement/unfished level while d0 sets depletion.
    B_unfished_vec <- B_initial_vec
    if (spinup_years > 0) {
      B_spin <- B_initial_vec
      for (s in seq_len(spinup_years)) {
        B_next <- RTMB::advector(numeric(n_areas))
        for (ia in seq_len(n_areas)) {
          Bt_ia <- B_spin[ia]
          # Standard Pella-Tomlinson production: P = r/(m-1) * B * (1 - (B/K)^(m-1)).
        # The AD path cannot branch on m, so the Fox limit (m = 1) is not
        # special-cased here; keep m bounded away from 1 when estimating it.
        prod_ia <- r * Bt_ia * (1 - (Bt_ia / K_vec[ia])^(m - 1)) / (m - 1)
          B_det_ia <- Bt_ia + prod_ia
          B_det_ia <- 0.5 * (B_det_ia + sqrt(B_det_ia * B_det_ia + 4e-8)) + 1e-8
          B_next[ia] <- B_det_ia
        }

        if (has_movement_inputs) {
          B_move <- RTMB::advector(numeric(n_areas))
          for (ia in seq_len(n_areas)) {
            moved <- 0
            # Inflow to area ia = sum over sources ib of the share moving
            # FROM ib TO ia, i.e. Kmat[ib, ia] (conserves total biomass).
            for (ib in seq_len(n_areas)) {
              moved <- moved + Kmat[ib, ia] * B_next[ib]
            }
            B_move[ia] <- (1 - move_rate) * B_next[ia] + move_rate * moved
            B_move[ia] <- 0.5 * (B_move[ia] + sqrt(B_move[ia] * B_move[ia] + 4e-8)) + 1e-8
          }
          B_spin <- B_move
        } else {
          B_spin <- B_next
        }
      }
      B_unfished_vec <- B_spin
    }
    B[[1]] <- d0 * B_unfished_vec

    nll <- 0

    for (t in seq_len(n_years - 1)) {
      # Element-wise Pella-Tomlinson production per area
      B_next <- RTMB::advector(numeric(n_areas))
      for (ia in seq_len(n_areas)) {
        Bt_ia <- B[[t]][ia]
        # Standard Pella-Tomlinson production: P = r/(m-1) * B * (1 - (B/K)^(m-1)).
        # The AD path cannot branch on m, so the Fox limit (m = 1) is not
        # special-cased here; keep m bounded away from 1 when estimating it.
        prod_ia <- r * Bt_ia * (1 - (Bt_ia / K_vec[ia])^(m - 1)) / (m - 1)
        B_det_ia <- Bt_ia + prod_ia - catch_mat[t, ia]
        # Soft lower bound: RTMB cannot branch on AD types, so use
        # a differentiable approximation to max(x, eps).
        # sqrt(x^2 + eps^2) ≈ |x| for |x| >> eps.
        B_det_ia <- 0.5 * (B_det_ia + sqrt(B_det_ia * B_det_ia + 4e-8)) + 1e-8

        env_effect_ia <- 0
        if (has_env) {
          for (ic in seq_len(n_cov)) {
            env_effect_ia <- env_effect_ia + env_array[t, ia, ic] * beta_cov[ic]
          }
        }

        if (has_proc) {
          log_B_next_ia <- log(B_det_ia) + env_effect_ia + proc_dev[t, ia]
          if (use_ar1) {
            if (t == 1L) {
              sd0 <- sigma_proc / sqrt(1 - rho * rho + 1e-8)
              nll <- nll - RTMB::dnorm(proc_dev[t, ia], 0, sd0, log = TRUE)
            } else {
              nll <- nll - RTMB::dnorm(proc_dev[t, ia], rho * proc_dev[t - 1, ia], sigma_proc, log = TRUE)
            }
          } else {
            nll <- nll - RTMB::dnorm(proc_dev[t, ia], 0, sigma_proc, log = TRUE)
          }
          B_next[ia] <- exp(log_B_next_ia)
        } else {
          if (has_env) {
            B_next[ia] <- exp(log(B_det_ia) + env_effect_ia)
          } else {
            B_next[ia] <- B_det_ia
          }
        }
      }

      # Apply gravity movement within the annual step so the redistributed
      # biomass becomes the state that enters the next year's production and
      # is seen by the observation model. This couples areas over time
      # (a genuine movement process), rather than reshaping a trajectory that
      # already evolved independently by area.
      if (has_movement_inputs) {
        B_moved <- RTMB::advector(numeric(n_areas))
        for (ia in seq_len(n_areas)) {
          moved <- 0
          # Inflow to area ia = sum over sources ib of the share moving
          # FROM ib TO ia, i.e. Kmat[ib, ia] (conserves total biomass).
          for (ib in seq_len(n_areas)) {
            moved <- moved + Kmat[ib, ia] * B_next[ib]
          }
          B_moved[ia] <- (1 - move_rate) * B_next[ia] + move_rate * moved
          # Differentiable lower bound (no if-branch on AD types)
          B_moved[ia] <- 0.5 * (B_moved[ia] + sqrt(B_moved[ia] * B_moved[ia] + 4e-8)) + 1e-8
        }
        B_next <- B_moved
      }

      B[[t + 1]] <- B_next
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
              obs_sd_ia <- if (!is.null(cpue_sd)) cpue_sd[t, ia, il] else NA_real_
              total_sd <- if (is.finite(obs_sd_ia) && obs_sd_ia >= 0) sqrt(sigma_obs * sigma_obs + obs_sd_ia * obs_sd_ia) else sigma_obs
              nll <- nll - RTMB::dnorm(log(cpue_val), log(cpue_pred), total_sd, log = TRUE)
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
            obs_sd_ia <- if (!is.null(cpue_sd)) cpue_sd[t, ia] else NA_real_
            total_sd <- if (is.finite(obs_sd_ia) && obs_sd_ia >= 0) sqrt(sigma_obs * sigma_obs + obs_sd_ia * obs_sd_ia) else sigma_obs
            nll <- nll - RTMB::dnorm(log(cpue_val), log(cpue_pred), total_sd, log = TRUE)
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
    RTMB::ADREPORT(K_vec)
    # Unfished spatial equilibrium (status baseline). Under movement this is the
    # joint redistributed equilibrium and differs from the per-area K.
    RTMB::ADREPORT(B_unfished_vec)
    RTMB::ADREPORT(m)
    RTMB::ADREPORT(sigma_obs)
    if (has_proc) RTMB::ADREPORT(sigma_proc)
    if (use_ar1) RTMB::ADREPORT(rho)
    if (has_movement_param) RTMB::ADREPORT(move_rate)
    if (has_env) RTMB::ADREPORT(beta_cov)

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

    nll <- .apply_parameter_priors(nll, parms, data_env$priors)

    nll
  }
}

.apply_parameter_priors <- function(nll, parms, priors) {
  if (is.null(priors) || length(priors) == 0) {
    return(nll)
  }

  for (prior in priors) {
    param_value <- parms[[prior$param]]
    if (prior$dist == "normal") {
      nll <- nll - RTMB::dnorm(param_value, mean = prior$mean, sd = prior$sd, log = TRUE)
    } else if (prior$dist == "lognormal") {
      nll <- nll - RTMB::dnorm(param_value, mean = prior$meanlog, sd = prior$sdlog, log = TRUE)
    } else if (prior$dist == "beta") {
      # Parameter is modelled on log scale, so apply beta prior on natural
      # scale x in (0, 1) using a logistic transform with tiny boundary offset.
      s <- 1 / (1 + exp(-param_value))
      eps <- 1e-12
      x <- eps + (1 - 2 * eps) * s
      jac <- (1 - 2 * eps) * s * (1 - s)
      nll <- nll - RTMB::dbeta(x, shape1 = prior$shape1, shape2 = prior$shape2, log = TRUE)
      nll <- nll - log(jac)
    } else if (prior$dist == "exponential") {
      # Parameter is modelled on log scale; apply Exponential(rate) prior on
      # natural scale x = exp(theta) with Jacobian |dx/dtheta| = x.
      x <- exp(param_value)
      nll <- nll - (log(prior$rate) - prior$rate * x + log(x))
    }
  }

  nll
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

  # Helper to extract per-area/label parameters.
  get_q_param <- function(par, areas, labels = NULL) {
    if (is.null(labels)) {
      vals <- numeric(length(areas))
      names(vals) <- areas
      for (a in areas) {
        key <- paste0("log_q.", a)
        if (has_key(par, key)) {
          vals[[a]] <- exp(par[[key]])
        } else if (length(areas) == 1L && has_key(par, "log_q")) {
          vals[[a]] <- exp(par[["log_q"]])
        } else {
          stop("Missing parameter: ", key)
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
          } else {
            stop("Missing parameter: ", key)
          }
        }
      }
      return(vals)
    }
  }
  # Return objective function
  function(par) {
    if (is.null(names(par))) names(par) <- param_names

    # Transform global parameters from log scale
    r <- exp(par[["log_r"]])
    m <- exp(par[["log_m"]])
    sigma_obs <- exp(par[["log_sigma_obs"]])

    n_years <- length(data$years)
    areas <- data$areas
    n_areas <- length(areas)
    has_label <- !is.null(data$labels)
    labels <- if (has_label) data$labels else NULL

    # Extract per-area/label parameters
    q_param <- get_q_param(par, areas, labels)
    if (!has_key(par, "log_d0")) {
      stop("Missing parameter: log_d0")
    }
    d0 <- exp(par[["log_d0"]])
    get_k_param <- function(par, areas) {
      vals <- numeric(length(areas))
      names(vals) <- areas
      for (a in areas) {
        key <- paste0("log_K.", a)
        if (has_key(par, key)) {
          vals[[a]] <- exp(par[[key]])
        } else if (length(areas) == 1L && has_key(par, "log_K")) {
          vals[[a]] <- exp(par[["log_K"]])
        } else {
          stop("Missing parameter: ", key)
        }
      }
      vals
    }
    K_vec_s <- get_k_param(par, areas)
    B_initial_vec <- d0 * K_vec_s

    # Initialize biomass trajectory matrix [year x area]
    B <- matrix(NA_real_,
      nrow = n_years, ncol = n_areas,
      dimnames = list(as.character(data$years), areas)
    )
    B[1, ] <- as.numeric(B_initial_vec)

    # Per-area carrying capacity
    names(K_vec_s) <- areas

    # Initialize negative log-likelihood
    nll <- 0

    # Catch as matrix
    catch_mat <- if (is.matrix(data$catch)) data$catch else matrix(data$catch, ncol = 1, dimnames = list(NULL, areas))

    # Calculate biomass dynamics per area, with optional gravity movement
    for (t in 1:(n_years - 1)) {
      Bt <- B[t, ]
      production <- .pt_production(Bt, r, K_vec_s, m)
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

    cpue_sd_data <- data$cpue_sd %||% NULL
    if (!has_label && !is.null(cpue_sd_data) && !is.matrix(cpue_sd_data)) {
      cpue_sd_data <- matrix(cpue_sd_data, ncol = 1, dimnames = list(NULL, areas))
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
              obs_sd <- if (!is.null(cpue_sd_data)) cpue_sd_data[t, a, l] else NA_real_
              total_sd <- if (is.finite(obs_sd) && obs_sd >= 0) sqrt(sigma_obs^2 + obs_sd^2) else sigma_obs
              log_obs <- log(cpue_obs)
              log_pred <- log(cpue_pred)
              res2 <- (log_obs - log_pred)^2
              nll <- nll + 0.5 * log(2 * pi * total_sd^2) + 0.5 * res2 / total_sd^2
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
          obs_sd <- if (!is.null(cpue_sd_data)) cpue_sd_data[t, ] else rep(NA_real_, length(cpue_obs))
          total_sd <- ifelse(is.finite(obs_sd[valid]) & obs_sd[valid] >= 0, sqrt(sigma_obs^2 + obs_sd[valid]^2), sigma_obs)
          log_obs <- log(cpue_obs[valid])
          log_pred <- log(cpue_pred[valid])
          res2 <- (log_obs - log_pred)^2
          nll <- nll + sum(0.5 * log(2 * pi * total_sd^2) + 0.5 * res2 / total_sd^2)
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
    bad_K <- !is.finite(K_vec_s) | K_vec_s <= 0 | K_vec_s > 1e9
    if (any(bad_K)) nll <- nll + 1000 * sum(bad_K)
    if (!is.finite(m) || m <= 0 || m > 10) nll <- nll + 1000
    if (!is.finite(sigma_obs) || sigma_obs <= 0 || sigma_obs > 2) nll <- nll + 1000
    # Per-area/label bounds (NA-safe)
    bad_q <- !is.finite(q_param) | q_param <= 0 | q_param > 1
    if (any(bad_q)) nll <- nll + 1000 * sum(bad_q)
    if (!is.finite(d0) || d0 <= 0 || d0 > 2) nll <- nll + 500

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
#' - d0: Initial depletion multiplier (B_initial = d0 x K)
#'
#' @keywords internal
generate_starting_values <- function(data) {
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
  cpue_total_vals <- as.numeric(rowSums(cpue_mat, na.rm = TRUE))
  valid <- is.finite(cpue_total_vals) & cpue_total_vals > 0 & is.finite(data$years)
  slope <- NA_real_
  if (sum(valid) >= 2) {
    x <- as.numeric(data$years[valid])
    y <- log(cpue_total_vals[valid])
    fit <- try(stats::lm.fit(cbind(1, x), y), silent = TRUE)
    if (!inherits(fit, "try-error") && length(fit$coefficients) >= 2) {
      slope <- fit$coefficients[2]
    }
  }
  r_est <- if (!is.finite(slope)) 0.2 else abs(slope) * 4
  r_est <- max(0.05, min(r_est, 1.0)) # Constrain to reasonable range

  estimate_k_from_catch <- function(catch_vec, r_est) {
    max_catch <- max(catch_vec, na.rm = TRUE)
    avg_catch <- mean(catch_vec, na.rm = TRUE)
    K_from_r <- 4 * max_catch / r_est
    K_from_catch <- avg_catch * 30
    max(100, sqrt(K_from_r * K_from_catch))
  }

  K_est <- estimate_k_from_catch(as.numeric(catch_mat), r_est)

  # Shape parameter (default to Schaefer)
  m_est <- 2.0

  # Process error - moderate default
  sigma_proc_est <- 0.2

  # Observation error from CPUE variability after removing known per-index SD where available.
  cpue_sd_mat <- data$cpue_sd %||% NULL
  if (!is.null(cpue_sd_mat) && has_label) {
    cpue_sd_mat <- apply(cpue_sd_mat, c(1, 2), function(x) mean(x, na.rm = TRUE))
  } else if (!is.null(cpue_sd_mat) && !is.matrix(cpue_sd_mat)) {
    cpue_sd_mat <- matrix(cpue_sd_mat, ncol = 1)
  }
  total_cpue_sd <- sd(log(as.numeric(cpue_mat)), na.rm = TRUE)
  known_sd2 <- if (!is.null(cpue_sd_mat)) mean(as.numeric(cpue_sd_mat)^2, na.rm = TRUE) else NA_real_
  sigma_obs_est <- if (is.finite(known_sd2)) sqrt(max(total_cpue_sd^2 - known_sd2, 0.01^2)) else total_cpue_sd
  sigma_obs_est <- max(0.01, min(ifelse(is.na(sigma_obs_est), 0.3, sigma_obs_est), 1.0))

  # Per-area catchability and initial depletion
  starting_values <- list(
    log_r = log(r_est),
    log_m = log(m_est),
    log_sigma_proc = log(sigma_proc_est),
    log_sigma_obs = log(sigma_obs_est),
    log_d0 = log(0.8)
  )

  if (n_areas == 1 && !has_label) {
    starting_values[["log_K"]] <- log(K_est)
    cpue_a <- cpue_mat[, 1]
    mean_cpue_a <- mean(cpue_a, na.rm = TRUE)
    q_est_a <- mean_cpue_a / (K_est / 2)
    q_est_a <- max(1e-8, min(ifelse(is.na(q_est_a) || !is.finite(q_est_a), 1e-3, q_est_a), 1.0))
    first_non_na <- which(!is.na(cpue_a))[1]
    cpue0 <- if (!is.na(first_non_na)) cpue_a[first_non_na] else mean_cpue_a
    if (is.na(cpue0) || !is.finite(cpue0)) cpue0 <- max(mean_cpue_a, 1e-6)
    starting_values[["log_q"]] <- log(q_est_a)
  } else if (!has_label) {
    for (j in seq_along(areas)) {
      a <- areas[[j]]
      K_est_a <- estimate_k_from_catch(catch_mat[, j], r_est)
      cpue_a <- cpue_mat[, j]
      mean_cpue_a <- mean(cpue_a, na.rm = TRUE)
      q_est_a <- mean_cpue_a / (K_est_a / 2)
      q_est_a <- max(1e-8, min(ifelse(is.na(q_est_a) || !is.finite(q_est_a), 1e-3, q_est_a), 1.0))
      first_non_na <- which(!is.na(cpue_a))[1]
      cpue0 <- if (!is.na(first_non_na)) cpue_a[first_non_na] else mean_cpue_a
      if (is.na(cpue0) || !is.finite(cpue0)) cpue0 <- max(mean_cpue_a, 1e-6)
      starting_values[[paste0("log_K.", a)]] <- log(K_est_a)
      starting_values[[paste0("log_q.", a)]] <- log(q_est_a)
    }
  } else {
    # per area and label q
    labels <- data$labels
    for (j in seq_along(areas)) {
      a <- areas[[j]]
      K_est_a <- estimate_k_from_catch(catch_mat[, j], r_est)
      starting_values[[paste0("log_K.", a)]] <- log(K_est_a)
      # now per-label q
      for (l in seq_along(labels)) {
        mean_cpue_al <- mean(data$cpue[, j, l], na.rm = TRUE)
        q_est_al <- mean_cpue_al / (K_est_a / 2)
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
#' @param r Productivity parameter (per year). The per-capita growth rate as
#'   \code{B} approaches zero is \code{r / (m - 1)} and \code{FMSY = r / m};
#'   \code{r} is the classical intrinsic growth rate only when \code{m = 2}.
#' @param K Carrying capacity (tonnes)
#' @param m Shape parameter (dimensionless)
#'
#' @return Production value (tonnes per year)
#'
#' @details
#' Implements the standard Pella-Tomlinson production function:
#' P(B) = r / (m - 1) * B * (1 - (B/K)^(m-1))
#'
#' Special cases:
#' - m = 1: Fox model, P(B) = r * B * log(K/B)
#' - m = 2: Schaefer model, P(B) = r * B * (1 - B/K)
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

  # Standard Pella-Tomlinson production function (with Fox limit at m = 1)
  production <- .pt_production(B, r, K, m)

  return(production)
}
