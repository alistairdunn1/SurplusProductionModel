#' Profile Likelihood Confidence Intervals
#'
#' Functions for computing profile likelihood confidence intervals for
#' model parameters and derived quantities (MSY, BMSY, FMSY).
#'
#' @name profiles
NULL

#' Profile Likelihood Confidence Intervals
#'
#' Computes profile likelihood confidence intervals by fixing each target
#' parameter at a grid of values and re-optimising the remaining parameters.
#' The profile log-likelihood ratio statistic
#' \eqn{2 (\ell_p(\theta) - \ell(\hat\theta))} is compared to a
#' \eqn{\chi^2_1} critical value to obtain the CI.
#'
#' @param model_fit A fitted \code{ProductionModel} object.
#' @param parameters Character vector of parameter names to profile.
#'   Supported names: \code{"r"}, \code{"K"}, \code{"m"}, \code{"sigma_obs"},
#'   \code{"sigma_proc"}, \code{"MSY"}, \code{"BMSY"}, \code{"FMSY"}, and
#'   per-area parameters such as \code{"q.A1"}, \code{"B0.A1"}.
#'   Defaults to \code{c("r", "K", "MSY", "BMSY")}.
#' @param ci_level Numeric, confidence level (default: 0.95).
#' @param n_points Integer, number of grid points per side of the MLE
#'   (total grid = \code{2 * n_points + 1}).  Default: 20.
#' @param range_factor Numeric, how many approximate standard errors away
#'   from the MLE to extend the profiling grid (default: 3).
#' @param delta Numeric, fractional step size used when SEs are unavailable
#'   (default: 0.5, meaning +/- 50 percent of the MLE).
#' @param biomass_target Optional numeric vector of target biomass fractions
#'   used to add depletion-based derived quantities such as `"B_40%K"`
#'   or `"F_40%B0"`.
#' @param baseline Character string indicating which biomass baseline to use
#'   for `biomass_target`: `"auto"` (default), `"B0"`, or `"K"`.
#' @param verbose Logical, print progress messages (default: FALSE).
#'
#' @return A list with class \code{"profile_likelihood"} containing:
#'   \describe{
#'     \item{profiles}{Named list.  Each element is a data.frame with columns
#'       \code{value} (natural-scale parameter value), \code{nll} (negative
#'       log-likelihood), \code{delta_nll} (profile deviance,
#'       \eqn{2(\text{nll} - \text{nll}_\text{min})}), and \code{converged}
#'       (logical).}
#'     \item{ci}{Named list of length-2 numeric vectors giving the profile
#'       CI bounds at \code{ci_level}.}
#'     \item{mle}{Named numeric vector of MLE values for each profiled
#'       parameter.}
#'     \item{ci_level}{Confidence level used.}
#'     \item{critical_value}{Chi-squared critical value used.}
#'   }
#'
#' @details
#' For model parameters (\code{r}, \code{K}, etc.) profiling is performed
#' on the log scale: the corresponding \code{log_*} parameter in the RTMB
#' objective is fixed at each grid value while all other parameters are
#' re-optimised.
#'
#' For derived quantities (\code{MSY}, \code{BMSY}, \code{FMSY}, and any
#' requested depletion-based targets), the function profiles over \code{r}
#' and \code{K}, records the derived quantity at each point, and then
#' inverts the profile deviance to find the CI.
#'
#' When standard errors from \code{sdreport} are available they are used
#' to set the grid range; otherwise the grid spans
#' \code{MLE * (1 +/- delta)}.
#'
#' @examples
#' \dontrun{
#' prof <- profile_likelihood(fitted_model)
#' print(prof)
#' plot(prof)
#'
#' # Profile specific parameters
#' prof_rK <- profile_likelihood(fitted_model,
#'                                parameters = c("r", "K"),
#'                                n_points = 30)
#'
#' # Depletion-based derived quantities (e.g., B_40%K and F_40%K)
#' prof_targets <- profile_likelihood(
#'   fitted_model,
#'   parameters = c("B_40%K", "F_40%K"),
#'   biomass_target = 0.4,
#'   baseline = "K"
#' )
#' }
#'
#' @export
profile_likelihood <- function(model_fit,
                               parameters = c("r", "K", "MSY", "BMSY"),
                               ci_level = 0.95,
                               n_points = 20L,
                               range_factor = 3,
                               delta = 0.5,
                               biomass_target = NULL,
                               baseline = c("auto", "B0", "K"),
                               verbose = FALSE) {
  defaults <- resolve_reference_point_defaults(
    biomass_target = biomass_target,
    baseline = baseline,
    biomass_target_missing = missing(biomass_target),
    baseline_missing = missing(baseline)
  )
  biomass_target <- defaults$biomass_target
  baseline <- defaults$baseline

  if (!inherits(model_fit, "ProductionModel")) {
    stop("Input must be a fitted ProductionModel object")
  }
  if (!model_fit$fitted) stop("Model must be fitted before profiling")

  obj <- model_fit$results$rtmb_obj
  if (is.null(obj)) {
    stop("rtmb_obj not found in results -- was the model fitted with RTMB?")
  }

  mle_nll   <- model_fit$results$likelihood
  mle_par   <- model_fit$results$opt_par   # optimised log-scale parameters
  if (is.null(mle_par)) {
    stop("Optimised parameters (opt_par) not found. Re-fit the model with the latest package version.")
  }
  se_vec    <- model_fit$results$std_errors
  nat_parms <- model_fit$parameters         # natural-scale named vector
  crit_val  <- qchisq(ci_level, df = 1)     # chi-sq(1) critical value
  ctrl      <- list(eval.max = 1000, iter.max = 500)

  # Derived quantity names that are not direct model parameters
  derived_names <- c("MSY", "BMSY", "FMSY")
  if (!is.null(biomass_target)) {
    target_ref <- calculate_reference_points_from_parameters(
      parameters = nat_parms,
      biomass_target = biomass_target,
      baseline = baseline
    )
    derived_names <- c(
      derived_names,
      target_ref$target_reference_points$biomass_name,
      target_ref$target_reference_points$f_name
    )
  }

  profiles <- list()
  ci_list  <- list()
  mle_list <- list()

  for (pname in parameters) {
    if (verbose) message("Profiling: ", pname)

    is_derived <- pname %in% derived_names || toupper(pname) %in% c("MSY", "BMSY", "FMSY")

    if (is_derived) {
      res <- .profile_derived(pname, model_fit, obj, mle_nll, mle_par,
                              se_vec, nat_parms, n_points, range_factor,
                              delta, crit_val, ctrl, verbose,
                              biomass_target, baseline)
    } else {
      res <- .profile_parameter(pname, model_fit, obj, mle_nll, mle_par,
                                se_vec, nat_parms, n_points, range_factor,
                                delta, crit_val, ctrl, verbose)
    }

    profiles[[pname]] <- res$profile_df
    ci_list[[pname]]  <- res$ci
    mle_list[[pname]] <- res$mle_val
  }

  out <- list(
    profiles       = profiles,
    ci             = ci_list,
    mle            = unlist(mle_list),
    ci_level       = ci_level,
    critical_value = crit_val
  )
  class(out) <- "profile_likelihood"
  out
}


# ---- Internal: profile a model parameter (log-scale) -------------------
.profile_parameter <- function(pname, model_fit, obj, mle_nll, mle_par,
                               se_vec, nat_parms, n_points, range_factor,
                               delta, crit_val, ctrl, verbose) {
  # Map natural-scale name to log-scale name in obj$par
  log_name <- .resolve_log_name(pname, names(mle_par), model_fit$data$areas)

  if (is.null(log_name)) {
    stop("Parameter '", pname, "' not found in the optimised parameter vector. ",
         "Available: ", paste(names(mle_par), collapse = ", "))
  }

  mle_log_val <- mle_par[[log_name]]
  mle_nat_val <- exp(mle_log_val)

  # Determine grid range on log scale
  se_log <- NULL
  if (!is.null(se_vec) && log_name %in% names(se_vec)) {
    se_cand <- se_vec[[log_name]]
    if (is.finite(se_cand) && se_cand > 0) se_log <- se_cand
  }

  if (!is.null(se_log)) {
    log_lo <- mle_log_val - range_factor * se_log
    log_hi <- mle_log_val + range_factor * se_log
  } else {
    # Use additive offset on log scale (equivalent to multiplicative on natural)
    log_lo <- mle_log_val - delta
    log_hi <- mle_log_val + delta
  }

  grid_log <- seq(log_lo, log_hi, length.out = 2 * n_points + 1)
  # Ensure MLE point is in the grid
  grid_log <- sort(unique(c(grid_log, mle_log_val)))

  idx_prof <- which(names(mle_par) == log_name)

  # Profile: fix one parameter, optimise the rest
  nll_vals  <- numeric(length(grid_log))
  conv_vals <- logical(length(grid_log))

  for (i in seq_along(grid_log)) {
    fixed_val <- grid_log[i]
    res <- .profile_one_point(obj, mle_par, idx_prof, fixed_val, ctrl)
    nll_vals[i]  <- res$nll
    conv_vals[i] <- res$converged
  }

  nat_vals  <- exp(grid_log)
  delta_nll <- 2 * (nll_vals - mle_nll)
  delta_nll[delta_nll < 0] <- 0  # numerical noise

  profile_df <- data.frame(
    value     = nat_vals,
    nll       = nll_vals,
    delta_nll = delta_nll,
    converged = conv_vals,
    stringsAsFactors = FALSE
  )

  ci <- .find_ci_from_profile(nat_vals, delta_nll, crit_val, mle_nat_val)

  list(profile_df = profile_df, ci = ci, mle_val = mle_nat_val)
}


# ---- Internal: profile a derived quantity (MSY, BMSY, FMSY) -----------
.profile_derived <- function(pname, model_fit, obj, mle_nll, mle_par,
                             se_vec, nat_parms, n_points, range_factor,
                             delta, crit_val, ctrl, verbose,
                             biomass_target, baseline) {
  # Strategy: profile over all estimable parameters one at a time
  # and compute the derived quantity at each profile point.
  # We profile r (the most influential on MSY/BMSY) unless it's fixed.

  # Compute MLE-derived quantities
  ref_mle <- calculate_reference_points_from_parameters(
    parameters = nat_parms,
    biomass_target = biomass_target,
    baseline = baseline,
    warn_on_invalid_m = FALSE
  )

  target <- if (toupper(pname) %in% c("MSY", "BMSY", "FMSY")) toupper(pname) else pname
  mle_dq <- extract_named_reference_value(ref_mle, target)

  # Profile over r to get the derived-quantity profile

  log_r_name <- .resolve_log_name("r", names(mle_par), model_fit$data$areas)
  if (is.null(log_r_name)) {
    stop("Cannot profile '", pname, "': parameter 'r' not found in optimised parameters")
  }
  mle_log_r <- mle_par[[log_r_name]]
  idx_r <- which(names(mle_par) == log_r_name)

  # Grid range for r
  se_log_r <- NULL
  if (!is.null(se_vec) && log_r_name %in% names(se_vec)) {
    se_cand <- se_vec[[log_r_name]]
    if (is.finite(se_cand) && se_cand > 0) se_log_r <- se_cand
  }

  if (!is.null(se_log_r)) {
    log_lo <- mle_log_r - range_factor * se_log_r
    log_hi <- mle_log_r + range_factor * se_log_r
  } else {
    log_lo <- mle_log_r - delta
    log_hi <- mle_log_r + delta
  }

  grid_log_r <- seq(log_lo, log_hi, length.out = 2 * n_points + 1)
  grid_log_r <- sort(unique(c(grid_log_r, mle_log_r)))

  # Also profile over K for better derived-quantity envelope
  log_K_name <- .resolve_log_name("K", names(mle_par), model_fit$data$areas)
  if (!is.null(log_K_name)) {
    mle_log_K <- mle_par[[log_K_name]]
    idx_K <- which(names(mle_par) == log_K_name)
    se_log_K <- NULL
    if (!is.null(se_vec) && log_K_name %in% names(se_vec)) {
      se_cand <- se_vec[[log_K_name]]
      if (is.finite(se_cand) && se_cand > 0) se_log_K <- se_cand
    }

    if (!is.null(se_log_K)) {
      log_K_lo <- mle_log_K - range_factor * se_log_K
      log_K_hi <- mle_log_K + range_factor * se_log_K
    } else {
      log_K_lo <- mle_log_K - delta
      log_K_hi <- mle_log_K + delta
    }
    grid_log_K <- seq(log_K_lo, log_K_hi, length.out = 2 * n_points + 1)
    grid_log_K <- sort(unique(c(grid_log_K, mle_log_K)))
  }

  # Collect from profiling r
  dq_vals_r   <- numeric(length(grid_log_r))
  nll_vals_r  <- numeric(length(grid_log_r))
  conv_vals_r <- logical(length(grid_log_r))

  for (i in seq_along(grid_log_r)) {
    res <- .profile_one_point(obj, mle_par, idx_r, grid_log_r[i], ctrl)
    nll_vals_r[i]  <- res$nll
    conv_vals_r[i] <- res$converged
    # Reconstruct natural-scale pars at this profile point
    nat_i <- .naturalise_profile_parameters(res$par_full)
    ref_i <- calculate_reference_points_from_parameters(
      parameters = nat_i,
      biomass_target = biomass_target,
      baseline = baseline,
      warn_on_invalid_m = FALSE
    )
    dq_vals_r[i] <- extract_named_reference_value(ref_i, target)
  }

  # Collect from profiling K (if available)
  dq_vals   <- dq_vals_r
  nll_vals  <- nll_vals_r
  conv_vals <- conv_vals_r

  if (!is.null(log_K_name)) {
    dq_vals_K   <- numeric(length(grid_log_K))
    nll_vals_K  <- numeric(length(grid_log_K))
    conv_vals_K <- logical(length(grid_log_K))

    for (i in seq_along(grid_log_K)) {
      res <- .profile_one_point(obj, mle_par, idx_K, grid_log_K[i], ctrl)
      nll_vals_K[i]  <- res$nll
      conv_vals_K[i] <- res$converged
      nat_i <- .naturalise_profile_parameters(res$par_full)
      ref_i <- calculate_reference_points_from_parameters(
        parameters = nat_i,
        biomass_target = biomass_target,
        baseline = baseline,
        warn_on_invalid_m = FALSE
      )
      dq_vals_K[i] <- extract_named_reference_value(ref_i, target)
    }

    dq_vals   <- c(dq_vals, dq_vals_K)
    nll_vals  <- c(nll_vals, nll_vals_K)
    conv_vals <- c(conv_vals, conv_vals_K)
  }

  delta_nll <- 2 * (nll_vals - mle_nll)
  delta_nll[delta_nll < 0] <- 0

  # Sort by derived quantity value for clean profiles
  ord <- order(dq_vals)
  profile_df <- data.frame(
    value     = dq_vals[ord],
    nll       = nll_vals[ord],
    delta_nll = delta_nll[ord],
    converged = conv_vals[ord],
    stringsAsFactors = FALSE
  )

  ci <- .find_ci_from_profile(dq_vals[ord], delta_nll[ord], crit_val, mle_dq)

  list(profile_df = profile_df, ci = ci, mle_val = mle_dq)
}


# ---- Internal: fix one parameter and re-optimise the rest --------------
.profile_one_point <- function(obj, mle_par, idx_fixed, fixed_val, ctrl) {
  n_par <- length(mle_par)
  bad_nll <- .Machine$double.xmax / 1e100
  # Indices of free parameters (all except the fixed one)
  idx_free <- setdiff(seq_len(n_par), idx_fixed)

  # Starting values for free parameters from the MLE
  start_free <- mle_par[idx_free]

  # Wrapper: insert fixed value, call obj$fn
  fn_wrapper <- function(par_free) {
    full_par <- numeric(n_par)
    full_par[idx_fixed] <- fixed_val
    full_par[idx_free]  <- par_free
    names(full_par) <- names(mle_par)
    val <- obj$fn(full_par)
    if (!is.finite(val)) {
      return(bad_nll)
    }
    val
  }

  gr_wrapper <- function(par_free) {
    full_par <- numeric(n_par)
    full_par[idx_fixed] <- fixed_val
    full_par[idx_free]  <- par_free
    names(full_par) <- names(mle_par)
    full_gr <- obj$gr(full_par)
    free_gr <- full_gr[idx_free]
    free_gr[!is.finite(free_gr)] <- 0
    free_gr
  }

  opt <- tryCatch(
    withCallingHandlers(
      nlminb(start_free, fn_wrapper, gr_wrapper, control = ctrl),
      warning = function(w) {
        if (grepl("NA/NaN function evaluation", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) list(objective = NA_real_, convergence = 99L,
                             par = start_free)
  )

  if (!is.finite(opt$objective)) {
    opt$objective <- bad_nll
    opt$convergence <- 98L
  }

  # Reconstruct full parameter vector
  full_par <- numeric(n_par)
  full_par[idx_fixed] <- fixed_val
  full_par[idx_free]  <- opt$par
  names(full_par) <- names(mle_par)

  list(
    nll       = opt$objective,
    converged = identical(opt$convergence, 0L),
    par_full  = full_par
  )
}


# ---- Internal: find CI bounds from profile deviance --------------------
.find_ci_from_profile <- function(values, delta_nll, crit_val, mle_val) {
  # Linearly interpolate the profile deviance to find where it crosses crit_val
  # Lower bound: search to the left of MLE
  # Upper bound: search to the right of MLE

  ci <- c(lower = NA_real_, upper = NA_real_)

  # Only use converged points (delta_nll should be finite)
  ok <- is.finite(delta_nll) & is.finite(values)
  if (sum(ok) < 3) return(ci)

  vals <- values[ok]
  dnll <- delta_nll[ok]

  # Find lower bound
  left <- which(vals <= mle_val)
  if (length(left) >= 2) {
    # Walk from MLE leftward
    left_vals <- vals[left]
    left_dnll <- dnll[left]
    # Re-order from MLE outward (decreasing value)
    ord_l <- order(left_vals, decreasing = TRUE)
    left_vals <- left_vals[ord_l]
    left_dnll <- left_dnll[ord_l]

    cross_idx <- which(left_dnll >= crit_val)
    if (length(cross_idx) > 0) {
      j <- cross_idx[1]
      if (j >= 2) {
        # Linear interpolation between points j-1 and j
        ci["lower"] <- .interp_crossing(left_vals[j - 1], left_dnll[j - 1],
                                        left_vals[j], left_dnll[j], crit_val)
      } else {
        ci["lower"] <- left_vals[j]
      }
    }
  }

  # Find upper bound
  right <- which(vals >= mle_val)
  if (length(right) >= 2) {
    right_vals <- vals[right]
    right_dnll <- dnll[right]
    ord_r <- order(right_vals)
    right_vals <- right_vals[ord_r]
    right_dnll <- right_dnll[ord_r]

    cross_idx <- which(right_dnll >= crit_val)
    if (length(cross_idx) > 0) {
      j <- cross_idx[1]
      if (j >= 2) {
        ci["upper"] <- .interp_crossing(right_vals[j - 1], right_dnll[j - 1],
                                        right_vals[j], right_dnll[j], crit_val)
      } else {
        ci["upper"] <- right_vals[j]
      }
    }
  }

  ci
}


# ---- Internal: linear interpolation to find crossing point -------------
.interp_crossing <- function(x1, y1, x2, y2, target_y) {
  if (abs(y2 - y1) < 1e-12) return((x1 + x2) / 2)
  x1 + (target_y - y1) * (x2 - x1) / (y2 - y1)
}


# ---- Internal: resolve natural-scale name to log-scale name -----------
.resolve_log_name <- function(pname, par_names, areas = NULL) {
  # Direct log-scale mapping
  # e.g. "r" -> "log_r", "K" -> "log_K", "sigma_obs" -> "log_sigma_obs"
  # Per-area: "q.A1" -> "log_q_A1", "B0.A1" -> "log_B0_A1"

  # Try straightforward conversion
  log_name <- paste0("log_", gsub("\\.", "_", pname))
  if (log_name %in% par_names) return(log_name)

  # For single-area models: "q" might map to "log_q_A1"
  if (!is.null(areas) && length(areas) == 1) {
    log_name_suffixed <- paste0("log_", gsub("\\.", "_", pname), "_", areas[1])
    if (log_name_suffixed %in% par_names) return(log_name_suffixed)
  }

  # Fuzzy match: try partial match
  candidates <- grep(paste0("^log_", gsub("\\.", "_", pname)), par_names, value = TRUE)
  if (length(candidates) == 1) return(candidates)

  NULL
}


# ---- Internal: naturalise profile parameters ---------------------------
.naturalise_profile_parameters <- function(log_par) {
  nat_par <- exp(log_par)
  nat_names <- sub("^log_", "", names(log_par))
  nat_names <- sub("^(q|B0|sigma_proc|sigma_obs)_([A-Z])", "\\1.\\2", nat_names)
  names(nat_par) <- nat_names
  nat_par
}


#' @export
print.profile_likelihood <- function(x, ...) {
  cat("Profile Likelihood Confidence Intervals\n")
  cat("========================================\n")
  cat("Confidence level:", sprintf("%.0f%%", x$ci_level * 100), "\n")
  cat("Critical value (chi-sq):", sprintf("%.3f", x$critical_value), "\n\n")

  for (pname in names(x$profiles)) {
    ci <- x$ci[[pname]]
    mle <- x$mle[[pname]]
    cat(sprintf("  %-12s  MLE = %10.4g   CI = [%10.4g, %10.4g]\n",
                pname, mle,
                ifelse(is.na(ci["lower"]), NA, ci["lower"]),
                ifelse(is.na(ci["upper"]), NA, ci["upper"])))
  }
  cat("\n")
  invisible(x)
}


#' Plot Profile Likelihood
#'
#' Plots the profile deviance curves for each profiled parameter with
#' the chi-squared critical value (CI threshold) overlaid.
#'
#' @param x A \code{profile_likelihood} object.
#' @param ... Additional arguments (not used).
#' @return A ggplot object (invisibly).
#' @export
plot.profile_likelihood <- function(x, ...) {
  n_prof <- length(x$profiles)
  if (n_prof == 0) {
    message("No profiles to plot.")
    return(invisible(NULL))
  }

  # Combine into long data.frame
  dfs <- list()
  for (pname in names(x$profiles)) {
    df_i <- x$profiles[[pname]]
    df_i$parameter <- pname
    dfs[[length(dfs) + 1]] <- df_i
  }
  df <- do.call(rbind, dfs)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = value, y = delta_nll)) +
    ggplot2::geom_line(colour = "#0072B2", linewidth = 0.8) +
    ggplot2::geom_point(colour = "#0072B2", size = 1.2, alpha = 0.6) +
    ggplot2::geom_hline(yintercept = x$critical_value,
                         linetype = "dashed", colour = "#D55E00") +
    ggplot2::annotate("text", x = -Inf, y = x$critical_value,
                       label = sprintf("  %.0f%% CI", x$ci_level * 100),
                       vjust = -0.5, hjust = 0,
                       colour = "#D55E00", size = 3) +
    ggplot2::facet_wrap(~parameter, scales = "free_x") +
    ggplot2::labs(title = "Profile Likelihood",
                  x = "Parameter Value",
                  y = expression(2 * Delta * NLL)) +
    ggplot2::theme_bw()

  # Add vertical lines for MLE and CI bounds
  for (pname in names(x$profiles)) {
    mle_val <- x$mle[[pname]]
    ci_vals <- x$ci[[pname]]
    p <- p + ggplot2::geom_vline(
      data = data.frame(parameter = pname, xint = mle_val),
      ggplot2::aes(xintercept = xint),
      linetype = "solid", colour = "grey50", linewidth = 0.4
    )
    if (!is.na(ci_vals["lower"])) {
      p <- p + ggplot2::geom_vline(
        data = data.frame(parameter = pname, xint = ci_vals["lower"]),
        ggplot2::aes(xintercept = xint),
        linetype = "dotted", colour = "#009E73", linewidth = 0.5
      )
    }
    if (!is.na(ci_vals["upper"])) {
      p <- p + ggplot2::geom_vline(
        data = data.frame(parameter = pname, xint = ci_vals["upper"]),
        ggplot2::aes(xintercept = xint),
        linetype = "dotted", colour = "#009E73", linewidth = 0.5
      )
    }
  }

  print(p)
  invisible(p)
}


# Declare global variables for R CMD check NSE
utils::globalVariables(c("value", "delta_nll", "parameter", "xint"))
