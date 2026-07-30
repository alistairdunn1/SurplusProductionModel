#' Forward Projection Under Assumed Harvest Control
#'
#' Projects a fitted Pella-Tomlinson model forward under assumed harvest
#' control, specified either as target fishing mortality
#' (\code{control_type = "F"}) or catch (\code{control_type = "catch"}).
#' Process uncertainty can be propagated by sampling from either a normal
#' process-error distribution or empirical historical process deviations.
#'
#' @param model_fit A fitted \code{ProductionModel} object.
#' @param horizon Integer projection length in years (default \code{10}).
#' @param control Numeric scalar or vector of length \code{horizon}. Interpreted
#'   as fishing mortality (\code{control_type = "F"}) or total catch
#'   (\code{control_type = "catch"}).
#' @param control_type Character, one of \code{"F"} or \code{"catch"}.
#' @param n_sim Number of Monte Carlo simulation trajectories (default
#'   \code{1000}).
#' @param process_error Character, one of \code{"historical"}, \code{"normal"},
#'   or \code{"none"}. Default is \code{"historical"}.
#' @param historical_window Number of most recent historical years to use when
#'   \code{process_error = "historical"} (default \code{10}).
#' @param bias_correction Logical. When \code{TRUE} (default) and
#'   \code{process_error = "normal"}, the lognormal process deviations include
#'   the \eqn{-\sigma^2/2} correction so that the projected biomass is
#'   mean-unbiased. Ignored for \code{"historical"} (empirical) and
#'   \code{"none"}.
#' @param seed Optional integer random seed for reproducible simulation.
#' @param probs Numeric vector of quantile probabilities used for summarized
#'   projection output.
#'
#' @return An object of class \code{"pt_projection"} with elements:
#'   \describe{
#'     \item{summary}{Data frame with yearly quantiles and status probabilities.}
#'     \item{years}{Projected years.}
#'     \item{control}{Expanded control vector used in projection.}
#'     \item{control_type}{Control interpretation (\code{"F"} or \code{"catch"}).}
#'     \item{simulations}{List containing simulated biomass, catch, and harvest-rate arrays.}
#'     \item{settings}{List of simulation settings used.}
#'   }
#'
#' @details
#' Historical process deviations are reconstructed from fitted biomass
#' trajectories as
#' \deqn{\epsilon_t = \log(B_{t+1}) - \log(\max(B_t + P(B_t) - C_t, \epsilon))}
#' and sampled with replacement from the most recent
#' \code{historical_window} years (default 10).
#'
#' @examples
#' \dontrun{
#' fit <- fit_pella_tomlinson_model(data_list)
#'
#' # 10-year projection under constant F with historical process error resampling
#' proj <- project_forward(
#'   fit,
#'   horizon = 10,
#'   control = 0.12,
#'   control_type = "F",
#'   process_error = "historical",
#'   historical_window = 10
#' )
#'
#' print(proj)
#' head(proj$summary)
#' }
#'
#' @export
project_forward <- function(model_fit,
                            horizon = 10L,
                            control,
                            control_type = c("F", "catch"),
                            n_sim = 1000L,
                            process_error = c("historical", "normal", "none"),
                            historical_window = 10L,
                            bias_correction = TRUE,
                            seed = NULL,
                            probs = c(0.05, 0.5, 0.95)) {
  if (!inherits(model_fit, "ProductionModel")) {
    stop("model_fit must be a ProductionModel object")
  }
  if (!isTRUE(model_fit$fitted)) {
    stop("model_fit must be fitted before projection")
  }

  control_type <- match.arg(control_type)
  process_error <- match.arg(process_error)

  horizon <- as.integer(horizon)
  n_sim <- as.integer(n_sim)
  historical_window <- as.integer(historical_window)

  if (!is.finite(horizon) || horizon < 1L) stop("horizon must be >= 1")
  if (!is.finite(n_sim) || n_sim < 1L) stop("n_sim must be >= 1")
  if (!is.finite(historical_window) || historical_window < 1L) {
    stop("historical_window must be >= 1")
  }
  if (!is.numeric(control) || any(!is.finite(control)) || any(control < 0)) {
    stop("control must be a non-negative numeric scalar or vector")
  }
  if (!(length(control) == 1L || length(control) == horizon)) {
    stop("control must be length 1 or length horizon")
  }
  if (!is.numeric(probs) || any(!is.finite(probs)) || any(probs <= 0 | probs >= 1)) {
    stop("probs must be numeric values strictly between 0 and 1")
  }
  if (!is.logical(bias_correction) || length(bias_correction) != 1L ||
    is.na(bias_correction)) {
    stop("bias_correction must be a single logical value")
  }

  if (!is.null(seed)) set.seed(seed)

  control_vec <- if (length(control) == 1L) rep(as.numeric(control), horizon) else as.numeric(control)

  # Core fitted values
  pars <- model_fit$parameters
  r <- unname(pars[["r"]])
  K <- unname(pars[["K"]])
  m <- unname(pars[["m"]])

  if (any(!is.finite(c(r, K, m))) || any(c(r, K, m) <= 0)) {
    stop("model_fit parameters must include positive finite r, K, and m")
  }

  b_hist <- model_fit$results$biomass
  c_hist <- model_fit$data$catch

  # Normalize historical biomass/catch to matrices [year x area]
  b_hist_mat <- if (is.matrix(b_hist)) b_hist else matrix(as.numeric(b_hist), ncol = 1)
  c_hist_mat <- if (is.matrix(c_hist)) c_hist else matrix(as.numeric(c_hist), ncol = 1)

  n_years <- nrow(b_hist_mat)
  n_areas <- ncol(b_hist_mat)
  if (n_years < 2L) stop("Need at least 2 historical years for projection")

  # Ensure catch history aligns in rows
  if (nrow(c_hist_mat) != n_years) {
    stop("Historical catch and biomass lengths do not align")
  }

  # Area shares for distributing total catch in multi-area models
  last_catch <- c_hist_mat[n_years, ]
  if (sum(last_catch, na.rm = TRUE) > 0) {
    area_share <- as.numeric(last_catch / sum(last_catch, na.rm = TRUE))
  } else {
    area_share <- rep(1 / n_areas, n_areas)
  }

  # Historical process deviation reconstruction by area
  hist_eps <- matrix(NA_real_, nrow = n_years - 1L, ncol = n_areas)
  for (t in seq_len(n_years - 1L)) {
    bt <- b_hist_mat[t, ]
    prod_t <- .pt_production(bt, r, K, m)
    prod_t[!is.finite(prod_t)] <- 0
    b_det_next <- pmax(bt + prod_t - c_hist_mat[t, ], 1e-8)
    b_next <- pmax(b_hist_mat[t + 1L, ], 1e-8)
    hist_eps[t, ] <- log(b_next) - log(b_det_next)
  }

  # Restrict to recent historical window
  use_rows <- seq.int(max(1L, nrow(hist_eps) - historical_window + 1L), nrow(hist_eps))
  hist_eps_recent <- hist_eps[use_rows, , drop = FALSE]
  good_rows <- apply(hist_eps_recent, 1, function(x) all(is.finite(x)))
  hist_eps_recent <- hist_eps_recent[good_rows, , drop = FALSE]

  sigma_proc <- pars[["sigma_proc"]]
  if (is.null(sigma_proc) || !is.finite(sigma_proc) || sigma_proc <= 0) {
    sigma_proc <- stats::sd(as.numeric(hist_eps_recent), na.rm = TRUE)
    if (!is.finite(sigma_proc) || sigma_proc <= 0) sigma_proc <- 0.1
  }

  if (process_error == "historical" && nrow(hist_eps_recent) == 0L) {
    warning("No finite historical process deviations available; falling back to normal process error")
    process_error <- "normal"
  }

  # Simulated arrays
  b_sim <- array(NA_real_, dim = c(horizon + 1L, n_areas, n_sim))
  c_sim <- array(NA_real_, dim = c(horizon, n_areas, n_sim))
  f_sim <- array(NA_real_, dim = c(horizon, n_areas, n_sim))

  b0 <- pmax(b_hist_mat[n_years, ], 1e-8)
  for (s in seq_len(n_sim)) {
    b_sim[1L, , s] <- b0

    for (h in seq_len(horizon)) {
      bt <- pmax(b_sim[h, , s], 1e-8)
      prod_h <- .pt_production(bt, r, K, m)
      prod_h[!is.finite(prod_h)] <- 0

      if (control_type == "F") {
        f_target <- control_vec[h]
        catch_h <- f_target * bt
      } else {
        catch_h <- control_vec[h] * area_share
      }

      # Prevent impossible removals
      catch_cap <- pmax(bt + prod_h - 0.01, 0)
      catch_h <- pmin(catch_h, catch_cap)

      b_det_next <- pmax(bt + prod_h - catch_h, 1e-8)

      # Normal deviations are drawn with the -sigma^2/2 lognormal bias
      # correction (when bias_correction = TRUE) so that E[exp(eps)] = 1 and
      # the projected biomass is mean-unbiased, consistent with the operating
      # model. Historical deviations are empirical realisations and are used
      # as-is (their sample mean already reflects the fitted dynamics).
      normal_mean <- if (isTRUE(bias_correction)) -0.5 * sigma_proc^2 else 0
      eps_h <- switch(process_error,
        none = rep(0, n_areas),
        normal = stats::rnorm(n_areas, mean = normal_mean, sd = sigma_proc),
        historical = {
          draw <- sample.int(nrow(hist_eps_recent), size = 1L)
          as.numeric(hist_eps_recent[draw, ])
        }
      )

      b_next <- pmax(b_det_next * exp(eps_h), 0.01)

      b_sim[h + 1L, , s] <- b_next
      c_sim[h, , s] <- catch_h
      f_sim[h, , s] <- catch_h / pmax(bt, 1e-8)
    }
  }

  # Aggregate status metrics (total biomass/catch across areas)
  b_tot <- apply(b_sim, c(1, 3), sum)
  c_tot <- apply(c_sim, c(1, 3), sum)
  f_tot <- c_tot / pmax(b_tot[-(horizon + 1L), , drop = FALSE], 1e-8)

  # Build summary by future year
  proj_years <- seq.int(max(model_fit$data$years) + 1L, by = 1L, length.out = horizon)

  bmsy <- model_fit$results$bmsy
  fmsy <- model_fit$results$fmsy

  summary_df <- data.frame(
    year = proj_years,
    biomass_q05 = NA_real_,
    biomass_q50 = NA_real_,
    biomass_q95 = NA_real_,
    harvest_q05 = NA_real_,
    harvest_q50 = NA_real_,
    harvest_q95 = NA_real_,
    prob_B_above_BMSY = NA_real_,
    prob_F_below_FMSY = NA_real_,
    stringsAsFactors = FALSE
  )

  for (h in seq_len(horizon)) {
    b_h <- b_tot[h + 1L, ]
    f_h <- f_tot[h, ]
    bq <- stats::quantile(b_h, probs = c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE)
    fq <- stats::quantile(f_h, probs = c(0.05, 0.5, 0.95), na.rm = TRUE, names = FALSE)

    summary_df$biomass_q05[h] <- bq[1]
    summary_df$biomass_q50[h] <- bq[2]
    summary_df$biomass_q95[h] <- bq[3]
    summary_df$harvest_q05[h] <- fq[1]
    summary_df$harvest_q50[h] <- fq[2]
    summary_df$harvest_q95[h] <- fq[3]

    summary_df$prob_B_above_BMSY[h] <- if (is.finite(bmsy) && bmsy > 0) mean(b_h >= bmsy) else NA_real_
    summary_df$prob_F_below_FMSY[h] <- if (is.finite(fmsy) && fmsy > 0) mean(f_h <= fmsy) else NA_real_
  }

  out <- list(
    summary = summary_df,
    years = proj_years,
    control = control_vec,
    control_type = control_type,
    simulations = list(
      biomass = b_sim,
      catch = c_sim,
      harvest_rate = f_sim
    ),
    settings = list(
      horizon = horizon,
      n_sim = n_sim,
      process_error = process_error,
      historical_window = historical_window,
      bias_correction = bias_correction,
      sigma_proc = sigma_proc,
      probs = probs,
      seed = seed
    )
  )

  class(out) <- "pt_projection"
  out
}

#' @export
print.pt_projection <- function(x, ...) {
  cat("Forward Projection (Pella-Tomlinson)\n")
  cat("===================================\n")
  cat("Horizon:", x$settings$horizon, "years\n")
  cat("Simulations:", x$settings$n_sim, "\n")
  cat("Control type:", x$control_type, "\n")
  cat("Process error:", x$settings$process_error, "\n\n")

  last_row <- x$summary[nrow(x$summary), , drop = FALSE]
  cat("Terminal Year:", last_row$year, "\n")
  cat("Biomass median (q05-q95):",
    sprintf("%.2f (%.2f-%.2f)",
      last_row$biomass_q50, last_row$biomass_q05, last_row$biomass_q95
    ), "\n"
  )
  cat("Harvest median (q05-q95):",
    sprintf("%.3f (%.3f-%.3f)",
      last_row$harvest_q50, last_row$harvest_q05, last_row$harvest_q95
    ), "\n"
  )

  if (is.finite(last_row$prob_B_above_BMSY)) {
    cat("Pr(B >= BMSY):", sprintf("%.3f", last_row$prob_B_above_BMSY), "\n")
  }
  if (is.finite(last_row$prob_F_below_FMSY)) {
    cat("Pr(F <= FMSY):", sprintf("%.3f", last_row$prob_F_below_FMSY), "\n")
  }

  invisible(x)
}
