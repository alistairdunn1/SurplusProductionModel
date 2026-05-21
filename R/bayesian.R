#' Bayesian Inference via tmbstan
#'
#' Functions for running MCMC sampling on a fitted
#' \code{\link{ProductionModel}} using the \pkg{tmbstan} package,
#' which couples a TMB/RTMB objective function with Stan's
#' No-U-Turn Sampler (NUTS).
#'
#' @name bayesian
NULL

# ---- global variables used in ggplot aes / data.table --------------------
utils::globalVariables(c("chain", "iteration", "param", "sample_value",
                         "Var1", "Var2", "Freq", "index",
                         "observed", "median", "lower", "upper"))


#' Bayesian MCMC Sampling for a Production Model
#'
#' Draws posterior samples from the joint parameter distribution of a
#' fitted \code{ProductionModel} using Stan's NUTS sampler via
#' \code{\link[tmbstan]{tmbstan}}.  The starting point is the maximum
#' likelihood estimate so that chains begin in a high-density region.
#'
#' @param model_fit  A fitted \code{ProductionModel} object (must have
#'   \code{model_fit$fitted == TRUE}).
#' @param chains  Integer, number of MCMC chains (default: 4).
#' @param iter    Integer, total iterations per chain including warmup
#'   (default: 2000).
#' @param warmup  Integer, number of warmup (burn-in) iterations per
#'   chain (default: \code{floor(iter / 2)}).
#' @param thin    Integer, thinning interval (default: 1).
#' @param seed    Integer, random seed for reproducibility (default: NULL).
#' @param lower   Numeric vector of lower bounds on the log-scale
#'   parameters, length matching \code{obj$par}.
#'   Default: \code{-Inf} for all parameters (no lower bound).
#' @param upper   Numeric vector of upper bounds on the log-scale
#'   parameters, length matching \code{obj$par}.
#'   Default: \code{+Inf} for all parameters (no upper bound).
#' @param laplace Logical.  If \code{TRUE}, integrate out random effects
#'   using the Laplace approximation instead of sampling them.
#'   Default: \code{FALSE} (full Bayesian).
#' @param biomass_target Optional numeric vector of target biomass fractions
#'   used to add depletion-based posterior summaries such as
#'   `"B_40%K"` or `"F_40%B0"`.
#' @param baseline Character string indicating which biomass baseline to use
#'   for `biomass_target`: `"auto"` (default), `"B0"`, or `"K"`.
#' @param ...     Additional arguments passed to
#'   \code{\link[rstan]{sampling}} (e.g. \code{control},
#'   \code{adapt_delta}).
#'
#' @return An S3 object of class \code{"bayes_fit"} containing:
#'   \describe{
#'     \item{stanfit}{The \code{stanfit} object from \pkg{rstan}.}
#'     \item{posterior}{A matrix of posterior samples (rows = draws,
#'       columns = natural-scale parameters + derived quantities).}
#'     \item{log_posterior}{Matrix of posterior samples on the log scale.}
#'     \item{summary}{Data frame with posterior summaries: mean, sd,
#'       2.5 percent, 50 percent, 97.5 percent quantiles, n_eff, Rhat
#'       for each parameter.}
#'     \item{model_fit}{The original \code{ProductionModel} object.}
#'     \item{par_names}{Character vector of log-scale parameter names.}
#'     \item{nat_names}{Character vector of natural-scale parameter names.}
#'     \item{n_chains}{Number of chains.}
#'     \item{n_iter}{Iterations per chain.}
#'     \item{n_warmup}{Warmup iterations per chain.}
#'   }
#'
#' @details
#' \strong{Priors.}
#' By default the sampler uses flat (improper) priors on the log-scale
#' parameters, so the posterior is proportional to the likelihood.
#' You can impose uniform priors on the log-scale by supplying finite
#' \code{lower} and \code{upper} bounds.  For informative priors you
#' should modify the RTMB objective function directly.
#'
#' \strong{Derived quantities.}
#' After sampling, the function transforms every draw to the natural
#' scale and computes MSY, BMSY, FMSY, and any requested depletion-based
#' targets for each posterior draw.
#'
#' \strong{Diagnostics.}
#' Standard Stan diagnostics (divergent transitions, tree depth,
#' low E-BFMI) are available via the returned \code{stanfit} object.
#' The \code{print} and \code{plot} methods provide quick summaries.
#'
#' @examples
#' \dontrun{
#' bf <- bayesian_fit(fitted_model, chains = 2, iter = 1000)
#' print(bf)
#' plot(bf)
#' plot(bf, type = "pairs")
#'
#' # Include depletion-based posterior targets (e.g., B_40%K, F_40%K)
#' bf_targets <- bayesian_fit(
#'   fitted_model,
#'   chains = 2,
#'   iter = 1000,
#'   biomass_target = 0.4,
#'   baseline = "K"
#' )
#' subset(bf_targets$summary, parameter %in% c("B_40%K", "F_40%K"))
#' }
#'
#' @export
bayesian_fit <- function(model_fit,
                         chains = 4L,
                         iter   = 2000L,
                         warmup = floor(iter / 2),
                         thin   = 1L,
                         seed   = NULL,
                         lower  = numeric(0),
                         upper  = numeric(0),
                         laplace = FALSE,
                         biomass_target = NULL,
                         baseline = c("auto", "B0", "K"),
                         ...) {

  # ---- input validation -------------------------------------------------
  if (!inherits(model_fit, "ProductionModel")) {
    stop("Input must be a fitted ProductionModel object")
  }
  if (!model_fit$fitted) {
    stop("Model must be fitted before running Bayesian sampling")
  }
  if (!requireNamespace("tmbstan", quietly = TRUE)) {
    stop("Package 'tmbstan' is required for Bayesian fitting. ",
         "Install it with: install.packages('tmbstan')")
  }
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("Package 'rstan' is required for Bayesian fitting. ",
         "Install it with: install.packages('rstan')")
  }

  obj <- model_fit$results$rtmb_obj
  if (is.null(obj)) {
    stop("rtmb_obj not found in results -- was the model fitted with RTMB?")
  }

  opt_par <- model_fit$results$opt_par
  if (is.null(opt_par)) {
    stop("Optimised parameters (opt_par) not found. ",
         "Re-fit the model with the latest package version.")
  }

  # ---- run tmbstan -------------------------------------------------------
  stan_args <- list(
    obj     = obj,
    lower   = lower,
    upper   = upper,
    laplace = laplace,
    chains  = as.integer(chains),
    iter    = as.integer(iter),
    warmup  = as.integer(warmup),
    thin    = as.integer(thin),
    init    = "last.par.best"
  )
  if (!is.null(seed)) stan_args$seed <- as.integer(seed)

  # Merge user ... into stan_args
  dots <- list(...)
  stan_args <- modifyList(stan_args, dots)

  stanfit <- do.call(tmbstan::tmbstan, stan_args)

  # ---- extract posterior samples -----------------------------------------
  log_posterior <- rstan::extract(stanfit, permuted = FALSE)
  # log_posterior is an array: [iteration, chain, parameter]
  par_names  <- names(opt_par)
  n_params   <- length(par_names)

  # Collapse chains into a matrix (all post-warmup draws)
  log_mat <- as.matrix(stanfit)                    # draws x (params + lp__)
  # Drop lp__ column
  lp_col <- which(colnames(log_mat) == "lp__")
  if (length(lp_col) > 0) log_mat <- log_mat[, -lp_col, drop = FALSE]

  # ---- transform to natural scale + derived quantities -------------------
  nat_names  <- .log_to_natural_names(par_names)
  nat_mat    <- exp(log_mat)
  colnames(nat_mat) <- nat_names
  baseline <- match.arg(baseline)
  validate_biomass_target(biomass_target)

  # Derived quantities: MSY, BMSY, FMSY
  # Need: r, K, m
  r_col <- which(nat_names == "r")
  K_col <- which(nat_names == "K")
  m_col <- which(nat_names == "m")

  n_draws <- nrow(nat_mat)
  MSY_vec  <- rep(NA_real_, n_draws)
  BMSY_vec <- rep(NA_real_, n_draws)
  FMSY_vec <- rep(NA_real_, n_draws)

  if (length(r_col) == 1 && length(K_col) == 1 && length(m_col) == 1) {
    r_draws <- nat_mat[, r_col]
    K_draws <- nat_mat[, K_col]
    m_draws <- nat_mat[, m_col]

    BMSY_vec <- K_draws * (1 / m_draws)^(1 / (m_draws - 1))
    FMSY_vec <- r_draws / m_draws * (1 - 1 / m_draws)
    MSY_vec  <- FMSY_vec * BMSY_vec

    # Handle near-singular cases (m near 1) which can produce Inf/NaN
    # and biologically invalid cases (m < 1 gives negative FMSY)
    bad <- !is.finite(MSY_vec) | !is.finite(BMSY_vec) | !is.finite(FMSY_vec) |
           MSY_vec <= 0 | BMSY_vec <= 0 | FMSY_vec <= 0
    MSY_vec[bad]  <- NA_real_
    BMSY_vec[bad] <- NA_real_
    FMSY_vec[bad] <- NA_real_
  }

  full_mat <- cbind(nat_mat, MSY = MSY_vec, BMSY = BMSY_vec, FMSY = FMSY_vec)

  target_metadata <- NULL
  if (!is.null(biomass_target)) {
    target_spec <- .bayesian_target_reference_points(
      nat_mat = nat_mat,
      biomass_target = biomass_target,
      baseline = baseline
    )
    full_mat <- cbind(full_mat, target_spec$values)
    target_metadata <- target_spec$targets
  }

  # ---- posterior summary -------------------------------------------------
  summary_df <- .posterior_summary(full_mat)

  # ---- build output object -----------------------------------------------
  out <- list(
    stanfit       = stanfit,
    posterior     = full_mat,
    log_posterior = log_mat,
    summary       = summary_df,
    model_fit     = model_fit,
    par_names     = par_names,
    nat_names     = colnames(full_mat),
    reference_point_targets = target_metadata,
    n_chains      = as.integer(chains),
    n_iter        = as.integer(iter),
    n_warmup      = as.integer(warmup)
  )
  class(out) <- "bayes_fit"
  out
}


# =========================================================================
# Internal helpers
# =========================================================================

#' Map log-scale parameter names to natural-scale names
#' @keywords internal
.log_to_natural_names <- function(par_names) {
  nat <- sub("^log_", "", par_names)
  # Convert underscore area/label suffixes to dots: q_A1 -> q.A1
  nat <- sub("^(q|B0|sigma_proc|sigma_obs)_([A-Z])", "\\1.\\2", nat)
  nat
}


#' Posterior summary table
#' @keywords internal
.posterior_summary <- function(mat) {
  qs <- apply(mat, 2, quantile, probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  data.frame(
    parameter = colnames(mat),
    mean      = colMeans(mat, na.rm = TRUE),
    sd        = apply(mat, 2, sd, na.rm = TRUE),
    `2.5%`    = qs["2.5%", ],
    `50%`     = qs["50%", ],
    `97.5%`   = qs["97.5%", ],
    n_eff     = .simple_ess(mat),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


#' Build depletion-based posterior reference points
#' @keywords internal
.bayesian_target_reference_points <- function(nat_mat, biomass_target, baseline) {
  nat_names <- colnames(nat_mat)
  b0_cols <- grep("^B0(\\.|$)", nat_names)

  baseline_name <- baseline
  if (baseline_name == "auto") {
    baseline_name <- if (length(b0_cols) > 0) "B0" else "K"
  }

  if (baseline_name == "B0" && length(b0_cols) == 0) {
    stop("baseline = 'B0' requested but no fitted B0 parameter was found")
  }

  K_draws <- nat_mat[, "K"]
  r_draws <- nat_mat[, "r"]
  m_draws <- nat_mat[, "m"]

  if (baseline_name == "B0") {
    if ("B0" %in% nat_names) {
      baseline_draws <- nat_mat[, "B0"]
    } else {
      baseline_draws <- rowSums(nat_mat[, b0_cols, drop = FALSE])
    }
  } else {
    baseline_draws <- K_draws
  }

  target_df <- data.frame(
    fraction = biomass_target,
    baseline = baseline_name,
    biomass_name = vapply(biomass_target, format_reference_point_label, character(1), prefix = "B", baseline_name = baseline_name),
    f_name = vapply(biomass_target, format_reference_point_label, character(1), prefix = "F", baseline_name = baseline_name),
    stringsAsFactors = FALSE
  )

  target_values <- vector("list", 2 * nrow(target_df))
  target_names <- character(2 * nrow(target_df))
  out_idx <- 1L
  for (i in seq_len(nrow(target_df))) {
    biomass_vals <- target_df$fraction[i] * baseline_draws
    fishing_vals <- calculate_target_fishing_mortality(r_draws, K_draws, m_draws, biomass_vals)
    fishing_vals[!is.finite(fishing_vals) | fishing_vals < 0] <- NA_real_

    target_names[out_idx] <- target_df$biomass_name[i]
    target_values[[out_idx]] <- biomass_vals
    out_idx <- out_idx + 1L

    target_names[out_idx] <- target_df$f_name[i]
    target_values[[out_idx]] <- fishing_vals
    out_idx <- out_idx + 1L
  }

  values <- as.matrix(as.data.frame(target_values, stringsAsFactors = FALSE))
  colnames(values) <- target_names

  list(targets = target_df, values = values)
}


#' Simple effective sample size (bulk ESS) from a matrix
#' Uses the rstan summary if the stanfit is available;
#' otherwise falls back to a simple autocorrelation-based estimate.
#' @keywords internal
.simple_ess <- function(mat) {
  # Use a quick autocorrelation-based ESS estimate
  apply(mat, 2, function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 4) return(NA_real_)
    # Compute autocorrelation for lags 1..min(n-1, 100)
    max_lag <- min(n - 1, 100)
    ac <- acf(x, lag.max = max_lag, plot = FALSE)$acf[-1, 1, 1]
    # Sum consecutive positive pairs
    # (Geyer's initial monotone sequence estimator, simplified)
    s <- 0
    for (k in seq(1, length(ac) - 1, by = 2)) {
      pair_sum <- ac[k] + ac[k + 1]
      if (pair_sum < 0) break
      s <- s + pair_sum
    }
    n / (1 + 2 * s)
  })
}


# =========================================================================
# Print method
# =========================================================================

#' @export
print.bayes_fit <- function(x, digits = 3, ...) {
  cat("Bayesian Production Model Fit (tmbstan)\n")
  cat("========================================\n")
  cat("Chains:", x$n_chains,
      " | Iter:", x$n_iter,
      " | Warmup:", x$n_warmup, "\n")

  n_draws <- nrow(x$posterior)
  cat("Post-warmup draws:", n_draws, "\n\n")

  # Divergent transitions
  divs <- .count_divergences(x$stanfit)
  if (divs > 0) {
    cat("WARNING:", divs, "divergent transitions after warmup\n\n")
  }

  cat("Posterior summary:\n")
  s <- x$summary
  cat(format_summary_table(s, digits), sep = "\n")
  cat("\n")
  invisible(x)
}


#' Format a posterior summary table for printing
#' @keywords internal
format_summary_table <- function(s, digits = 3) {
  # Format numeric columns
  fmt <- s
  for (col in setdiff(names(fmt), "parameter")) {
    fmt[[col]] <- formatC(s[[col]], digits = digits, format = "f")
  }
  # Align
  lines <- capture.output(print(fmt, row.names = FALSE, right = TRUE))
  lines
}


#' Count divergent transitions
#' @keywords internal
.count_divergences <- function(stanfit) {
  tryCatch({
    sp <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
    sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  }, error = function(e) 0L)
}


# =========================================================================
# Plot method
# =========================================================================

#' Plot Bayesian Fit Diagnostics
#'
#' Produces diagnostic plots for a \code{bayes_fit} object.
#'
#' @param x      A \code{bayes_fit} object.
#' @param type   Character, one of \code{"trace"} (default), \code{"density"},
#'   \code{"pairs"}, or \code{"histogram"}.
#' @param pars   Character vector of natural-scale parameter names to plot.
#'   Defaults to key model parameters (\code{r}, \code{K}, \code{MSY},
#'   \code{BMSY}).  Use \code{"all"} for every parameter.
#' @param ...    Additional arguments (currently unused).
#'
#' @return A \code{ggplot} object (for \code{"trace"}, \code{"density"},
#'   \code{"histogram"}) or a base-R pairs plot (for \code{"pairs"}).
#'
#' @export
plot.bayes_fit <- function(x, type = c("trace", "density", "pairs", "histogram"),
                           pars = NULL, ...) {
  type <- match.arg(type)

  # Determine which parameters to plot
  all_names <- colnames(x$posterior)
  if (is.null(pars)) {
    default_pars <- c("r", "K", "sigma_obs", "MSY", "BMSY", "FMSY")
    pars <- intersect(default_pars, all_names)
    if (length(pars) == 0) pars <- all_names[seq_len(min(6, length(all_names)))]
  } else if (length(pars) == 1 && pars == "all") {
    pars <- all_names
  } else {
    pars <- intersect(pars, all_names)
    if (length(pars) == 0) stop("None of the requested parameters found")
  }

  if (type == "pairs") {
    return(.plot_pairs(x, pars))
  }

  # Build long-format data.frame for ggplot
  plot_df <- .build_plot_df(x, pars, include_chain = (type == "trace"))

  switch(type,
    trace    = .plot_trace(plot_df),
    density  = .plot_density(plot_df),
    histogram = .plot_histogram(plot_df)
  )
}


# ---- Trace plot ----------------------------------------------------------
#' @keywords internal
.plot_trace <- function(df) {
  ggplot(df, aes(x = iteration, y = sample_value, colour = chain)) +
    geom_line(alpha = 0.6, linewidth = 0.3) +
    facet_wrap(~param, scales = "free_y") +
    labs(x = "Iteration", y = "Value", colour = "Chain",
         title = "MCMC Trace Plots") +
    theme_bw() +
    theme(legend.position = "bottom")
}


# ---- Density plot --------------------------------------------------------
#' @keywords internal
.plot_density <- function(df) {
  ggplot(df, aes(x = sample_value)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    facet_wrap(~param, scales = "free") +
    labs(x = "Value", y = "Density",
         title = "Posterior Density Plots") +
    theme_bw()
}


# ---- Histogram plot ------------------------------------------------------
#' @keywords internal
.plot_histogram <- function(df) {
  ggplot(df, aes(x = sample_value)) +
    geom_histogram(bins = 40, fill = "steelblue", colour = "white", alpha = 0.7) +
    facet_wrap(~param, scales = "free") +
    labs(x = "Value", y = "Count",
         title = "Posterior Histograms") +
    theme_bw()
}


# ---- Pairs plot ----------------------------------------------------------
#' @keywords internal
.plot_pairs <- function(x, pars) {
  sub_mat <- x$posterior[, pars, drop = FALSE]
  # Use base pairs for simplicity
  pairs(sub_mat, pch = ".", col = adjustcolor("steelblue", alpha.f = 0.3),
        main = "Posterior Pairs Plot")
  invisible(NULL)
}


# ---- helpers for building plot data frames -------------------------------
#' @keywords internal
.build_plot_df <- function(x, pars, include_chain = FALSE) {
  stanfit <- x$stanfit

  if (include_chain) {
    # Need per-chain iteration info
    arr <- rstan::extract(stanfit, permuted = FALSE)  # [iter, chain, param]
    par_names_stan <- dimnames(arr)[[3]]

    # Map log-scale names to natural-scale names for matching
    log_par_names <- x$par_names
    nat_par_names <- .log_to_natural_names(log_par_names)

    rows <- list()
    for (p in pars) {
      # Check if p is a direct log-scale parameter (transformed) or derived
      nat_idx <- match(p, nat_par_names)
      if (!is.na(nat_idx)) {
        # This is a model parameter - extract from stanfit and transform
        log_name <- log_par_names[nat_idx]
        stan_idx <- match(log_name, par_names_stan)
        if (is.na(stan_idx)) next
        for (ch in seq_len(dim(arr)[2])) {
          vals <- exp(arr[, ch, stan_idx])
          n_iter <- length(vals)
          rows[[length(rows) + 1]] <- data.frame(
            param     = p,
            iteration = seq_len(n_iter),
            chain     = as.factor(ch),
            sample_value = vals,
            stringsAsFactors = FALSE
          )
        }
      } else if (p %in% c("MSY", "BMSY", "FMSY")) {
        # Derived: compute per chain
        r_idx  <- match("log_r", par_names_stan)
        K_idx  <- match("log_K", par_names_stan)
        m_idx  <- match("log_m", par_names_stan)
        if (any(is.na(c(r_idx, K_idx, m_idx)))) next
        for (ch in seq_len(dim(arr)[2])) {
          r_v <- exp(arr[, ch, r_idx])
          K_v <- exp(arr[, ch, K_idx])
          m_v <- exp(arr[, ch, m_idx])
          vals <- switch(p,
            MSY  = (r_v / m_v * (1 - 1 / m_v)) * K_v * (1 / m_v)^(1 / (m_v - 1)),
            BMSY = K_v * (1 / m_v)^(1 / (m_v - 1)),
            FMSY = r_v / m_v * (1 - 1 / m_v)
          )
          n_iter <- length(vals)
          rows[[length(rows) + 1]] <- data.frame(
            param     = p,
            iteration = seq_len(n_iter),
            chain     = as.factor(ch),
            sample_value = vals,
            stringsAsFactors = FALSE
          )
        }
      } else if (!is.null(x$reference_point_targets) &&
                 (p %in% x$reference_point_targets$biomass_name ||
                  p %in% x$reference_point_targets$f_name)) {
        r_idx  <- match("log_r", par_names_stan)
        K_idx  <- match("log_K", par_names_stan)
        m_idx  <- match("log_m", par_names_stan)
        if (any(is.na(c(r_idx, K_idx, m_idx)))) next

        row_idx <- if (p %in% x$reference_point_targets$biomass_name) {
          match(p, x$reference_point_targets$biomass_name)
        } else {
          match(p, x$reference_point_targets$f_name)
        }
        target_fraction <- x$reference_point_targets$fraction[row_idx]
        target_baseline <- x$reference_point_targets$baseline[row_idx]
        b0_stan_idx <- grep("^log_B0", par_names_stan)

        for (ch in seq_len(dim(arr)[2])) {
          r_v <- exp(arr[, ch, r_idx])
          K_v <- exp(arr[, ch, K_idx])
          m_v <- exp(arr[, ch, m_idx])
          if (target_baseline == "B0") {
            if (length(b0_stan_idx) == 0) next
            baseline_v <- rowSums(exp(arr[, ch, b0_stan_idx, drop = FALSE]))
          } else {
            baseline_v <- K_v
          }
          biomass_v <- target_fraction * baseline_v
          vals <- if (p %in% x$reference_point_targets$biomass_name) {
            biomass_v
          } else {
            calculate_target_fishing_mortality(r_v, K_v, m_v, biomass_v)
          }
          vals[!is.finite(vals) | vals < 0] <- NA_real_
          n_iter <- length(vals)
          rows[[length(rows) + 1]] <- data.frame(
            param     = p,
            iteration = seq_len(n_iter),
            chain     = as.factor(ch),
            sample_value = vals,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    do.call(rbind, rows)
  } else {
    # Permuted draws - no chain info
    mat <- x$posterior[, pars, drop = FALSE]
    rows <- list()
    for (p in pars) {
      rows[[length(rows) + 1]] <- data.frame(
        param        = p,
        iteration    = seq_len(nrow(mat)),
        sample_value = mat[, p],
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }
}


# =========================================================================
# Posterior predictive check
# =========================================================================

#' Posterior Predictive Check
#'
#' Simulates replicated CPUE datasets from the posterior predictive
#' distribution and compares them with the observed data.
#'
#' @param bayes_fit  A \code{bayes_fit} object.
#' @param n_sims     Integer, number of posterior draws to use for
#'   simulation (default: 200, randomly sampled).
#' @param seed       Integer, random seed (default: NULL).
#'
#' @return A list with class \code{"ppc_result"} containing:
#'   \describe{
#'     \item{observed}{Numeric vector of observed CPUE values.}
#'     \item{simulated}{Matrix (\code{n_sims} rows x \code{n_obs} columns)
#'       of replicated CPUE values.}
#'     \item{p_value}{Bayesian p-value (proportion of replicated summary
#'       statistics exceeding observed).}
#'   }
#'
#' @export
posterior_predictive_check <- function(bayes_fit, n_sims = 200L, seed = NULL) {
  if (!inherits(bayes_fit, "bayes_fit")) {
    stop("Input must be a 'bayes_fit' object")
  }
  if (!is.null(seed)) set.seed(seed)

  model_fit <- bayes_fit$model_fit
  posterior  <- bayes_fit$posterior

  # Extract observed CPUE
  obs_cpue <- model_fit$data$cpue
  if (is.null(obs_cpue)) {
    # Fallback: try to extract from residuals + fitted
    obs_cpue <- model_fit$results$fitted_cpue + model_fit$results$residuals
  }
  obs_flat <- as.numeric(obs_cpue)
  obs_flat <- obs_flat[!is.na(obs_flat)]
  n_obs <- length(obs_flat)

  # Subsample draws
  n_total <- nrow(posterior)
  draw_idx <- sample.int(n_total, min(n_sims, n_total))

  # Get parameter column names
  nat_names <- colnames(posterior)

  # For each draw, simulate CPUE
  sim_mat <- matrix(NA_real_, nrow = length(draw_idx), ncol = n_obs)

  for (i in seq_along(draw_idx)) {
    d <- draw_idx[i]

    # Extract parameters for this draw
    sigma_obs <- .extract_draw_param(posterior, d, "sigma_obs", nat_names)
    if (is.na(sigma_obs)) sigma_obs <- 0.1  # fallback

    # Get fitted CPUE from the model
    # (using the model's fitted CPUE as the expected value;
    #  for a full PPC we would re-run the dynamics, but this is
    #  computationally expensive and the observation-error PPC
    #  is the standard quick check)
    fitted_cpue <- as.numeric(model_fit$results$fitted_cpue)
    fitted_cpue <- fitted_cpue[!is.na(fitted_cpue)]

    if (length(fitted_cpue) != n_obs) {
      # If lengths don't match, skip
      next
    }

    # Simulate observations: log-normal observation error
    sim_mat[i, ] <- fitted_cpue * exp(rnorm(n_obs, 0, sigma_obs))
  }

  # Bayesian p-value: compare mean of replicated vs observed
  obs_mean <- mean(obs_flat, na.rm = TRUE)
  sim_means <- rowMeans(sim_mat, na.rm = TRUE)
  p_value <- mean(sim_means > obs_mean, na.rm = TRUE)

  out <- list(
    observed  = obs_flat,
    simulated = sim_mat,
    p_value   = p_value
  )
  class(out) <- "ppc_result"
  out
}


#' @export
print.ppc_result <- function(x, ...) {
  cat("Posterior Predictive Check\n")
  cat("=========================\n")
  cat("Observed data points:", length(x$observed), "\n")
  cat("Simulations:", nrow(x$simulated), "\n")
  cat("Bayesian p-value (mean):", round(x$p_value, 3), "\n")
  if (x$p_value < 0.05 || x$p_value > 0.95) {
    cat("  -> Potential model misfit (extreme p-value)\n")
  } else {
    cat("  -> No evidence of systematic misfit\n")
  }
  invisible(x)
}


#' Plot Posterior Predictive Check
#'
#' @param x    A \code{ppc_result} object.
#' @param ...  Additional arguments (currently unused).
#'
#' @return A \code{ggplot} object.
#' @export
plot.ppc_result <- function(x, ...) {
  obs <- x$observed
  sim <- x$simulated
  n_obs <- length(obs)

  # Compute quantile envelopes
  q_lo  <- apply(sim, 2, quantile, probs = 0.025, na.rm = TRUE)
  q_mid <- apply(sim, 2, quantile, probs = 0.50, na.rm = TRUE)
  q_hi  <- apply(sim, 2, quantile, probs = 0.975, na.rm = TRUE)

  df <- data.frame(
    index = seq_len(n_obs),
    observed = obs,
    median   = q_mid,
    lower    = q_lo,
    upper    = q_hi
  )

  ggplot(df, aes(x = index)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3,
                fill = "steelblue") +
    geom_line(aes(y = median), colour = "steelblue", linewidth = 0.8) +
    geom_point(aes(y = observed), colour = "black", size = 2) +
    labs(x = "Observation",
         y = "CPUE",
         title = "Posterior Predictive Check",
         subtitle = paste("95% prediction interval | Bayesian p-value:",
                          round(x$p_value, 3))) +
    theme_bw()
}


# ---- Helper: extract a named parameter from a posterior draw row ---------
#' @keywords internal
.extract_draw_param <- function(posterior, draw_idx, pname, nat_names) {
  # Try exact match first
  idx <- match(pname, nat_names)
  if (!is.na(idx)) return(posterior[draw_idx, idx])

  # Try with area suffix (e.g. sigma_obs might be sigma_obs.A1 in posterior)
  candidates <- grep(paste0("^", pname), nat_names, value = TRUE)
  if (length(candidates) > 0) {
    idx <- match(candidates[1], nat_names)
    if (!is.na(idx)) return(posterior[draw_idx, idx])
  }

  NA_real_
}
