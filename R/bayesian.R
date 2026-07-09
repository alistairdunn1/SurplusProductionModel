#' Bayesian Inference via SparseNUTS
#'
#' Functions for running MCMC sampling on a fitted
#' \code{\link{ProductionModel}} using the \pkg{SparseNUTS} package,
#' which couples a TMB/RTMB objective function with Stan's
#' sparse No-U-Turn Sampler (SNUTS).
#'
#' @name bayesian
NULL

# ---- global variables used in ggplot aes / data.table --------------------
utils::globalVariables(c(
  "chain", "iteration", "param", "sample_value",
  "Var1", "Var2", "Freq", "index",
  "observed", "median", "lower", "upper"
))


#' Bayesian MCMC Sampling for a Production Model
#'
#' Draws posterior samples from the joint parameter distribution of a
#' fitted \code{ProductionModel} using Stan's sparse NUTS sampler via
#' \code{\link[SparseNUTS]{sample_snuts}}.  The starting point is the maximum
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
#'   `"B_40%K"` or `"F_40%K"`.
#' @param baseline Character string indicating which biomass baseline to use
#'   for `biomass_target`: `"auto"` (default), `"B_initial"`, or `"K"`.
#' @param globals  Named list of objects to pass to parallel R sessions
#'   (required when RTMB model closures reference external data or functions).
#'   Passed to \code{\link[SparseNUTS]{sample_snuts}} as-is.  Default
#'   \code{NULL} is suitable for self-contained models.
#' @param ...     Additional arguments passed to
#'   \code{\link[SparseNUTS]{sample_snuts}} (e.g. \code{metric},
#'   \code{model_name}).
#'
#' @return An S3 object of class \code{"bayes_fit"} containing:
#'   \describe{
#'     \item{snuts_fit}{The \code{tmbfit} object from \pkg{SparseNUTS}.}
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
#' Divergences and sampler parameters are available via
#' \code{SparseNUTS::extract_sampler_params(x$snuts_fit)}.
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
                         iter = 2000L,
                         warmup = floor(iter / 2),
                         thin = 1L,
                         seed = NULL,
                         lower = numeric(0),
                         upper = numeric(0),
                         laplace = FALSE,
                         biomass_target = NULL,
                         baseline = c("auto", "B_initial", "K"),
                         globals = NULL,
                         ...) {
  defaults <- resolve_reference_point_defaults(
    biomass_target = biomass_target,
    baseline = baseline,
    biomass_target_missing = missing(biomass_target),
    baseline_missing = missing(baseline)
  )
  biomass_target <- defaults$biomass_target
  baseline <- defaults$baseline

  # ---- input validation -------------------------------------------------
  if (!inherits(model_fit, "ProductionModel")) {
    stop("Input must be a fitted ProductionModel object")
  }
  if (!model_fit$fitted) {
    stop("Model must be fitted before running Bayesian sampling")
  }
  if (!requireNamespace("SparseNUTS", quietly = TRUE)) {
    stop(
      "Package 'SparseNUTS' is required for Bayesian fitting. ",
      "Install it with: remotes::install_github('noaa-afsc/SparseNUTS')"
    )
  }

  obj <- model_fit$results$rtmb_obj
  if (is.null(obj)) {
    stop("rtmb_obj not found in results -- was the model fitted with RTMB?")
  }

  opt_par <- model_fit$results$opt_par
  if (is.null(opt_par)) {
    stop(
      "Optimised parameters (opt_par) not found. ",
      "Re-fit the model with the latest package version."
    )
  }

  # ---- run SparseNUTS ----------------------------------------------------
  num_samples <- as.integer(iter) - as.integer(warmup)
  snuts_args <- list(
    obj          = obj,
    num_samples  = num_samples,
    num_warmup   = as.integer(warmup),
    chains       = as.integer(chains),
    thin         = as.integer(thin),
    laplace      = laplace,
    init         = "last.par.best"
  )
  if (!is.null(seed)) snuts_args$seed <- as.integer(seed)
  if (!is.null(globals)) snuts_args$globals <- globals

  # Merge user ... into snuts_args
  dots <- list(...)
  snuts_args <- modifyList(snuts_args, dots)

  snuts_fit <- do.call(SparseNUTS::sample_snuts, snuts_args)

  # ---- extract posterior samples -----------------------------------------
  par_names <- names(opt_par)

  # Collapse chains into a matrix (all post-warmup draws, log-scale parameters)
  log_mat <- as.matrix(as.data.frame(snuts_fit))
  # Drop lp__ column if present
  lp_col <- which(colnames(log_mat) == "lp__")
  if (length(lp_col) > 0) log_mat <- log_mat[, -lp_col, drop = FALSE]

  # When random effects are sampled (full Bayesian, laplace = FALSE) the SNUTS
  # posterior includes the random-effect columns (e.g. proc_dev) alongside the
  # fixed-effect scalar parameters. Retain only the fixed-effect columns that
  # correspond to opt_par; the random effects are nuisance parameters and are
  # excluded from the natural-scale transform and posterior summary.
  fixed_cols <- which(colnames(log_mat) %in% par_names)
  if (length(fixed_cols) > 0 && length(fixed_cols) < ncol(log_mat)) {
    log_mat <- log_mat[, fixed_cols, drop = FALSE]
  }

  # ---- transform to natural scale + derived quantities -------------------
  # Derive names from the retained columns; their order need not match opt_par.
  par_names <- colnames(log_mat)
  nat_names <- .log_to_natural_names(par_names)
  # Transform each column from its estimation scale to the natural scale. Most
  # parameters are on the log scale, but movement_rate is on the logit scale
  # (log_movement_rate = qlogis(movement_rate)) and the AR(1) coefficient is on
  # the atanh scale (theta_rho = atanh(rho)); these need the matching inverse
  # transforms rather than exp().
  nat_mat <- log_mat
  for (j in seq_along(par_names)) {
    nat_mat[, j] <- .natural_scale_value(par_names[j], log_mat[, j])
  }
  colnames(nat_mat) <- nat_names

  # Derived quantities: MSY, BMSY, FMSY
  # Need: r, K, m
  r_col <- which(nat_names == "r")
  K_col <- which(nat_names == "K")
  m_col <- which(nat_names == "m")

  n_draws <- nrow(nat_mat)
  MSY_vec <- rep(NA_real_, n_draws)
  BMSY_vec <- rep(NA_real_, n_draws)
  FMSY_vec <- rep(NA_real_, n_draws)

  if (length(r_col) == 1 && length(K_col) == 1 && length(m_col) == 1) {
    r_draws <- nat_mat[, r_col]
    K_draws <- nat_mat[, K_col]
    m_draws <- nat_mat[, m_col]

    # Standard Pella-Tomlinson: BMSY = K m^(-1/(m-1)), FMSY = r/m.
    BMSY_vec <- K_draws * (1 / m_draws)^(1 / (m_draws - 1))
    FMSY_vec <- r_draws / m_draws
    MSY_vec <- FMSY_vec * BMSY_vec

    # Handle near-singular cases (m near 1) which can produce Inf/NaN
    # and biologically invalid cases (m < 1 gives negative FMSY)
    bad <- !is.finite(MSY_vec) | !is.finite(BMSY_vec) | !is.finite(FMSY_vec) |
      MSY_vec <= 0 | BMSY_vec <= 0 | FMSY_vec <= 0
    MSY_vec[bad] <- NA_real_
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

  # ---- chain identifiers for convergence diagnostics ---------------------
  # Draws in full_mat are row-aligned with the SNUTS sampler output; recover
  # the per-draw chain index so split-Rhat and bulk-ESS can be computed.
  chain_id <- tryCatch(
    {
      sp <- SparseNUTS::extract_sampler_params(snuts_fit)
      if (!is.null(sp) && "chain" %in% names(sp) && nrow(sp) == nrow(full_mat)) {
        as.integer(sp$chain)
      } else {
        NULL
      }
    },
    error = function(e) NULL
  )
  if (is.null(chain_id)) {
    # Fall back to an even split across the requested number of chains.
    n_draws_total <- nrow(full_mat)
    n_ch <- max(1L, as.integer(chains))
    per <- ceiling(n_draws_total / n_ch)
    chain_id <- rep(seq_len(n_ch), each = per)[seq_len(n_draws_total)]
  }

  # ---- posterior summary -------------------------------------------------
  summary_df <- .posterior_summary(full_mat, chain_id)

  # ---- build output object -----------------------------------------------
  out <- list(
    snuts_fit = snuts_fit,
    posterior = full_mat,
    log_posterior = log_mat,
    summary = summary_df,
    model_fit = model_fit,
    par_names = par_names,
    nat_names = colnames(full_mat),
    reference_point_targets = target_metadata,
    n_chains = as.integer(chains),
    n_iter = as.integer(iter),
    n_warmup = as.integer(warmup)
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
  # Estimation-scale parameters that are not simple log transforms.
  nat <- ifelse(par_names == "log_movement_rate", "movement_rate",
    ifelse(par_names == "theta_rho", "rho", sub("^log_", "", par_names))
  )
  # Convert underscore area/label suffixes to dots: q_A1 -> q.A1, K_A1 -> K.A1
  nat <- sub("^(q|K|B_initial|sigma_proc|sigma_obs)_([A-Z])", "\\1.\\2", nat)
  nat
}


#' Transform a draw (or vector of draws) from the estimation scale to the
#' natural scale, honouring each parameter's link function.
#'
#' Most parameters use a log link, but \code{log_movement_rate} uses a logit
#' link (\code{movement_rate = plogis(theta)}) and \code{theta_rho} uses an
#' atanh link (\code{rho = tanh(theta)}). Mirrors
#' \code{transform_parameters_to_natural} for the MCMC draws.
#'
#' @param name Estimation-scale parameter name.
#' @param x Numeric draw(s) on the estimation scale.
#' @keywords internal
.natural_scale_value <- function(name, x) {
  if (identical(name, "log_movement_rate")) {
    stats::plogis(x)
  } else if (identical(name, "theta_rho")) {
    tanh(x)
  } else if (grepl("^log_", name)) {
    exp(x)
  } else {
    x
  }
}


#' Posterior summary table
#'
#' Includes rank-normalised split-Rhat and bulk effective sample size
#' (Vehtari et al. 2021) when per-draw chain identifiers are supplied.
#'
#' @param mat Matrix of posterior draws (rows = draws, columns = quantities).
#' @param chain_id Optional integer vector of chain identifiers, one per row of
#'   \code{mat}. When \code{NULL} or a single chain, between-chain diagnostics
#'   are returned as \code{NA}.
#' @keywords internal
.posterior_summary <- function(mat, chain_id = NULL) {
  qs <- apply(mat, 2, quantile, probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  diag <- .mcmc_diagnostics(mat, chain_id)
  data.frame(
    parameter = colnames(mat),
    mean = colMeans(mat, na.rm = TRUE),
    sd = apply(mat, 2, sd, na.rm = TRUE),
    `2.5%` = qs["2.5%", ],
    `50%` = qs["50%", ],
    `97.5%` = qs["97.5%", ],
    n_eff = diag$ess_bulk,
    Rhat = diag$rhat,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


#' MCMC convergence diagnostics (rank-normalised split-Rhat and bulk ESS)
#'
#' Computes, per column of \code{mat}, the rank-normalised split-Rhat and the
#' bulk effective sample size following Vehtari et al. (2021). When chain
#' structure is unavailable (single chain, or unequal/short chains) Rhat is
#' returned as \code{NA} and a pooled autocorrelation ESS is used as a
#' fallback.
#'
#' @param mat Matrix of draws (rows = draws, cols = quantities).
#' @param chain_id Integer vector of chain identifiers, one per row.
#' @return List with numeric vectors \code{rhat} and \code{ess_bulk}.
#' @keywords internal
.mcmc_diagnostics <- function(mat, chain_id = NULL) {
  n_col <- ncol(mat)
  rhat <- rep(NA_real_, n_col)
  ess_bulk <- rep(NA_real_, n_col)

  usable <- !is.null(chain_id) && length(chain_id) == nrow(mat) &&
    length(unique(chain_id)) >= 2L
  if (usable) {
    split_lengths <- tapply(seq_along(chain_id), chain_id, length)
    n_per <- min(split_lengths)
    usable <- is.finite(n_per) && n_per >= 4L
  }

  if (!usable) {
    # Fallback: pooled autocorrelation-based ESS, no between-chain Rhat.
    ess_bulk <- .simple_ess(mat)
    return(list(rhat = rhat, ess_bulk = ess_bulk))
  }

  chains <- sort(unique(chain_id))
  for (j in seq_len(n_col)) {
    # Assemble an [iterations x chains] matrix, truncated to the common length.
    draws <- vapply(chains, function(cc) {
      x <- mat[chain_id == cc, j]
      x[seq_len(n_per)]
    }, numeric(n_per))
    if (!is.matrix(draws)) draws <- matrix(draws, nrow = n_per)

    if (any(!is.finite(draws))) {
      next
    }
    # Bulk diagnostics operate on rank-normalised draws.
    ranks <- .rank_normalise(as.numeric(draws))
    z <- matrix(ranks, nrow = n_per)

    rhat[j] <- .split_rhat(z)
    ess_bulk[j] <- .ess_multichain(z)
  }

  list(rhat = rhat, ess_bulk = ess_bulk)
}


#' Rank-normalise a numeric vector (Blom transform)
#' @keywords internal
.rank_normalise <- function(x) {
  n <- length(x)
  r <- rank(x, ties.method = "average")
  stats::qnorm((r - 0.375) / (n + 0.25))
}


#' Split-Rhat from an [iterations x chains] matrix
#' @keywords internal
.split_rhat <- function(z) {
  n <- nrow(z)
  m <- ncol(z)
  half <- floor(n / 2)
  if (half < 2L) {
    return(NA_real_)
  }
  # Split each chain in half -> 2m sub-chains of length `half`.
  splits <- cbind(z[seq_len(half), , drop = FALSE], z[(n - half + 1L):n, , drop = FALSE])
  chain_means <- colMeans(splits)
  chain_vars <- apply(splits, 2, stats::var)
  w <- mean(chain_vars)
  b <- half * stats::var(chain_means)
  if (!is.finite(w) || w <= 0) {
    return(NA_real_)
  }
  var_plus <- ((half - 1) * w + b) / half
  sqrt(var_plus / w)
}


#' Bulk effective sample size from an [iterations x chains] matrix
#'
#' Uses the multi-chain autocovariance estimator with Geyer's initial
#' monotone positive sequence truncation (Vehtari et al. 2021).
#' @keywords internal
.ess_multichain <- function(z) {
  n <- nrow(z)
  m <- ncol(z)
  if (n < 4L) {
    return(NA_real_)
  }

  chain_means <- colMeans(z)
  chain_vars <- apply(z, 2, stats::var)
  w <- mean(chain_vars)
  if (!is.finite(w) || w <= 0) {
    return(NA_real_)
  }
  b <- n * stats::var(chain_means)
  var_plus <- ((n - 1) * w + b) / n

  # Mean (over chains) autocorrelation at each lag via each chain's acf.
  max_lag <- n - 1L
  acov <- matrix(0, nrow = max_lag + 1L, ncol = m)
  for (cc in seq_len(m)) {
    a <- stats::acf(z[, cc], lag.max = max_lag, plot = FALSE, demean = TRUE)$acf[, 1, 1]
    # Convert autocorrelation to autocovariance using the chain variance.
    acov[seq_along(a), cc] <- a * chain_vars[cc]
  }
  mean_acov <- rowMeans(acov)
  rho <- 1 - (w - mean_acov[-1]) / var_plus # rho_t for t = 1..max_lag

  # Geyer initial positive/monotone sequence on paired lags.
  rho_sum <- 0
  t <- 1L
  prev_pair <- Inf
  while (t + 1L <= length(rho)) {
    pair <- rho[t] + rho[t + 1L]
    if (pair < 0) break
    pair <- min(pair, prev_pair) # enforce monotone non-increasing
    rho_sum <- rho_sum + pair
    prev_pair <- pair
    t <- t + 2L
  }

  tau <- 1 + 2 * rho_sum
  if (!is.finite(tau) || tau <= 0) {
    return(NA_real_)
  }
  (n * m) / tau
}


#' Build depletion-based posterior reference points
#' @keywords internal
.bayesian_target_reference_points <- function(nat_mat, biomass_target, baseline) {
  nat_names <- colnames(nat_mat)
  has_d0 <- "d0" %in% nat_names

  baseline_name <- baseline
  if (baseline_name == "auto") {
    baseline_name <- if (has_d0) "B_initial" else "K"
  }

  if (baseline_name == "B_initial" && !has_d0) {
    stop("baseline = 'B_initial' requested but no fitted d0 parameter was found")
  }

  K_draws <- nat_mat[, "K"]
  r_draws <- nat_mat[, "r"]
  m_draws <- nat_mat[, "m"]

  # B_initial = d0 * K (initial-depletion parameter times carrying capacity).
  if (baseline_name == "B_initial") {
    baseline_draws <- nat_mat[, "d0"] * K_draws
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
#' Uses an autocorrelation-based estimate.
#' @keywords internal
.simple_ess <- function(mat) {
  # Use a quick autocorrelation-based ESS estimate
  apply(mat, 2, function(x) {
    x <- x[is.finite(x)]
    n <- length(x)
    if (n < 4) {
      return(NA_real_)
    }
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
  cat("Bayesian Production Model Fit (SparseNUTS)\n")
  cat("============================================\n")
  cat(
    "Chains:", x$n_chains,
    " | Iter:", x$n_iter,
    " | Warmup:", x$n_warmup, "\n"
  )

  n_draws <- nrow(x$posterior)
  cat("Post-warmup draws:", n_draws, "\n\n")

  # Divergent transitions
  divs <- .count_divergences(x$snuts_fit)
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
.count_divergences <- function(snuts_fit) {
  tryCatch(
    {
      sp <- SparseNUTS::extract_sampler_params(snuts_fit)
      sum(sp[["divergent__"]], na.rm = TRUE)
    },
    error = function(e) 0L
  )
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
    trace = .plot_trace(plot_df),
    density = .plot_density(plot_df),
    histogram = .plot_histogram(plot_df)
  )
}


# ---- Trace plot ----------------------------------------------------------
#' @keywords internal
.plot_trace <- function(df) {
  ggplot(df, aes(x = iteration, y = sample_value, colour = chain)) +
    geom_line(alpha = 0.6, linewidth = 0.3) +
    facet_wrap(~param, scales = "free_y") +
    labs(
      x = "Iteration", y = "Value", colour = "Chain",
      title = "MCMC Trace Plots"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}


# ---- Density plot --------------------------------------------------------
#' @keywords internal
.plot_density <- function(df) {
  ggplot(df, aes(x = sample_value)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    facet_wrap(~param, scales = "free") +
    labs(
      x = "Value", y = "Density",
      title = "Posterior Density Plots"
    ) +
    theme_bw()
}


# ---- Histogram plot ------------------------------------------------------
#' @keywords internal
.plot_histogram <- function(df) {
  ggplot(df, aes(x = sample_value)) +
    geom_histogram(bins = 40, fill = "steelblue", colour = "white", alpha = 0.7) +
    facet_wrap(~param, scales = "free") +
    labs(
      x = "Value", y = "Count",
      title = "Posterior Histograms"
    ) +
    theme_bw()
}


# ---- Pairs plot ----------------------------------------------------------
#' @keywords internal
.plot_pairs <- function(x, pars) {
  sub_mat <- x$posterior[, pars, drop = FALSE]
  # Use base pairs for simplicity
  pairs(sub_mat,
    pch = ".", col = adjustcolor("steelblue", alpha.f = 0.3),
    main = "Posterior Pairs Plot"
  )
  invisible(NULL)
}


# ---- helpers for building plot data frames -------------------------------
#' @keywords internal
.build_plot_df <- function(x, pars, include_chain = FALSE) {
  if (include_chain) {
    # Retrieve per-chain iteration info from SparseNUTS
    sp <- tryCatch(
      SparseNUTS::extract_sampler_params(x$snuts_fit),
      error = function(e) NULL
    )
    # Retrieve log-scale samples (rows match sp rows, same order)
    log_samp <- as.matrix(as.data.frame(x$snuts_fit))
    lp_col <- which(colnames(log_samp) == "lp__")
    if (length(lp_col) > 0) log_samp <- log_samp[, -lp_col, drop = FALSE]

    # Determine chain and within-chain iteration for each draw
    if (!is.null(sp) && "chain" %in% names(sp) && nrow(sp) == nrow(log_samp)) {
      chain_id <- as.factor(sp$chain)
      iter_id <- sp$iteration
    } else {
      # Fallback: assign draws sequentially across chains
      n_draws <- nrow(log_samp)
      n_chains <- x$n_chains
      n_per <- ceiling(n_draws / n_chains)
      chain_id <- as.factor(rep(seq_len(n_chains), each = n_per)[seq_len(n_draws)])
      iter_id <- rep(seq_len(n_per), times = n_chains)[seq_len(n_draws)]
    }

    log_par_names <- x$par_names
    nat_par_names <- .log_to_natural_names(log_par_names)

    rows <- list()
    for (p in pars) {
      nat_idx <- match(p, nat_par_names)
      if (!is.na(nat_idx)) {
        # Direct model parameter — transform the estimation-scale samples using
        # the parameter's link (log, logit for movement_rate, atanh for rho).
        log_name <- log_par_names[nat_idx]
        samp_idx <- match(log_name, colnames(log_samp))
        if (is.na(samp_idx)) next
        vals <- .natural_scale_value(log_name, log_samp[, samp_idx])
      } else if (p %in% c("MSY", "BMSY", "FMSY")) {
        # Derived reference points — compute from sampled parameters
        r_col <- match(grep("^log_r(\\.|$)", colnames(log_samp), value = TRUE)[1], colnames(log_samp))
        K_col <- match(grep("^log_K(\\.|$)", colnames(log_samp), value = TRUE)[1], colnames(log_samp))
        m_col <- match(grep("^log_m(\\.|$)", colnames(log_samp), value = TRUE)[1], colnames(log_samp))
        if (any(is.na(c(r_col, K_col, m_col)))) next
        r_v <- exp(log_samp[, r_col])
        K_v <- exp(log_samp[, K_col])
        m_v <- exp(log_samp[, m_col])
        vals <- switch(p,
          MSY  = (r_v / m_v) * K_v * (1 / m_v)^(1 / (m_v - 1)),
          BMSY = K_v * (1 / m_v)^(1 / (m_v - 1)),
          FMSY = r_v / m_v
        )
      } else if (!is.null(x$reference_point_targets) &&
        (p %in% x$reference_point_targets$biomass_name ||
          p %in% x$reference_point_targets$f_name)) {
        r_col <- match(grep("^log_r(\\.|$)", colnames(log_samp), value = TRUE)[1], colnames(log_samp))
        K_col <- match(grep("^log_K(\\.|$)", colnames(log_samp), value = TRUE)[1], colnames(log_samp))
        m_col <- match(grep("^log_m(\\.|$)", colnames(log_samp), value = TRUE)[1], colnames(log_samp))
        if (any(is.na(c(r_col, K_col, m_col)))) next

        row_idx <- if (p %in% x$reference_point_targets$biomass_name) {
          match(p, x$reference_point_targets$biomass_name)
        } else {
          match(p, x$reference_point_targets$f_name)
        }
        target_fraction <- x$reference_point_targets$fraction[row_idx]
        target_baseline <- x$reference_point_targets$baseline[row_idx]

        r_v <- exp(log_samp[, r_col])
        K_v <- exp(log_samp[, K_col])
        m_v <- exp(log_samp[, m_col])
        if (target_baseline == "B_initial") {
          d0_col <- grep("^log_d0$", colnames(log_samp))
          if (length(d0_col) == 0) next
          baseline_v <- exp(log_samp[, d0_col[1]]) * K_v
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
      } else {
        next
      }

      rows[[length(rows) + 1]] <- data.frame(
        param = p,
        iteration = iter_id,
        chain = chain_id,
        sample_value = vals,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  } else {
    # Permuted draws — no chain info needed
    mat <- x$posterior[, pars, drop = FALSE]
    rows <- list()
    for (p in pars) {
      rows[[length(rows) + 1]] <- data.frame(
        param = p,
        iteration = seq_len(nrow(mat)),
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
  posterior <- bayes_fit$posterior

  # Extract observed CPUE
  obs_cpue <- model_fit$data$cpue
  if (is.null(obs_cpue)) {
    stop("Observed CPUE must be present in model_fit$data$cpue for posterior predictive checks")
  }
  obs_all <- as.numeric(obs_cpue)
  obs_keep <- which(!is.na(obs_all))
  obs_flat <- obs_all[obs_keep]
  n_obs <- length(obs_flat)

  # Subsample draws
  n_total <- nrow(posterior)
  draw_idx <- sample.int(n_total, min(n_sims, n_total))

  # Get parameter column names
  nat_names <- colnames(posterior)

  # Rebuild a data list that exposes movement inputs at the top level so the
  # per-draw recomputation reproduces the fitted dynamics (including movement).
  recompute_data <- model_fit$data
  if (!is.null(model_fit$data$movement)) {
    mv <- model_fit$data$movement
    recompute_data$distance_matrix <- mv$distance_matrix
    recompute_data$attractiveness <- mv$attractiveness
    recompute_data$decay <- mv$decay
    recompute_data$movement_rate <- mv$movement_rate
  }

  # Parameters to vary by draw: the natural-scale model parameters present in
  # both the MLE fit and the posterior (correctly back-transformed, including
  # movement_rate on the logit scale and rho on the atanh scale).
  mle_par <- model_fit$parameters
  varying <- intersect(names(mle_par), nat_names)

  # For each draw, recompute the biomass/CPUE trajectory from that draw's
  # parameters (propagating parameter uncertainty), then add observation error.
  sim_mat <- matrix(NA_real_, nrow = length(draw_idx), ncol = n_obs)

  for (i in seq_along(draw_idx)) {
    d <- draw_idx[i]

    # Per-draw parameter vector: start from the MLE (to supply any fixed or
    # non-sampled parameters) and overwrite with this posterior draw.
    par_d <- mle_par
    par_d[varying] <- posterior[d, varying]

    sigma_obs <- .extract_draw_param(posterior, d, "sigma_obs", nat_names)
    if (is.na(sigma_obs)) {
      sigma_obs <- unname(mle_par[["sigma_obs"]])
    }

    res_d <- tryCatch(
      calculate_model_results(par_d, recompute_data),
      error = function(e) NULL
    )
    if (is.null(res_d)) next

    fit_all <- as.numeric(res_d$fitted_cpue)
    if (length(fit_all) < max(obs_keep)) next
    fit_flat <- fit_all[obs_keep]

    # Simulate observations: log-normal observation error around the draw's
    # own fitted CPUE.
    sim_mat[i, ] <- fit_flat * exp(rnorm(n_obs, 0, sigma_obs))
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
  q_lo <- apply(sim, 2, quantile, probs = 0.025, na.rm = TRUE)
  q_mid <- apply(sim, 2, quantile, probs = 0.50, na.rm = TRUE)
  q_hi <- apply(sim, 2, quantile, probs = 0.975, na.rm = TRUE)

  df <- data.frame(
    index = seq_len(n_obs),
    observed = obs,
    median = q_mid,
    lower = q_lo,
    upper = q_hi
  )

  ggplot(df, aes(x = index)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
      alpha = 0.3,
      fill = "steelblue"
    ) +
    geom_line(aes(y = median), colour = "steelblue", linewidth = 0.8) +
    geom_point(aes(y = observed), colour = "black", size = 2) +
    labs(
      x = "Observation",
      y = "CPUE",
      title = "Posterior Predictive Check",
      subtitle = paste(
        "95% prediction interval | Bayesian p-value:",
        round(x$p_value, 3)
      )
    ) +
    theme_bw()
}


# ---- Helper: extract a named parameter from a posterior draw row ---------
#' @keywords internal
.extract_draw_param <- function(posterior, draw_idx, pname, nat_names) {
  idx <- match(pname, nat_names)
  if (!is.na(idx)) {
    return(posterior[draw_idx, idx])
  }

  NA_real_
}
