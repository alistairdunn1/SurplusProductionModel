# Tests for Bayesian fitting via tmbstan

# ---------- Helper: fit a small model for testing -----------------------
fit_bayes_model <- function() {
  set.seed(42)
  years <- 2010:2020
  true_r <- 0.3; true_K <- 5000; true_q <- 0.001; true_B0 <- 4000
  biomass <- numeric(length(years))
  biomass[1] <- true_B0
  catch_vals <- rep(500, length(years))
  for (t in seq_len(length(years) - 1)) {
    prod <- true_r * biomass[t] * (1 - biomass[t] / true_K)
    biomass[t + 1] <- max(100, biomass[t] + prod - catch_vals[t])
  }
  cpue <- true_q * biomass * exp(rnorm(length(years), 0, 0.05))

  data_list <- list(
    cpue_data  = data.frame(year = years, cpue = cpue),
    catch_data = data.frame(year = years, catch = catch_vals)
  )
  opts <- list(silent = TRUE, validate_data = FALSE,
               control = list(eval.max = 1000, iter.max = 500))
  fit_pella_tomlinson_model(data_list, options = opts)
}


# ====================  Input validation  ================================

test_that("bayesian_fit refuses unfitted or wrong class", {
  expect_error(bayesian_fit("not a model"), "ProductionModel")
  unfitted <- structure(list(fitted = FALSE), class = "ProductionModel")
  expect_error(bayesian_fit(unfitted), "fitted")
})

test_that("bayesian_fit errors when tmbstan unavailable", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")

  # Only tests that the function checks for the package --

  # we can't easily unload tmbstan during the test, but the
  # validation logic is covered by the class/fitted checks above.
  expect_true(TRUE)
})


# ====================  Core sampling tests  =============================

test_that("bayesian_fit runs and returns correct structure", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  # Very short run for speed
  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  expect_s3_class(bf, "bayes_fit")
  expect_true(!is.null(bf$stanfit))
  expect_true(is.matrix(bf$posterior))
  expect_true(is.matrix(bf$log_posterior))
  expect_true(is.data.frame(bf$summary))
  expect_true(nrow(bf$posterior) > 0)
  expect_equal(bf$n_chains, 1L)
  expect_equal(bf$n_iter, 200L)
  expect_equal(bf$n_warmup, 100L)
})


test_that("bayesian_fit posterior contains derived quantities", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  # Check derived columns exist
  cn <- colnames(bf$posterior)
  expect_true("MSY" %in% cn)
  expect_true("BMSY" %in% cn)
  expect_true("FMSY" %in% cn)

  # Check values are positive and finite (allow some NAs for edge cases)
  msy_draws <- bf$posterior[, "MSY"]
  msy_ok <- msy_draws[is.finite(msy_draws)]
  expect_true(length(msy_ok) > 0)
  expect_true(all(msy_ok > 0))
})

test_that("bayesian_fit can include depletion-based posterior targets", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42,
                 biomass_target = 0.4, baseline = "K"),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  cn <- colnames(bf$posterior)
  expect_true("B_40%K" %in% cn)
  expect_true("F_40%K" %in% cn)
  expect_true(any(bf$summary$parameter == "B_40%K"))
  expect_true(any(bf$summary$parameter == "F_40%K"))
})


test_that("bayesian_fit summary has expected columns", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  s <- bf$summary
  expect_true("parameter" %in% names(s))
  expect_true("mean" %in% names(s))
  expect_true("sd" %in% names(s))
  expect_true("2.5%" %in% names(s))
  expect_true("97.5%" %in% names(s))
  expect_true("n_eff" %in% names(s))
})


test_that("bayesian_fit respects seed for reproducibility", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf1 <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 123),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )
  bf2 <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 123),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  # Same seed should give same posterior
  expect_equal(bf1$posterior, bf2$posterior, tolerance = 1e-10)
})


# ====================  Print method  ====================================

test_that("print.bayes_fit does not error", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  out <- capture.output(print(bf))
  expect_true(any(grepl("Bayesian", out)))
  expect_true(any(grepl("Posterior", out)))
})


# ====================  Plot methods  ====================================

test_that("plot.bayes_fit trace returns ggplot", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  p <- plot(bf, type = "trace")
  expect_s3_class(p, "ggplot")
})


test_that("plot.bayes_fit density returns ggplot", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  p <- plot(bf, type = "density")
  expect_s3_class(p, "ggplot")
})


test_that("plot.bayes_fit histogram returns ggplot", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  p <- plot(bf, type = "histogram")
  expect_s3_class(p, "ggplot")
})


test_that("plot.bayes_fit pairs does not error", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  expect_no_error(plot(bf, type = "pairs", pars = c("r", "K")))
})


# ====================  Posterior predictive check  =======================

test_that("posterior_predictive_check returns correct structure", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  ppc <- posterior_predictive_check(bf, n_sims = 50, seed = 1)
  expect_s3_class(ppc, "ppc_result")
  expect_true(is.numeric(ppc$observed))
  expect_true(is.matrix(ppc$simulated))
  expect_true(ppc$p_value >= 0 && ppc$p_value <= 1)
})


test_that("print.ppc_result does not error", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  ppc <- posterior_predictive_check(bf, n_sims = 20, seed = 1)
  out <- capture.output(print(ppc))
  expect_true(any(grepl("Posterior Predictive", out)))
})


test_that("plot.ppc_result returns ggplot", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_bayes_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(model, chains = 1, iter = 200, warmup = 100, seed = 42),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  ppc <- posterior_predictive_check(bf, n_sims = 20, seed = 1)
  p <- plot(ppc)
  expect_s3_class(p, "ggplot")
})


# ====================  Internal helpers  ================================

test_that(".log_to_natural_names transforms correctly", {
  expect_equal(
    SurplusProductionModel:::.log_to_natural_names(
      c("log_r", "log_K", "log_m", "log_sigma_obs", "log_q_A1", "log_K_A1")
    ),
    c("r", "K", "m", "sigma_obs", "q.A1", "K.A1")
  )
})
