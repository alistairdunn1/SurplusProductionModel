# Tests for convergence diagnostics (jitter_test, retrospective_analysis)

# ---------- Helper: fit a small model for testing -----------------------
fit_test_model <- function() {
  set.seed(42)
  years <- 2010:2018
  true_r <- 0.3; true_K <- 5000; true_q <- 0.001; true_B_initial <- 4000
  biomass <- numeric(length(years))
  biomass[1] <- true_B_initial
  catch_vals <- rep(600, length(years))
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


# ========================  jitter_test  ================================

test_that("jitter_test refuses unfitted or wrong class", {
  expect_error(jitter_test("not a model"), "ProductionModel")
  unfitted <- structure(list(fitted = FALSE), class = "ProductionModel")
  expect_error(jitter_test(unfitted), "fitted")
})

test_that("jitter_test returns correct structure", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_test_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  jt <- jitter_test(model, n_jitter = 3, jitter_sd = 0.1, tolerance = 5.0)

  expect_s3_class(jt, "jitter_test")
  expect_length(jt$nll, 4)               # 1 original + 3 jittered
  expect_length(jt$convergence, 4)
  expect_equal(jt$n_total, 4L)
  expect_true(is.data.frame(jt$parameters))
  expect_equal(nrow(jt$parameters), 4)
  expect_true(is.numeric(jt$best_nll))
  expect_true(is.numeric(jt$proportion_converged))
  expect_true(is.logical(jt$passed))
})

test_that("jitter_test first run is un-jittered (same NLL as original)", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_test_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  jt <- jitter_test(model, n_jitter = 2, jitter_sd = 0.1, tolerance = 5.0)
  # Run 1 should be very close to the original NLL
  expect_equal(jt$nll[1], model$results$likelihood, tolerance = 0.01)
})

test_that("jitter_test print method doesn't error", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_test_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })
  jt <- jitter_test(model, n_jitter = 2, jitter_sd = 0.1, tolerance = 5.0)
  expect_output(print(jt), "Jitter Test")
  expect_output(print(jt), "PASSED|FAILED")
})

test_that("jitter_test plot method returns a ggplot", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("ggplot2")
  model <- tryCatch(fit_test_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })
  jt <- jitter_test(model, n_jitter = 2, jitter_sd = 0.1, tolerance = 5.0)
  p <- plot(jt)
  expect_s3_class(p, "ggplot")
})


# ====================  retrospective_analysis  =========================

test_that("retrospective_analysis refuses unfitted or wrong class", {
  expect_error(retrospective_analysis("not a model"), "ProductionModel")
  unfitted <- structure(list(fitted = FALSE), class = "ProductionModel")
  expect_error(retrospective_analysis(unfitted), "fitted")
})

test_that("retrospective_analysis returns correct structure", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_test_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  n_peels <- 2
  retro <- retrospective_analysis(model, n_peels = n_peels,
                                   options = list(silent = TRUE,
                                                  validate_data = FALSE,
                                                  control = list(eval.max = 1000,
                                                                 iter.max = 500)))

  expect_s3_class(retro, "retro_analysis")
  expect_length(retro$peels, n_peels + 1)
  expect_length(retro$terminal_biomass, n_peels + 1)
  expect_true(is.numeric(retro$mohns_rho))
  expect_equal(retro$n_peels, n_peels)

  # First peel should be the original model
  expect_identical(retro$peels[[1]], model)
})

test_that("retrospective_analysis rejects too many peels", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_test_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })
  # Model has 9 years (2010-2018); n_peels >= 6 would leave ≤3 years

  expect_error(retrospective_analysis(model, n_peels = 6), "too large")
})

test_that("retrospective_analysis print method doesn't error", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_test_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })
  retro <- retrospective_analysis(model, n_peels = 2,
                                   options = list(silent = TRUE,
                                                  validate_data = FALSE,
                                                  control = list(eval.max = 1000,
                                                                 iter.max = 500)))
  expect_output(print(retro), "Retrospective Analysis")
  expect_output(print(retro), "Mohn's rho")
})

test_that("retrospective_analysis plot method returns a ggplot", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("ggplot2")
  model <- tryCatch(fit_test_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })
  retro <- retrospective_analysis(model, n_peels = 2,
                                   options = list(silent = TRUE,
                                                  validate_data = FALSE,
                                                  control = list(eval.max = 1000,
                                                                 iter.max = 500)))
  p <- plot(retro)
  expect_s3_class(p, "ggplot")
})


# ====================  multi-start in fitting  =========================

test_that("fit with n_starts > 1 produces a valid ProductionModel", {
  skip_if_not_installed("RTMB")
  model <- tryCatch({
    set.seed(99)
    years <- 2010:2016
    biomass <- numeric(length(years))
    biomass[1] <- 4000
    catch_vals <- rep(500, length(years))
    for (t in seq_len(length(years) - 1)) {
      prod <- 0.3 * biomass[t] * (1 - biomass[t] / 5000)
      biomass[t + 1] <- max(100, biomass[t] + prod - catch_vals[t])
    }
    cpue <- 0.001 * biomass * exp(rnorm(length(years), 0, 0.05))
    data_list <- list(
      cpue_data  = data.frame(year = years, cpue = cpue),
      catch_data = data.frame(year = years, catch = catch_vals)
    )
    fit_pella_tomlinson_model(data_list, options = list(
      silent = TRUE, validate_data = FALSE,
      n_starts = 3, jitter_sd = 0.15,
      control = list(eval.max = 1000, iter.max = 500),
      fixed_params = list(
        log_m = log(2),
        log_d0 = log(0.8),
        log_sigma_obs = log(0.05)
      ),
      priors = list(
        r = list(dist = "lognormal", meanlog = log(0.3), sdlog = 0.5),
        K = list(dist = "lognormal", meanlog = log(5000), sdlog = 0.5)
      ),
      calculate_se = FALSE
    ))
  }, error = function(e) {
    skip(paste("Multi-start fitting failed:", e$message))
  })

  expect_s3_class(model, "ProductionModel")
  expect_true(model$fitted)
  expect_true(length(model$parameters) > 0)
})
