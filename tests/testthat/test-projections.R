make_projection_test_model <- function() {
  years <- 2000:2014
  catch <- rep(350, length(years))
  cpue <- seq(1.8, 1.2, length.out = length(years))

  model <- SurplusProductionModel::ProductionModel(
    years = years,
    catch = catch,
    cpue = cpue,
    effort = catch / pmax(cpue, 1e-8),
    parameters = c(r = 0.35, K = 6000, m = 2, q = 3e-4, sigma_proc = 0.15, sigma_obs = 0.25, d0 = 5000 / 6000)
  )

  res <- SurplusProductionModel:::calculate_model_results(model$parameters, model$data)
  model$results <- list(
    biomass = res$biomass,
    harvest_rate = res$harvest_rate,
    fitted_cpue = res$fitted_cpue,
    residuals = res$residuals,
    bmsy = 3000,
    fmsy = 0.175
  )
  model$fitted <- TRUE
  model
}

test_that("project_forward validates input", {
  expect_error(project_forward("not-a-model", control = 0.1), "ProductionModel")

  model <- make_projection_test_model()
  model$fitted <- FALSE
  expect_error(project_forward(model, control = 0.1), "fitted")
})

test_that("project_forward supports F control with historical process sampling", {
  model <- make_projection_test_model()

  proj <- project_forward(
    model,
    horizon = 8,
    control = 0.12,
    control_type = "F",
    process_error = "historical",
    historical_window = 10,
    n_sim = 200,
    seed = 42
  )

  expect_s3_class(proj, "pt_projection")
  expect_true(all(c("summary", "simulations", "settings") %in% names(proj)))
  expect_equal(nrow(proj$summary), 8)
  expect_equal(dim(proj$simulations$biomass), c(9, 1, 200))
  expect_equal(length(proj$control), 8)
  expect_true(all(proj$summary$prob_B_above_BMSY >= 0 & proj$summary$prob_B_above_BMSY <= 1, na.rm = TRUE))
})

test_that("project_forward supports catch control with vector control path", {
  model <- make_projection_test_model()
  catches <- seq(300, 500, length.out = 6)

  proj <- project_forward(
    model,
    horizon = 6,
    control = catches,
    control_type = "catch",
    process_error = "normal",
    n_sim = 150,
    seed = 99
  )

  expect_equal(proj$control, catches)
  expect_equal(nrow(proj$summary), 6)
  expect_true(all(is.finite(proj$summary$biomass_q50)))
  expect_true(all(is.finite(proj$summary$harvest_q50)))
})

test_that("project_forward with no process error is deterministic across sims", {
  model <- make_projection_test_model()

  proj <- project_forward(
    model,
    horizon = 5,
    control = 0.1,
    control_type = "F",
    process_error = "none",
    n_sim = 40,
    seed = 7
  )

  b <- proj$simulations$biomass[, 1, ]
  b_sd <- apply(b, 1, stats::sd)
  expect_true(all(b_sd < 1e-12))
})

test_that("print.pt_projection runs without error", {
  model <- make_projection_test_model()
  proj <- project_forward(model, horizon = 3, control = 0.12, n_sim = 20, seed = 11)
  expect_no_error(print(proj))
})
