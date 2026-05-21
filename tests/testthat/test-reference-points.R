test_that("calculate_reference_points input validation works", {
  # Test non-ProductionModel object
  expect_error(calculate_reference_points("not a model"))
  expect_error(calculate_reference_points(list(fake = "model")))

  # Create unfitted model
  unfitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6)
  )

  expect_error(
    calculate_reference_points(unfitted_model),
    "Model must be fitted"
  )

  # Create fitted model with missing parameters
  fitted_model <- unfitted_model
  fitted_model$fitted <- TRUE
  fitted_model$parameters <- c(r = 0.3, K = 5000) # Missing 'm'

  expect_error(
    calculate_reference_points(fitted_model),
    "Missing required parameters"
  )
})

test_that("calculate_reference_points works for Schaefer model", {
  # Create fitted Schaefer model (m = 2)
  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6),
    parameters = c(r = 0.3, K = 5000, m = 2.0, q = 0.001, sigma_proc = 0.2, sigma_obs = 0.3)
  )
  fitted_model$fitted <- TRUE

  # Add biomass results for current status calculation
  fitted_model$results <- list(
    biomass = c(4000, 3800, 3600, 3400, 3200, 3000),
    harvest_rate = rep(0.25, 6)
  )

  # Calculate reference points
  ref_points <- calculate_reference_points(fitted_model)

  # Check structure
  expect_true(is.list(ref_points))
  expect_true(all(c("msy", "bmsy", "fmsy", "current_status") %in% names(ref_points)))
  expect_s3_class(ref_points, "pt_reference_points")

  # Check Schaefer analytical solutions
  r <- 0.3
  K <- 5000
  expected_msy <- r * K / 4
  expected_bmsy <- K / 2
  expected_fmsy <- expected_msy / expected_bmsy

  expect_equal(ref_points$msy, expected_msy, tolerance = 1e-10)
  expect_equal(ref_points$bmsy, expected_bmsy, tolerance = 1e-10)
  expect_equal(ref_points$fmsy, expected_fmsy, tolerance = 1e-10)

  # Check special case identification
  expect_equal(ref_points$special_case, "Schaefer model (m ~ 2)")

  # Check current status
  expect_true(!is.null(ref_points$current_status))
  expect_equal(ref_points$current_status$current_biomass, 3000)
  expect_equal(ref_points$current_status$b_bmsy_ratio, 3000 / expected_bmsy)
})

test_that("calculate_reference_points works for Fox model", {
  # Create fitted Fox model (m = 1)
  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6),
    parameters = c(r = 0.3, K = 5000, m = 1.0, q = 0.001, sigma_proc = 0.2, sigma_obs = 0.3)
  )
  fitted_model$fitted <- TRUE

  # Calculate reference points
  ref_points <- calculate_reference_points(fitted_model)

  # Check Fox analytical solutions
  r <- 0.3
  K <- 5000
  expected_msy <- r * K / exp(1) # r * K / e
  expected_bmsy <- K / exp(1) # K / e
  expected_fmsy <- expected_msy / expected_bmsy

  expect_equal(ref_points$msy, expected_msy, tolerance = 1e-10)
  expect_equal(ref_points$bmsy, expected_bmsy, tolerance = 1e-10)
  expect_equal(ref_points$fmsy, expected_fmsy, tolerance = 1e-10)

  # Check special case identification
  expect_equal(ref_points$special_case, "Fox model (m ~ 1)")
})

test_that("calculate_reference_points works for general Pella-Tomlinson model", {
  # Create fitted model with m = 3
  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6),
    parameters = c(r = 0.3, K = 5000, m = 3.0, q = 0.001, sigma_proc = 0.2, sigma_obs = 0.3)
  )
  fitted_model$fitted <- TRUE

  # Calculate reference points
  ref_points <- calculate_reference_points(fitted_model)

  # Check general Pella-Tomlinson formulation
  r <- 0.3
  K <- 5000
  m <- 3.0
  expected_msy <- r * K * (m - 1)^((m - 1) / m) / m
  expected_bmsy <- K * (m - 1)^(1 / m) / m
  expected_fmsy <- expected_msy / expected_bmsy

  expect_equal(ref_points$msy, expected_msy, tolerance = 1e-10)
  expect_equal(ref_points$bmsy, expected_bmsy, tolerance = 1e-10)
  expect_equal(ref_points$fmsy, expected_fmsy, tolerance = 1e-10)

  # Check special case identification
  expect_equal(ref_points$special_case, "General Pella-Tomlinson model")
})

test_that("calculate_reference_points handles edge cases", {
  # Test very low m value (edge case)
  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6),
    parameters = c(r = 0.3, K = 5000, m = 0.5, q = 0.001, sigma_proc = 0.2, sigma_obs = 0.3)
  )
  fitted_model$fitted <- TRUE

  # Should give warning and use Schaefer approximation
  expect_warning(ref_points <- calculate_reference_points(fitted_model))

  # Should fall back to Schaefer values
  expect_equal(ref_points$msy, 0.3 * 5000 / 4)
  expect_equal(ref_points$bmsy, 5000 / 2)

  # Test negative parameter values
  fitted_model$parameters["r"] <- -0.3
  expect_error(calculate_reference_points(fitted_model), "must be positive")
})

test_that("calculate_reference_points current status calculations work", {
  # Create model with declining biomass
  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1500, 6), # High catch
    cpue = rep(1.0, 6),
    effort = rep(1500, 6),
    parameters = c(r = 0.3, K = 5000, m = 2.0, q = 0.001, sigma_proc = 0.2, sigma_obs = 0.3)
  )
  fitted_model$fitted <- TRUE

  # Add biomass trajectory showing decline to below BMSY
  bmsy <- 2500 # K/2 for Schaefer
  fitted_model$results <- list(
    biomass = c(3000, 2800, 2600, 2400, 2200, 1800), # Declining to below BMSY
    harvest_rate = c(0.5, 0.54, 0.58, 0.62, 0.68, 0.83) # Increasing harvest rate
  )

  ref_points <- calculate_reference_points(fitted_model)

  # Check current status calculations
  expect_equal(ref_points$current_status$current_biomass, 1800)
  expect_equal(ref_points$current_status$b_bmsy_ratio, 1800 / 2500)
  expect_equal(ref_points$current_status$status, "Below BMSY but above half BMSY")

  # Check harvest rate status
  expect_equal(ref_points$current_status$current_harvest_rate, 0.83)
  fmsy <- 375 / 2500 # MSY / BMSY
  expect_equal(ref_points$current_status$f_fmsy_ratio, 0.83 / fmsy)
  expect_equal(ref_points$current_status$harvest_status, "Overfishing occurring")

  # Test overfished status
  fitted_model$results$biomass[6] <- 1000 # Below half BMSY
  ref_points2 <- calculate_reference_points(fitted_model)
  expect_equal(ref_points2$current_status$status, "Below half BMSY (overfished)")

  # Test above BMSY status
  fitted_model$results$biomass[6] <- 3000 # Above BMSY
  fitted_model$results$harvest_rate[6] <- 0.1 # Low harvest rate
  ref_points3 <- calculate_reference_points(fitted_model)
  expect_equal(ref_points3$current_status$status, "Above BMSY")
  expect_equal(ref_points3$current_status$harvest_status, "No overfishing")
})

test_that("calculate_reference_points calculates user-defined biomass targets analytically", {
  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6),
    parameters = c(r = 0.3, K = 5000, m = 2.0, B0 = 6000, q = 0.001, sigma_proc = 0.2, sigma_obs = 0.3)
  )
  fitted_model$fitted <- TRUE

  ref_points <- calculate_reference_points(fitted_model, biomass_target = c(0.4, 0.35))

  expect_true("target_reference_points" %in% names(ref_points))
  expect_equal(nrow(ref_points$target_reference_points), 2)
  expect_equal(ref_points$target_reference_points$baseline, c("B0", "B0"))
  expect_equal(ref_points$target_reference_points$biomass_name, c("B_40%B0", "B_35%B0"))
  expect_equal(ref_points$target_reference_points$f_name, c("F_40%B0", "F_35%B0"))

  expected_b40 <- 0.4 * 6000
  expected_f40 <- 0.3 / 2 * (1 - expected_b40 / 5000)
  expect_equal(ref_points$target_reference_points$biomass[1], expected_b40, tolerance = 1e-10)
  expect_equal(ref_points$target_reference_points$fishing_mortality[1], expected_f40, tolerance = 1e-10)
})

test_that("calculate_reference_points can use K as the depletion baseline", {
  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6),
    parameters = c(r = 0.3, K = 5000, m = 2.0, B0 = 6000, q = 0.001, sigma_proc = 0.2, sigma_obs = 0.3)
  )
  fitted_model$fitted <- TRUE

  ref_points <- calculate_reference_points(fitted_model, biomass_target = 0.4, baseline = "K")

  expect_equal(ref_points$target_reference_points$baseline, "K")
  expect_equal(ref_points$target_reference_points$biomass_name, "B_40%K")
  expect_equal(ref_points$target_reference_points$f_name, "F_40%K")
  expect_equal(ref_points$target_reference_points$biomass, 0.4 * 5000, tolerance = 1e-10)
  expect_equal(ref_points$target_reference_points$fishing_mortality, 0.3 / 2 * (1 - 0.4), tolerance = 1e-10)
})

test_that("calculate_reference_points handles Fox depletion targets analytically", {
  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6),
    parameters = c(r = 0.3, K = 5000, m = 1.0, B0 = 6000, q = 0.001, sigma_proc = 0.2, sigma_obs = 0.3)
  )
  fitted_model$fitted <- TRUE

  ref_points <- calculate_reference_points(fitted_model, biomass_target = 0.4, baseline = "K")
  expected_biomass <- 0.4 * 5000
  expected_f <- 0.3 * log(5000 / expected_biomass)

  expect_equal(ref_points$target_reference_points$biomass, expected_biomass, tolerance = 1e-10)
  expect_equal(ref_points$target_reference_points$fishing_mortality, expected_f, tolerance = 1e-10)
})

test_that("calculate_reference_points validates biomass target inputs", {
  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6),
    parameters = c(r = 0.3, K = 5000, m = 2.0, q = 0.001, sigma_proc = 0.2, sigma_obs = 0.3)
  )
  fitted_model$fitted <- TRUE

  expect_error(calculate_reference_points(fitted_model, biomass_target = 0), "must be > 0 and <= 1")
  expect_error(calculate_reference_points(fitted_model, biomass_target = 1.2), "must be > 0 and <= 1")
  expect_error(calculate_reference_points(fitted_model, biomass_target = "0.4"), "must be a numeric vector")
  expect_error(calculate_reference_points(fitted_model, biomass_target = 0.4, baseline = "B0"), "no fitted B0 parameter")
})

test_that("print method for reference points works", {
  # Create reference points object
  ref_points <- list(
    msy = 375,
    bmsy = 2500,
    fmsy = 0.15,
    current_status = list(
      current_biomass = 2200,
      b_bmsy_ratio = 0.88,
      status = "Below BMSY but above half BMSY",
      current_harvest_rate = 0.18,
      f_fmsy_ratio = 1.2,
      harvest_status = "Overfishing occurring"
    ),
    model_type = "Pella-Tomlinson",
    shape_parameter = 2.0,
    special_case = "Schaefer model (m ≈ 2)",
    calculation_date = Sys.time()
  )
  class(ref_points) <- c("pt_reference_points", "list")

  # Test print method doesn't error
  expect_output(print(ref_points), "Pella-Tomlinson Reference Points")
  expect_output(print(ref_points), "MSY:")
  expect_output(print(ref_points), "BMSY:")
  expect_output(print(ref_points), "Current Stock Status:")
})

test_that("estimate_biomass function works correctly", {
  # Create fitted model with biomass results
  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6)
  )
  fitted_model$fitted <- TRUE
  fitted_model$results <- list(
    biomass = c(4000, 3800, 3600, 3400, 3200, 3000),
    harvest_rate = c(0.25, 0.26, 0.28, 0.29, 0.31, 0.33)
  )

  # Test full extraction
  biomass_est <- estimate_biomass(fitted_model)

  expect_true(is.data.frame(biomass_est))
  expect_true(all(c("year", "biomass") %in% names(biomass_est)))
  expect_equal(nrow(biomass_est), 6)
  expect_equal(biomass_est$year, 2000:2005)
  expect_equal(biomass_est$biomass, c(4000, 3800, 3600, 3400, 3200, 3000))

  # Should include harvest rates
  expect_true("harvest_rate" %in% names(biomass_est))

  # Test year subset
  biomass_subset <- estimate_biomass(fitted_model, years = 2002:2004)
  expect_equal(nrow(biomass_subset), 3)
  expect_equal(biomass_subset$year, 2002:2004)
  expect_equal(biomass_subset$biomass, c(3600, 3400, 3200))

  # Test error conditions
  unfitted_model <- fitted_model
  unfitted_model$fitted <- FALSE
  expect_error(estimate_biomass(unfitted_model))

  # Test no biomass results
  no_biomass_model <- fitted_model
  no_biomass_model$results <- list(other = "stuff")
  expect_error(estimate_biomass(no_biomass_model))

  # Test years not available
  expect_error(estimate_biomass(fitted_model, years = 1990:1995))

  # Test partial year availability
  expect_warning(biomass_partial <- estimate_biomass(fitted_model, years = 1999:2002))
  expect_equal(biomass_partial$year, 2000:2002)
})
