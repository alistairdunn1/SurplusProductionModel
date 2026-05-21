# Test Data Validation Functions
# Tests for CPUE, catch, and parameter validation functions

test_that("validate_cpue_data works with valid data", {
  # Create valid test data
  cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
  years <- 2000:2005
  effort <- c(500, 550, 600, 450, 500, 650)
  catch <- cpue * effort

  # Basic validation should pass
  expect_true(validate_cpue_data(cpue))

  # Detailed validation
  result <- validate_cpue_data(cpue,
    years = years, effort = effort, catch = catch,
    return_details = TRUE
  )
  expect_true(result$valid)
  expect_equal(result$summary$data_quality, "PASS")
  expect_equal(result$quality_metrics$n_observations, 6)
  expect_equal(result$quality_metrics$n_valid, 6)
  expect_equal(result$quality_metrics$n_missing, 0)
})

test_that("validate_cpue_data catches invalid data", {
  # Test with negative values
  expect_error(
    validate_cpue_data(c(-1, 2, 3)),
    "CPUE data must be positive"
  )

  # Test with zero values
  expect_error(
    validate_cpue_data(c(0, 2, 3)),
    "CPUE data must be positive"
  )

  # Test with missing values
  expect_error(
    validate_cpue_data(c(NA, 2, 3)),
    "CPUE data contains missing values"
  )

  # Test with infinite values
  expect_error(
    validate_cpue_data(c(Inf, 2, 3)),
    "CPUE data contains infinite values"
  )

  # Test with non-numeric data
  expect_error(
    validate_cpue_data(c("a", "b", "c")),
    "CPUE data must be numeric"
  )

  # Test with empty data
  expect_error(
    validate_cpue_data(numeric(0)),
    "CPUE data cannot be empty"
  )
})

test_that("validate_cpue_data handles mismatched lengths", {
  cpue <- c(2.1, 2.0, 1.9)
  years <- 2000:2005 # Different length
  effort <- c(500, 550) # Different length

  expect_error(
    validate_cpue_data(cpue, years = years),
    "Years length.*must match CPUE length"
  )

  expect_error(
    validate_cpue_data(cpue, effort = effort),
    "Effort length.*must match CPUE length"
  )
})

test_that("validate_cpue_data detects outliers", {
  # Data with clear outlier
  cpue <- c(2.0, 2.1, 1.9, 20.0, 2.0, 1.8) # 20.0 is clear outlier

  result <- validate_cpue_data(cpue,
    outlier_method = "zscore", outlier_threshold = 2.0,
    return_details = TRUE
  )
  expect_true(length(result$outliers) > 0)
  expect_true(4 %in% result$outliers) # Position of outlier

  # IQR method
  result_iqr <- validate_cpue_data(cpue,
    outlier_method = "iqr", outlier_threshold = 1.5,
    return_details = TRUE
  )
  expect_true(length(result_iqr$outliers) > 0)
})

test_that("validate_cpue_data validates cross-relationships", {
  cpue <- c(2.0, 2.0, 2.0)
  effort <- c(500, 500, 500)
  catch <- c(1000, 1000, 1100) # Last value inconsistent

  # Should generate warning about inconsistency
  expect_warning(
    validate_cpue_data(cpue, effort = effort, catch = catch),
    "CPUE values inconsistent with catch/effort"
  )
})

test_that("validate_catch_data works with valid data", {
  # Create valid test data
  catch <- c(1000, 1100, 1200, 950, 1050, 1300)
  years <- 2000:2005
  cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
  effort <- catch / cpue

  # Basic validation should pass
  expect_true(validate_catch_data(catch))

  # Detailed validation
  result <- validate_catch_data(catch,
    years = years, cpue = cpue, effort = effort,
    return_details = TRUE
  )
  expect_true(result$valid)
  expect_equal(result$summary$data_quality, "PASS")
  expect_equal(result$quality_metrics$n_observations, 6)
  expect_equal(result$quality_metrics$n_valid, 6)
  expect_true(result$quality_metrics$total_catch > 0)
})

test_that("validate_catch_data catches invalid data", {
  # Test with negative values
  expect_error(
    validate_catch_data(c(-100, 1000, 1200)),
    "Catch data must be non-negative"
  )

  # Test with missing values
  expect_error(
    validate_catch_data(c(NA, 1000, 1200)),
    "Catch data contains missing values"
  )

  # Test with infinite values
  expect_error(
    validate_catch_data(c(Inf, 1000, 1200)),
    "Catch data contains infinite values"
  )

  # Test with non-numeric data
  expect_error(
    validate_catch_data(c("high", "medium", "low")),
    "Catch data must be numeric"
  )

  # Test with empty data
  expect_error(
    validate_catch_data(numeric(0)),
    "Catch data cannot be empty"
  )
})

test_that("validate_catch_data detects extreme variations", {
  # Data with extreme year-to-year change
  catch <- c(1000, 1100, 10000, 1200) # Large jump to 10000

  result <- validate_catch_data(catch, return_details = TRUE)
  expect_true(length(result$warnings) > 0)
  expect_true(any(grepl("Extreme year-to-year", result$warnings)))
})

test_that("validate_catch_data validates cross-relationships", {
  catch <- c(1000, 1000, 1000)
  cpue <- c(2.0, 2.0, 2.0)
  effort <- c(500, 500, 550) # Last value inconsistent

  # Should generate warning about inconsistency
  expect_warning(
    validate_catch_data(catch, cpue = cpue, effort = effort),
    "Catch values inconsistent with CPUE\\*effort"
  )
})

test_that("validate_pt_parameters works with valid parameters", {
  # Valid default parameters
  params <- c(
    r = 0.1, K = 50000, m = 2.0, q = 0.001,
    sigma_proc = 0.1, sigma_obs = 0.2
  )

  expect_true(validate_pt_parameters(params))

  # Detailed validation
  result <- validate_pt_parameters(params, return_details = TRUE)
  expect_true(result$valid)
  expect_equal(result$summary$parameter_check, "PASS")
  expect_true(length(result$biological_metrics) > 0)
  expect_true(result$biological_metrics$msy > 0)
  expect_true(result$biological_metrics$bmsy > 0)
  expect_true(result$biological_metrics$bmsy < params["K"])
})

test_that("validate_pt_parameters catches invalid parameters", {
  # Test missing names
  expect_error(
    validate_pt_parameters(c(0.1, 50000, 2.0, 0.001, 0.1, 0.2)),
    "Parameters must have names"
  )

  # Test missing required parameters
  expect_error(
    validate_pt_parameters(c(r = 0.1, K = 50000)),
    "Missing required parameters"
  )

  # Test negative parameters
  expect_error(
    validate_pt_parameters(c(
      r = -0.1, K = 50000, m = 2.0, q = 0.001,
      sigma_proc = 0.1, sigma_obs = 0.2
    )),
    "Parameter r must be positive"
  )

  # Test infinite parameters
  expect_error(
    validate_pt_parameters(c(
      r = Inf, K = 50000, m = 2.0, q = 0.001,
      sigma_proc = 0.1, sigma_obs = 0.2
    )),
    "Parameter r must be finite"
  )

  # Test zero parameters
  expect_error(
    validate_pt_parameters(c(
      r = 0, K = 50000, m = 2.0, q = 0.001,
      sigma_proc = 0.1, sigma_obs = 0.2
    )),
    "Parameter r must be positive"
  )
})

test_that("validate_pt_parameters warns about extreme values", {
  # Parameters outside recommended bounds
  params <- c(
    r = 5.0, K = 50000, m = 2.0, q = 0.001, # r too high
    sigma_proc = 0.1, sigma_obs = 0.2
  )

  expect_warning(
    validate_pt_parameters(params),
    "Parameter r is above recommended maximum"
  )

  # Very low K
  params2 <- c(
    r = 0.1, K = 100, m = 2.0, q = 0.001, # K too low
    sigma_proc = 0.1, sigma_obs = 0.2
  )

  expect_warning(
    validate_pt_parameters(params2),
    "Parameter K is below recommended minimum"
  )
})

test_that("validate_pt_parameters calculates biological metrics correctly", {
  # Schaefer model (m = 2)
  params_schaefer <- c(
    r = 0.2, K = 10000, m = 2.0, q = 0.001,
    sigma_proc = 0.1, sigma_obs = 0.2
  )

  result <- validate_pt_parameters(params_schaefer, return_details = TRUE)

  # Unified canonical Schaefer special-case formulas (m = 2)
  expected_msy <- 0.2 * 10000 / 4
  expected_bmsy <- 10000 / 2

  expect_equal(result$biological_metrics$msy, expected_msy, tolerance = 0.1)
  expect_equal(result$biological_metrics$bmsy, expected_bmsy, tolerance = 0.1)
  expect_equal(result$biological_metrics$bmsy_k_ratio, 1 / 2, tolerance = 0.01)
})

test_that("validate_pt_parameters handles Fox model (m ≈ 1)", {
  # Fox-like model parameters (m close to 1)
  params_fox <- c(
    r = 0.2, K = 10000, m = 1.1, q = 0.001,
    sigma_proc = 0.1, sigma_obs = 0.2
  )

  result <- validate_pt_parameters(params_fox, return_details = TRUE)

  expect_true(result$valid)
  expect_true(result$biological_metrics$msy > 0)
  expect_true(result$biological_metrics$bmsy > 0)
  expect_true(result$biological_metrics$bmsy < params_fox["K"])
})

test_that("validate_pt_parameters warns about unusual shape parameters", {
  # Very low m should be valid but warn as biologically unusual
  params_invalid_m <- c(
    r = 0.1, K = 50000, m = 0.4, q = 0.001,
    sigma_proc = 0.1, sigma_obs = 0.2
  )

  result <- validate_pt_parameters(params_invalid_m, return_details = TRUE)
  expect_true(result$valid)
  expect_true(any(grepl("Shape parameter m < 0.5", result$warnings)))
  expect_true(result$biological_metrics$msy > 0)
  expect_true(result$biological_metrics$bmsy > 0)

  # Test very high m (should be valid but generate warning)
  params_high_m <- c(
    r = 0.1, K = 50000, m = 7.0, q = 0.001,
    sigma_proc = 0.1, sigma_obs = 0.2
  )

  expect_warning(
    validate_pt_parameters(params_high_m),
    "Shape parameter m > 5 may indicate unrealistic production dynamics"
  )
})

test_that("detect_outliers function works correctly", {
  # Test data with clear outliers
  x <- c(1, 2, 3, 2, 1, 10, 3, 2) # 10 is outlier

  # Z-score method
  outliers_z <- SurplusProductionModel:::detect_outliers(x, method = "zscore", threshold = 2)
  expect_true(6 %in% outliers_z) # Position of outlier

  # IQR method
  outliers_iqr <- SurplusProductionModel:::detect_outliers(x, method = "iqr", threshold = 1.5)
  expect_true(6 %in% outliers_iqr) # Position of outlier

  # Test with no outliers
  x_normal <- c(1, 2, 3, 2, 1, 3, 2)
  outliers_none <- SurplusProductionModel:::detect_outliers(x_normal, method = "zscore", threshold = 3)
  expect_equal(length(outliers_none), 0)

  # Test with insufficient data
  x_small <- c(1, 2)
  outliers_small <- SurplusProductionModel:::detect_outliers(x_small, method = "zscore")
  expect_equal(length(outliers_small), 0)
})

test_that("validation functions handle edge cases", {
  # Single observation
  expect_true(validate_cpue_data(2.5))
  expect_true(validate_catch_data(1000))

  # Very large datasets (performance test)
  large_cpue <- rep(2.0, 1000)
  large_catch <- rep(1000, 1000)

  expect_true(validate_cpue_data(large_cpue))
  expect_true(validate_catch_data(large_catch))

  # Zero catch (valid for some years)
  catch_with_zero <- c(0, 1000, 1200)
  expect_true(validate_catch_data(catch_with_zero))

  result <- validate_catch_data(catch_with_zero, return_details = TRUE)
  expect_equal(result$quality_metrics$n_zero_catch, 1)
})
