# Test S4 Classes and Constructors
# Tests for ProductionModel class definition, validation, and constructor functions

test_that("ProductionModel class can be created with valid data", {
  # Create test data
  years <- 2000:2010
  catch <- c(1000, 1100, 1200, 950, 1050, 1300, 1150, 1250, 1400, 1350, 1200)
  cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8, 2.0, 1.9, 1.7, 1.8, 1.9)
  effort <- catch / cpue

  # Test constructor
  model <- ProductionModel(years = years, catch = catch, cpue = cpue, effort = effort)

  # Basic checks
  expect_s4_class(model, "ProductionModel")
  expect_true(validObject(model))
  expect_false(fitted(model))
  expect_equal(model@model_type, "pella_tomlinson")
  expect_true(inherits(model@creation_date, "POSIXct"))
})

test_that("ProductionModel validates data consistency", {
  years <- 2000:2005
  catch <- c(1000, 1100, 1200, 950, 1050, 1300)
  cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
  effort <- catch / cpue

  # Test mismatched data lengths
  expect_error(
    ProductionModel(years = years[1:5], catch = catch, cpue = cpue, effort = effort),
    "catch data length.*must match years length"
  )

  expect_error(
    ProductionModel(years = years, catch = catch[1:5], cpue = cpue, effort = effort),
    "catch data length.*must match years length"
  )

  expect_error(
    ProductionModel(years = years, catch = catch, cpue = cpue[1:5], effort = effort),
    "cpue data length.*must match years length"
  )

  expect_error(
    ProductionModel(years = years, catch = catch, cpue = cpue, effort = effort[1:5]),
    "effort data length.*must match years length"
  )
})

test_that("ProductionModel validates data types and ranges", {
  years <- 2000:2005
  catch <- c(1000, 1100, 1200, 950, 1050, 1300)
  cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
  effort <- catch / cpue

  # Test invalid years (non-integer)
  expect_error(
    ProductionModel(
      years = c(2000.5, 2001.5, 2002.5, 2003.5, 2004.5, 2005.5),
      catch = catch, cpue = cpue, effort = effort
    ),
    "Must be of type 'integerish'"
  )

  # Test non-increasing years
  expect_error(
    ProductionModel(
      years = c(2000, 2001, 2000, 2003, 2004, 2005),
      catch = catch, cpue = cpue, effort = effort
    ),
    "years must be in strictly increasing order"
  )

  # Test negative catch
  expect_error(
    ProductionModel(
      years = years, catch = c(-100, 1100, 1200, 950, 1050, 1300),
      cpue = cpue, effort = effort
    ),
    "Element 1 is not >= 0"
  )

  # Test non-positive CPUE
  expect_error(
    ProductionModel(
      years = years, catch = catch,
      cpue = c(0, 2.0, 1.9, 2.2, 2.1, 1.8), effort = effort
    ),
    "cpue data must be positive finite numbers"
  )

  # Test non-positive effort
  expect_error(
    ProductionModel(
      years = years, catch = catch, cpue = cpue,
      effort = c(-100, effort[-1])
    ),
    "Element 1 is not >= 0"
  )

  # Test missing values
  expect_error(
    ProductionModel(years = years, catch = c(NA, catch[-1]), cpue = cpue, effort = effort),
    "Contains missing values"
  )

  expect_error(
    ProductionModel(years = years, catch = catch, cpue = c(Inf, cpue[-1]), effort = effort),
    "Must be finite"
  )
})

test_that("ProductionModel accepts custom parameters", {
  years <- 2000:2005
  catch <- c(1000, 1100, 1200, 950, 1050, 1300)
  cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
  effort <- catch / cpue

  custom_params <- c(
    r = 0.15,
    K = 60000,
    m = 1.5,
    q = 0.0015,
    sigma_proc = 0.12,
    sigma_obs = 0.18
  )

  model <- ProductionModel(
    years = years, catch = catch, cpue = cpue,
    effort = effort, parameters = custom_params
  )

  expect_equal(parameters(model), custom_params)
  expect_true(validObject(model))
})

test_that("ProductionModel validates parameter constraints", {
  years <- 2000:2005
  catch <- c(1000, 1100, 1200, 950, 1050, 1300)
  cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
  effort <- catch / cpue

  # Test missing parameter names
  expect_error(
    ProductionModel(
      years = years, catch = catch, cpue = cpue, effort = effort,
      parameters = c(0.1, 50000, 2.0, 0.001, 0.1, 0.2)
    ),
    "Must have names"
  )

  # Test missing required parameters
  expect_error(
    ProductionModel(
      years = years, catch = catch, cpue = cpue, effort = effort,
      parameters = c(r = 0.1, K = 50000, m = 2.0)
    ),
    "must include the elements"
  )

  # Test negative parameters
  expect_error(
    ProductionModel(
      years = years, catch = catch, cpue = cpue, effort = effort,
      parameters = c(
        r = -0.1, K = 50000, m = 2.0, q = 0.001,
        sigma_proc = 0.1, sigma_obs = 0.2
      )
    ),
    "Parameter 'r' must be positive"
  )

  expect_error(
    ProductionModel(
      years = years, catch = catch, cpue = cpue, effort = effort,
      parameters = c(
        r = 0.1, K = -50000, m = 2.0, q = 0.001,
        sigma_proc = 0.1, sigma_obs = 0.2
      )
    ),
    "Parameter 'K' must be positive"
  )
})

test_that("ProductionModel accessor methods work correctly", {
  years <- 2000:2005
  catch <- c(1000, 1100, 1200, 950, 1050, 1300)
  cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
  effort <- catch / cpue

  model <- ProductionModel(years = years, catch = catch, cpue = cpue, effort = effort)

  # Test data accessor
  model_data <- model_data(model)
  expect_type(model_data, "list")
  expect_equal(model_data$years, years)
  expect_equal(model_data$catch, catch)
  expect_equal(model_data$cpue, cpue)
  expect_equal(model_data$effort, effort)

  # Test parameters accessor
  model_params <- parameters(model)
  expect_type(model_params, "double")
  expect_true(all(c("r", "K", "m", "q", "sigma_proc", "sigma_obs") %in% names(model_params)))
  expect_true(all(model_params > 0))

  # Test results accessor (should be empty for unfitted model)
  model_results <- results(model)
  expect_type(model_results, "list")
  expect_equal(length(model_results), 0)

  # Test fitted accessor
  expect_false(fitted(model))
})

test_that("ProductionModel print method works", {
  years <- 2000:2005
  catch <- c(1000, 1100, 1200, 950, 1050, 1300)
  cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
  effort <- catch / cpue

  model <- ProductionModel(years = years, catch = catch, cpue = cpue, effort = effort)

  # Test that print doesn't error and returns invisibly
  output <- capture.output(result <- print(model))
  expect_identical(result, model)
  expect_true(length(output) > 0)
  expect_true(any(grepl("Production Model Object", output)))
  expect_true(any(grepl("pella_tomlinson", output)))
  expect_true(any(grepl("Fitted: No", output)))
})

test_that("ProductionModel show method works", {
  years <- 2000:2005
  catch <- c(1000, 1100, 1200, 950, 1050, 1300)
  cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
  effort <- catch / cpue

  model <- ProductionModel(years = years, catch = catch, cpue = cpue, effort = effort)

  # Test that show works (it calls print)
  output <- capture.output(result <- show(model))
  expect_identical(result, model)
  expect_true(length(output) > 0)
})

test_that("ProductionModel validity function catches invalid objects", {
  # Create a valid model first
  years <- 2000:2005
  catch <- c(1000, 1100, 1200, 950, 1050, 1300)
  cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
  effort <- catch / cpue

  model <- ProductionModel(years = years, catch = catch, cpue = cpue, effort = effort)

  # Test that manually creating invalid objects fails validation

  # Invalid model_type
  invalid_model <- model
  invalid_model@model_type <- "invalid_type"
  expect_error(validObject(invalid_model), "model_type must be 'pella_tomlinson'")

  # Invalid fitted status
  invalid_model <- model
  invalid_model@fitted <- c(TRUE, FALSE)
  expect_error(validObject(invalid_model), "fitted must be a single logical value")

  # Invalid creation_date (test by creating new object with invalid slot)
  expect_error(
    {
      new("ProductionModel",
        parameters = model@parameters,
        data = model@data,
        results = model@results,
        fitted = model@fitted,
        model_type = model@model_type,
        creation_date = c(Sys.time(), Sys.time())
      ) # Multiple dates
    },
    "creation_date must be a single POSIXct value"
  )
})

test_that("ProductionModel handles edge cases", {
  # Single data point
  model_single <- ProductionModel(
    years = 2000L,
    catch = 1000,
    cpue = 2.0,
    effort = 500
  )
  expect_true(validObject(model_single))
  expect_equal(length(model_data(model_single)$years), 1)

  # Large dataset
  years_large <- 1950:2023
  n <- length(years_large)
  catch_large <- rnorm(n, 1000, 100)
  catch_large <- pmax(catch_large, 0) # Ensure non-negative
  cpue_large <- exp(rnorm(n, log(2), 0.3))
  effort_large <- catch_large / cpue_large

  model_large <- ProductionModel(
    years = years_large,
    catch = catch_large,
    cpue = cpue_large,
    effort = effort_large
  )
  expect_true(validObject(model_large))
  expect_equal(length(model_data(model_large)$years), n)
})
