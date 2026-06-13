# Tests for package-level reference point defaults

test_that("set/get reference-point defaults work", {
  old_target <- getOption("SurplusProductionModel.biomass_target_default")
  old_baseline <- getOption("SurplusProductionModel.baseline_default")
  on.exit({
    do.call(options, setNames(list(old_target), "SurplusProductionModel.biomass_target_default"))
    do.call(options, setNames(list(old_baseline), "SurplusProductionModel.baseline_default"))
  }, add = TRUE)

  set_reference_point_defaults(biomass_target = 0.4, baseline = "K")
  defaults <- get_reference_point_defaults()

  expect_equal(defaults$biomass_target, 0.4)
  expect_equal(defaults$baseline, "K")

  # Explicit NULL should clear biomass target default
  set_reference_point_defaults(biomass_target = NULL)
  defaults2 <- get_reference_point_defaults()
  expect_null(defaults2$biomass_target)
})

test_that("calculate_reference_points uses package defaults when args are omitted", {
  old_target <- getOption("SurplusProductionModel.biomass_target_default")
  old_baseline <- getOption("SurplusProductionModel.baseline_default")
  on.exit({
    do.call(options, setNames(list(old_target), "SurplusProductionModel.biomass_target_default"))
    do.call(options, setNames(list(old_baseline), "SurplusProductionModel.baseline_default"))
  }, add = TRUE)

  set_reference_point_defaults(biomass_target = 0.4, baseline = "K")

  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6),
    parameters = c(r = 0.3, K = 5000, m = 2.0, d0 = 1.2, q = 0.001, sigma_proc = 0.2, sigma_obs = 0.3)
  )
  fitted_model$fitted <- TRUE

  ref_points <- calculate_reference_points(fitted_model)

  expect_equal(ref_points$target_reference_points$baseline, "K")
  expect_equal(ref_points$target_reference_points$biomass_name, "B_40%K")
  expect_equal(ref_points$target_reference_points$f_name, "F_40%K")
})

test_that("explicit arguments override package defaults", {
  old_target <- getOption("SurplusProductionModel.biomass_target_default")
  old_baseline <- getOption("SurplusProductionModel.baseline_default")
  on.exit({
    do.call(options, setNames(list(old_target), "SurplusProductionModel.biomass_target_default"))
    do.call(options, setNames(list(old_baseline), "SurplusProductionModel.baseline_default"))
  }, add = TRUE)

  set_reference_point_defaults(biomass_target = 0.4, baseline = "K")

  fitted_model <- ProductionModel(
    years = 2000:2005,
    catch = rep(1000, 6),
    cpue = rep(1.5, 6),
    effort = rep(1000 / 1.5, 6),
    parameters = c(r = 0.3, K = 5000, m = 2.0, d0 = 1.2, q = 0.001, sigma_proc = 0.2, sigma_obs = 0.3)
  )
  fitted_model$fitted <- TRUE

  ref_points <- calculate_reference_points(
    fitted_model,
    biomass_target = 0.35,
    baseline = "B_initial"
  )

  expect_equal(ref_points$target_reference_points$baseline, "B_initial")
  expect_equal(ref_points$target_reference_points$biomass_name, "B_35%B_initial")
  expect_equal(ref_points$target_reference_points$f_name, "F_35%B_initial")

  # Explicit NULL should suppress default target output
  ref_points_no_target <- calculate_reference_points(fitted_model, biomass_target = NULL)
  expect_null(ref_points_no_target$target_reference_points)
})
