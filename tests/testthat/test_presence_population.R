context("PresencePopulation")

test_that("initializes with region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  expect_silent(population <- PresencePopulation(region))
  expect_is(population, "PresencePopulation")
  expect_s3_class(population, "Population")
  expect_equal(population$get_type(), "presence_only")
})

test_that("makes populations with spread delay", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  expect_silent(population <- PresencePopulation(region, spread_delay = 2))
  expect_silent(initial_pop <- population$make(initial = TRUE))
  expect_named(attributes(initial_pop), c("type", "spread_delay"))
  expect_length(attr(initial_pop, "spread_delay"), region$get_locations())
  expect_true(all(is.na(attr(initial_pop, "spread_delay"))))
})

test_that("grows populations with spread delay", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- PresencePopulation(region, spread_delay = 2)
  idx <- which(template[region$get_indices()] > 0.5)
  incursion <- rep(FALSE, region$get_locations())
  incursion[idx] <- TRUE
  expect_silent(n <- population$make(incursion = incursion))
  expect_true(all(is.na(attr(n, "spread_delay"))))
  expect_silent(n <- population$grow(n))
  expect_true(all(attr(n, "spread_delay")[idx] == 2))
  n <- population$grow(n)
  expect_true(all(attr(n, "spread_delay")[idx] == 1))
  n <- population$grow(n)
  expect_true(all(attr(n, "spread_delay")[idx] == 0))
  n <- population$grow(n)
  expect_true(all(attr(n, "spread_delay")[idx] == 0)) # again
  expect_true(all(is.na(attr(n, "spread_delay")[-idx])))
})
