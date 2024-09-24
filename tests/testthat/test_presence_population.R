context("PresencePopulation")

test_that("initializes with region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  expect_silent(population <- PresencePopulation(region))
  expect_is(population, "PresencePopulation")
  expect_s3_class(population, "Population")
  expect_is(population$get_region(), "Region")
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
  initial_n <- template[region$get_indices()][,1] > 0.5
  idx <- which(initial_n)
  initial_age <- initial_n*0
  initial_age[idx[1:150]] <- c(rep(3,50), rep(2,50), rep(1,50))
  attr(initial_n, "age") <- initial_age
  expect_silent(initial_pop <- population$make(initial = initial_n))
  expect_true(all(is.na(attr(initial_pop, "spread_delay")[-idx])))
  expect_equal(attr(initial_pop, "spread_delay")[idx],
               c(rep(0,50), rep(1,50), rep(2,50), rep(3,155)))
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
  initial_age <- incursion*0
  initial_age[idx[1:150]] <- c(rep(3,50), rep(2,50), rep(1,50))
  attr(incursion, "age") <- initial_age
  expect_silent(n <- population$make(initial = incursion))
  expect_equal(attr(n, "spread_delay")[idx],
               c(rep(0,50), rep(1,50), rep(2,50), rep(3,155)))
  expect_silent(n <- population$grow(n))
  expect_equal(attr(n, "spread_delay")[idx],
               c(rep(0,100), rep(1,50), rep(2,155)))
  n <- population$grow(n) # silent
  expect_equal(attr(n, "spread_delay")[idx], c(rep(0,150), rep(1,155)))
  n <- population$grow(n)
  expect_true(all(attr(n, "spread_delay")[idx] == 0))
  n <- population$grow(n)
  expect_true(all(attr(n, "spread_delay")[idx] == 0)) # again
  expect_true(all(is.na(attr(n, "spread_delay")[-idx])))
})
