context("UnstructPopulation")

test_that("initializes with region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  expect_error(population <- UnstructPopulation(region, growth = c(1.1, 1.2)),
               "Population growth should be a single value (e.g. 1.2).",
               fixed = TRUE)
  expect_error(population <- UnstructPopulation(region, growth = 1.2,
                                                capacity = 30),
               paste("Population capacity should be a vector with a value for",
                     "each region location."))
  expect_silent(population <- UnstructPopulation(region, growth = 1.2))
  expect_is(population, "UnstructPopulation")
  expect_s3_class(population, "Population")
  expect_equal(population$get_type(), "unstructured")
})

test_that("grows populations without capacity", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  expect_silent(population <- UnstructPopulation(region, growth = 1.2,
                                                 incursion_mean = 10))
  idx <- which(template[region$get_indices()] > 0)
  incursion <- rep(FALSE, region$get_locations())
  incursion[idx] <- TRUE
  expect_silent(n <- population$make(incursion = incursion))
  idx <- which(n > 0)
  expect_equal(round(mean(population$grow(n)[idx]/n[idx]), 1), 1.2)
})

test_that("grows populations with capacity", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  capacity <- template[region$get_indices()][,1]*30
  expect_silent(population <- UnstructPopulation(region, growth = 1.2,
                                                 capacity = capacity,
                                                 incursion_mean = 10))
  idx <- which(template[region$get_indices()] > 0)
  incursion <- rep(FALSE, region$get_locations())
  incursion[idx] <- TRUE
  expect_silent(n <- population$make(incursion = incursion))
  idx <- which(n > 0)
  expect_equal(round(mean(population$grow(n)[idx]/n[idx]), 1), 1.0)
})
