context("UnstructPopulation")

test_that("initializes with region and other parameters", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  expect_error(population <- UnstructPopulation(region, growth = c(1.1, 1.2)),
               "Population growth should be a single value (e.g. 1.2).",
               fixed = TRUE)
  expect_error(population <- UnstructPopulation(region, growth = 1.2,
                                                capacity = 30),
               paste("Population capacity should be a vector or matrix",
                     "with a value or row for each region location."))
  expect_silent(population <- UnstructPopulation(region, growth = 1.2))
  expect_is(population, "UnstructPopulation")
  expect_s3_class(population, "Population")
  expect_is(population$get_region(), "Region")
  expect_equal(population$get_type(), "unstructured")
  region <- Region()
  expect_error(population <- UnstructPopulation(region, growth = 1.2,
                                                capacity = 30),
               paste("Population capacity area is required when capacity is",
                     "specified and the region is spatially implicit"))
  expect_error(population <- UnstructPopulation(region, growth = 1.2,
                                                capacity = 30,
                                                capacity_area = 0),
               paste("Population capacity area should be a numeric value",
                     "> 0."))
  expect_silent(population <- UnstructPopulation(region, growth = 1.2,
                                                 capacity = 30,
                                                 capacity_area = 1e+06))
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
  expect_equal(round(mean(population$grow(n, 1)[idx]/n[idx]), 1), 1.2)
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
  set.seed(1243)
  expect_silent(n <- population$make(incursion = incursion))
  idx <- which(n > 0)
  set.seed(1243)
  expect_silent(n1 <- population$grow(n, 1))
  expect_equal(round(mean(n1[idx]/n[idx]), 1), 1.0)
  temp_capacity <- cbind(capacity, capacity*0.8, capacity*1.2)
  expect_silent(population <- UnstructPopulation(region, growth = 1.2,
                                                 capacity = temp_capacity,
                                                 incursion_mean = 10))
  set.seed(1243)
  expect_silent(n1_temp <- population$grow(n, 1))
  expect_equal(n1_temp[idx], n1[idx])
  mean_growth_1 <- mean(n1_temp[idx]/n[idx])
  set.seed(1243)
  expect_silent(mean_growth_2 <- mean(population$grow(n, 2)[idx]/n[idx]))
  expect_true(mean_growth_2 < mean_growth_1)
  set.seed(1243)
  expect_silent(mean_growth_3 <- mean(population$grow(n, 3)[idx]/n[idx]))
  expect_true(mean_growth_3 > mean_growth_1)
  set.seed(1243)
  expect_silent(n4_temp <- population$grow(n, 4))
  expect_equal(n4_temp[idx], n1[idx])
  region <- Region()
  population <- UnstructPopulation(region, growth = 1.2, capacity = 300,
                                   capacity_area = 1e+06)
  expect_equal(attr(population$get_capacity(), "area"), 1e+06)
  n <- 100
  set.seed(1243); new_n <- stats::rpois(1, 1.2*100)
  set.seed(1243)
  expect_equal(population$grow(n, 1), new_n)
  region$set_max_implicit_area(1e+06)
  r <- exp(log(1.2)*(1 - n/300))
  set.seed(1243); new_n <- stats::rpois(1, r*100)
  set.seed(1243)
  expect_equal(population$grow(n, 1), new_n)
  attr(n, "diffusion_rate") <- 2000
  attr(n, "diffusion_radius") <- 1000
  r <- exp(log(1.2)*(1 - n/(300*pi*3000^2/1e+06)))
  set.seed(1243); new_n <- stats::rpois(1, r*100)
  attr(new_n, "diffusion_rate") <- 2000
  attr(new_n, "diffusion_radius") <- 1000
  set.seed(1243)
  expect_equal(population$grow(n, 1), new_n)
})
