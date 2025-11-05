context("UnstructPopulation")

test_that("initializes with region and other parameters", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  expect_error(population <- UnstructPopulation(region, growth = c(1.1, 1.2)),
               paste("Population growth should be a single value (e.g. 1.2)",
                     "or a matrix with a single row or a row for each region",
                     "location."), fixed = TRUE)
  expect_error(population <- UnstructPopulation(region, growth = -1),
               "Population growth values should be >= 0.")
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
  set.seed(1243)
  expect_silent(n <- population$make(incursion = incursion))
  idx <- which(n > 0)
  set.seed(1243)
  expect_silent(n1 <- population$grow(n, 1))
  expect_equal(round(mean(n1[idx]/n[idx]), 1), 1.2)
  # growth variation
  growth <- rep(1.2, region$get_locations())
  growth <- cbind(growth, growth*0.8, growth*0.6)
  expect_silent(population <- UnstructPopulation(region, growth = growth,
                                                 incursion_mean = 10))
  set.seed(1243)
  expect_silent(n1 <- population$grow(n, 1))
  expect_equal(round(sum(n1[idx])/sum(n[idx]), 3), 1.2)
  set.seed(1243)
  expect_silent(n2 <- population$grow(n, 2))
  expect_equal(round(sum(n2[idx])/sum(n[idx]), 2), 1.2*0.8)
  set.seed(1243)
  expect_silent(n3 <- population$grow(n, 3))
  expect_equal(round(sum(n3[idx])/sum(n[idx]), 2), 1.2*0.6)
  set.seed(1243)
  expect_equal(population$grow(n, 5), n2)
  # growth control
  expect_silent(population <- UnstructPopulation(region, growth = 1.2,
                                                 incursion_mean = 10))
  set.seed(1243)
  expect_silent(n2_no_control <- population$grow(n, 2))
  attr(n, "control_growth") <- c(rep(0, 4000), rep(0.5, 4000), rep(1, region$get_locations() - 8000))
  set.seed(1243)
  expect_silent(n2_control_growth <- population$grow(n, 2))
  expect_equal(n2_control_growth[1:4000], rep(0, 4000))
  expect_true(
    abs(sum(n2_control_growth[4001:8000])/
          (sum(n2_no_control[4001:8000])*0.5) - 1) < 0.01)
  expect_true(
    abs(sum(n2_control_growth[8001:region$get_locations()])/
          sum(n2_no_control[8001:region$get_locations()]) - 1) < 0.01)
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
  # dynamic capacity multiplier
  mult <- (capacity > 10)*0.5
  expect_silent(population$set_capacity_mult(mult))
  expect_equal(population$get_capacity(), capacity*mult)
  expect_silent(population$set_capacity_mult(2))
  expect_equal(population$get_capacity(), capacity*2)
  attr(n, "dynamic_mult") <- list(NULL, list(mult))
  attr(attr(n, "dynamic_mult")[[2]], "links") <- c("capacity", "suitability")
  expect_silent(population$set_capacity_mult(n))
  expect_equal(population$get_capacity(), capacity*mult)
  # growth variation
  growth <- rep(1.2, region$get_locations())
  growth <- cbind(growth, growth*0.8, growth*0.6)
  expect_silent(population <- UnstructPopulation(region, growth = growth,
                                                 capacity = capacity,
                                                 incursion_mean = 10))
  set.seed(1243)
  expect_equal(mean(population$grow(n, 1)[idx]/n[idx]), mean(n1[idx]/n[idx]))
  expect_true(abs(round(mean(population$grow(n, 2)[idx]/n[idx]), 2) -
                    0.96) <= 0.02)
  expect_true(abs(round(mean(population$grow(n, 3)[idx]/n[idx]), 2) -
                    0.72) <= 0.02)
  expect_true(abs(round(mean(population$grow(n, 5)[idx]/n[idx]), 2) -
                    0.96) <= 0.02)
  # temporal capacity
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
