context("Population")

test_that("initializes with region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  expect_silent(population <- Population(region))
  expect_is(population, "Population")
  expect_is(population$get_region(), "Region")
  expect_equal(population$get_type(), "presence_only")
})

test_that("makes populations with initial values", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  expect_equal(population$make(initial = 10), rep(10, region$get_locations()))
  expect_error(population$make(initial = 0:20),
               paste("Initial population values must be consistent with the",
                     "number of region locations."))
  initial_pop <- round(template[region$get_indices()]*30)[,1]
  expect_equal(population$make(initial = initial_pop), initial_pop)
  attr(initial_pop, "age") <- 2
  expect_equal(population$make(initial = initial_pop), initial_pop)
})

test_that("makes populations with incursions", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  incursion <- template[region$get_indices()][,1] > 0
  expect_equal(population$make(incursion = incursion), incursion)
  expect_error(population <- Population(region, incursion_mean = 0),
               paste("Incursion mean population size should be a numeric",
                     "value greater than zero."))
  expect_silent(population <- Population(region, incursion_mean = 10))
  expect_silent(n <- population$make(incursion = incursion))
  expect_is(n, "integer")
  expect_true(sum(n > 0) <= sum(incursion))
  expect_equal(round(mean(n[which(incursion)])), 10)
  expect_silent(population <- Population(region))
  expect_silent(population$set_incursion_mean(15))
  expect_silent(n <- population$make(incursion = incursion))
  expect_equal(round(mean(n[which(incursion)])), 15)
})

test_that("makes populations with establishment prob", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  establish <- template[region$get_indices()][,1]
  expect_error(Population(region, establish_pr = matrix(1:12, ncol = 3)/20),
               paste("Establishment probability should be a vector or matrix",
                     "with a value or row for each region location."))

  expect_silent(population <- Population(region, establish_pr = establish))
  incursion <- rep(TRUE, region$get_locations())
  attr(incursion, "type") <- "weight"
  expect_is(population$make(incursion = incursion), "logical")
  expect_equal(sum(population$make(incursion = incursion)),
               region$get_locations())
  attr(incursion, "type") <- "prob"
  set.seed(1234)
  expect_silent(expected_t1 <- population$make(incursion = incursion))
  sum(expected_t1) < region$get_locations() # TRUE
  sum(expected_t1) <= sum(establish > 0) # TRUE
  set.seed(1234)
  all(population$make(incursion = incursion, tm = 2) == expected_t1)
  temp_establish <- cbind(establish, establish*0.8, establish*0.6)
  expect_silent(population <- Population(region,
                                         establish_pr = temp_establish))
  set.seed(1234)
  expect_equal(population$make(incursion = incursion, tm = 0), expected_t1)
  set.seed(1234)
  expect_silent(expected_t2 <- population$make(incursion = incursion, tm = 2))
  expect_true(sum(expected_t2) < sum(expected_t1))
  set.seed(1234)
  expect_true(sum(population$make(incursion = incursion, tm = 3)) <
                sum(expected_t2))
  set.seed(1234)
  expect_equal(population$make(incursion = incursion, tm = 4), expected_t1)
})
