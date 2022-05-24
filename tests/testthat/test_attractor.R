context("Attractor")

test_that("initializes with raster layer", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  attr_rast <- template*2
  expect_error(attractor <- Attractor(attr_rast, 1:10),
               "Region model must be a 'Region' or inherited class object.")
  expect_error(attractor <- Attractor(
    terra::crop(attr_rast, terra::ext(attr_rast)*0.9), region),
    "The spatial object x should be compatible with that defining the region.")
  expect_silent(attractor <- Attractor(attr_rast, region))
  expect_is(attractor, "Attractor")
  expect_equal(attractor$get_type(), "source")
})

test_that("gets values from raster attractor layer", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  region$set_aggr(aggr_factor = 5, inner_radius = 10000)
  attr_rast <- template*2
  attractor <- Attractor(attr_rast, region, multiplier = 0.3)
  expect_is(attractor$get_rast(), "SpatRaster")
  expect_silent(attr_values <- attractor$get_rast()[region$get_indices()][,1])
  expect_equal(attr_values, attr_rast[region$get_indices()][,1]*0.3)
  idx <- which(attr_values > 0.3)
  expect_equal(attractor$get_values(), attr_values)
  expect_equal(attractor$get_values(idx), attr_values[idx])
  expect_equal(attractor$get_values(idx, na.incl = TRUE),
               attr_rast[idx][,1]*0.3)
  aggr <- region$get_aggr()
  attr_aggr_rast <- terra::aggregate(attr_rast, fact = aggr$factor,
                                     fun = "mean", na.rm = TRUE)
  idx <- which(attr_aggr_rast[aggr$indices][,1] > 0.9)
  expect_equal(attractor$get_aggr_values(),
               attr_aggr_rast[aggr$indices][,1]*0.3)
  expect_equal(attractor$get_aggr_values(idx),
               attr_aggr_rast[aggr$indices[idx]][,1]*0.3)
  expect_equal(attractor$get_aggr_values(idx, na.incl = TRUE),
               attr_aggr_rast[idx][,1]*0.3)
})

test_that("initializes with numeric vector", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  expect_error(attractor <- Attractor(14:1, 1:10),
               "Region model must be a 'Region' or inherited class object.")
  expect_error(attractor <- Attractor(1:10, region),
               paste("The vector x should have values for each location",
                     "defined in the region."))
  expect_silent(attractor <- Attractor(14:1, region, type = "both",
                                       multiplier = 0.3))
  expect_is(attractor, "Attractor")
  expect_equal(attractor$get_type(), "both")
})

test_that("gets values from numeric attractor vector with grid region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  attractor <- Attractor(1:region$get_locations(), region, multiplier = 0.3)
  expect_is(attractor$get_rast(), "SpatRaster")
  expect_silent(attr_values <- attractor$get_rast()[region$get_indices()][,1])
  expect_equal(attr_values, (1:region$get_locations())*0.3)
  idx <- which(attr_values > 3500)
  expect_equal(attractor$get_values(), attr_values)
  expect_equal(attractor$get_values(idx), attr_values[idx])
})

test_that("gets values from numeric attractor vector with patch region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  attractor <- Attractor(14:1, region, type = "both", multiplier = 0.3)
  expect_null(attractor$get_rast())
  expect_equal(attractor$get_values(), (14:1)*0.3)
  expect_equal(attractor$get_values(1:5), (14:10)*0.3)
})
