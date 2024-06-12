context("Incursions")

test_that("initializes with raster layer", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template*0)
  incursion_rast <- template*10
  expect_silent(incursions <- Incursions(incursion_rast, region,
                                         incursion_mean = 5,
                                         incursion_stages = 2:3))
  expect_is(incursions, "Incursions")
  expect_equal(incursions$get_type(), "weight")
  expect_equal(incursions$get_incursion_mean(), 5)
  expect_equal(incursions$get_incursion_stages(), 2:3)
  expect_error(incursions <- Incursions(incursion_rast, 1:10),
               "Region model must be a 'Region' or inherited class object.")
  expect_error(incursions <- Incursions(1:10, region),
               paste("Vector x length must be equal to the number of region",
                     "locations."))
})

test_that("generates weighted incursion location", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template*0)
  incursion_rast <- template*10
  incursions <- Incursions(incursion_rast, region, type = "weight")
  expect_silent(generated_incursion <- incursions$generate())
  expect_is(generated_incursion, "logical")
  expect_equal(length(which(generated_incursion)), 1)
  expect_equal(attr(generated_incursion, "type"), "weight")
})

test_that("generates probability incursions", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template*0)
  incursion_rast <- template
  incursions <- Incursions(incursion_rast, region, type = "prob")
  mean_incursions <- (mean(template[region$get_indices()][,1])*
                        region$get_locations())
  expect_silent(generated_incursions <- incursions$generate())
  expect_is(generated_incursions, "logical")
  expect_true(abs(length(which(generated_incursions)) -
                    mean_incursions)/mean_incursions < 0.1)
  expect_equal(attr(generated_incursions, "type"), "prob")
})

