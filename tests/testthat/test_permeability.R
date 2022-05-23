context("Permeability")

test_that("initializes with raster layer", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  perm_rast <- region$get_template()
  perm_rast[region$get_indices()] <- round(region$get_indices()/1560)/10
  expect_silent(perm <- Permeability(perm_rast, region))
  expect_is(perm, "Permeability")
  expect_is(perm$get_rast(), "SpatRaster")
  perm_data <- matrix(c( 1,  2, 0.8,
                         1,  3, 0.7,
                         1,  4, 0.7,
                         10, 11, 0.8), ncol = 3, byrow = TRUE)
  expect_error(perm <- Permeability(perm_rast, perm_data),
               "Region model must be a 'Region' or inherited class object.")
  expect_error(perm <- Permeability(perm_data, region),
               paste("The spatial object x should be compatible with that",
                     "defining the region."))
})

test_that("initializes with patch/network data", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  perm_data <- matrix(c( 1,  2, 0.8,
                         1,  3, 0.7,
                         1,  4, 0.7,
                         1,  6, 0.6,
                         1, 10, 0.8,
                         2,  5, 0.6,
                         3, 12, 0.6,
                         3, 14, 0.5,
                         4, 13, 0.5,
                         5,  8, 0.6,
                         6,  9, 0.6,
                         7, 10, 0.5,
                         10, 11, 0.8), ncol = 3, byrow = TRUE)
  colnames(perm_data) <- c("i", "j", "weight")
  expect_silent(perm <- Permeability(perm_data, region))
  expect_is(perm$get_data(), "data.frame")
  perm_rast <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  expect_error(perm <- Permeability(perm_data, perm_rast),
               "Region model must be a 'Region' or inherited class object.")
  expect_error(perm <- Permeability(perm_rast, region),
               paste("The spatial object x should be compatible with that",
                     "defining the region."))
})
