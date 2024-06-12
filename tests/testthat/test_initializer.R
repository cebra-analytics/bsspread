context("Initializer")

test_that("initializes with raster layer", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  initial_rast <- round(template*10)
  expect_error(initializer <- Initializer(1:10, region),
               paste("Vector (or array) x length (or number of rows) must be",
                     "equal to the number of region locations."), fixed = TRUE)
  expect_error(initializer <- Initializer(initial_rast, 1:10),
               "Region model must be a 'Region' or inherited class object.")
  expect_error(initializer <- Initializer(initial_rast, region,
                                          population_model = 1:10),
               "Population model must be a 'Population' or inherited class object.")
  pop_model <- Population(region)
  expect_silent(initializer <- Initializer(initial_rast, region,
                                           population_model = pop_model))
  expect_is(initializer, "Initializer")
  expect_silent(n <- initializer$initialize())
  expect_equal(n, initial_rast[region$get_indices()][,1])
})

test_that("initializes with incursion object", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  incursion_rast <- template/10
  incursions <- Incursions(incursion_rast, region, type = "prob",
                           continued = TRUE, incursion_mean = 5,
                           incursion_stages = 2:3)
  pop_model <- Population(region)
  expect_silent(initializer <- Initializer(incursions, region,
                                           population_model = pop_model))
  expect_silent(n <- initializer$initialize())
  idx <- which(n > 0)
  expect_true(length(idx) > 1)
  expect_equal(round(mean(n[idx])), 5)
  expect_is(initializer$continued_incursions, "function")
  expect_silent(continued_incursions <- initializer$continued_incursions())
  expect_is(continued_incursions, "function")
  expect_silent(n <- continued_incursions(tm = 2, n))
  expect_true(length(which(n > 0)) > length(idx))
  expect_true(all(idx %in% which(n > 0)))
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  pop_model <- StagedPopulation(region, growth = stage_matrix)
  expect_silent(initializer <- Initializer(incursions, region,
                                           population_model = pop_model))
  expect_silent(n <- initializer$initialize())
  expect_true(all(n[,1] == 0))
  idx <- which(rowSums(n) > 0)
  expect_equal(round(mean(rowSums(n[idx,]))), 5)
})
