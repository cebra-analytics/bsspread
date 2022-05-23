context("StagedPopulation")

test_that("initializes with region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)

  expect_error(population <- StagedPopulation(region, growth = c(1.1, 1.2)),
               "Population growth should be a square matrix (at least 2 x 2).",
               fixed = TRUE)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  attr(stage_matrix, "survivals") <- 0.5
  expect_error(population <- StagedPopulation(region, growth = stage_matrix),
               paste("Growth survivals attribute should be a numeric matrix",
                     "with the same dimensions as the growth (stage/age)",
                     "matrix."), fixed = TRUE)
  attr(stage_matrix, "survivals") <-
    matrix(c(0.0, 0.0, 0.0,
             0.3, 0.0, 0.0,
             0.0, 0.6, 0.8),
           nrow = 3, ncol = 3, byrow = TRUE)
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix))
  expect_is(population, "StagedPopulation")
  expect_s3_class(population, "Population")
  expect_equal(population$get_type(), "stage_structured")
})

test_that("makes populations with incursions", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  attr(stage_matrix, "survivals") <-
    matrix(c(0.0, 0.0, 0.0,
             0.3, 0.0, 0.0,
             0.0, 0.6, 0.8),
           nrow = 3, ncol = 3, byrow = TRUE)
  expect_error(population <- StagedPopulation(region, growth = stage_matrix,
                                 incursion_stages = 2:4),
               "Incursion stages should specify index values between 1 and 3.")
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                               incursion_stages = 2:3,
                                               incursion_mean = 10))
  incursion <- template[region$get_indices()][,1] > 0
  expect_silent(n <- population$make(incursion = incursion))
  expect_equal(dim(n), c(region$get_locations(), 3))
  expect_true(all(n[,1] == 0))
  pop_tot <- rowSums(n)
  expect_equal(round(mean(pop_tot[which(pop_tot > 0)])), 10)
})

test_that("grows populations without capacity", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  attr(stage_matrix, "survivals") <-
    matrix(c(0.0, 0.0, 0.0,
             0.3, 0.0, 0.0,
             0.0, 0.6, 0.8),
           nrow = 3, ncol = 3, byrow = TRUE)
  population <- StagedPopulation(region, growth = stage_matrix)
  idx <- which(template[region$get_indices()] > 0)
  initial <- rep(0, region$get_locations())
  initial[idx] <- stats::rpois(length(idx), 20)
  expect_silent(n <- population$make(initial = initial))
  idx <- which(rowSums(n) > 0)
  expect_true(abs(mean(rowSums(population$grow(n)[idx,]))/20 -
                    Re((eigen(stage_matrix,)$values)[1])) < 0.05)
})

test_that("grows populations with capacity", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  attr(stage_matrix, "survivals") <-
    matrix(c(0.0, 0.0, 0.0,
             0.3, 0.0, 0.0,
             0.0, 0.6, 0.8),
           nrow = 3, ncol = 3, byrow = TRUE)
  capacity <- template[region$get_indices()][,1]*10
  expect_error(population <- StagedPopulation(region, growth = stage_matrix,
                                              capacity = capacity,
                                              capacity_stages = 2:4),
               "Capacity stages should specify index values between 1 and 3.")
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                               capacity = capacity,
                                               capacity_stages = 2:3))

  idx <- which(template[region$get_indices()] > 0)
  initial <- rep(0, region$get_locations())
  initial[idx] <- stats::rpois(length(idx), capacity[idx]*3)
  expect_silent(n <- population$make(initial = initial))
  idx <- which(rowSums(n) > 0)
  expect_equal(round(mean(rowSums(population$grow(n)[idx,]))/
                       mean(rowSums(n[idx,])), 1), 1.0)
})
