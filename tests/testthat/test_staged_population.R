context("StagedPopulation")

test_that("initializes with region and other parameters", {
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
  expect_equal(attr(population$get_growth(), "labels"),
               c("stage 1", "stage 2", "stage 3"))
  attr(stage_matrix, "labels") <- c("a", "b")
  expect_error(population <- StagedPopulation(region, growth = stage_matrix),
               paste("Stage/age labels attribute should be a character vector",
                     "compatible with the dimensions of the growth",
                     "(stage/age) matrix."), fixed = TRUE)
  attr(stage_matrix, "labels") <- c("a", "b", "c")
  expect_error(population <- StagedPopulation(region, growth = stage_matrix,
                                                capacity = 30),
               paste("Population capacity should be a vector with a value for",
                     "each region location."))

  expect_silent(population <- StagedPopulation(region, growth = stage_matrix))
  expect_equal(attr(population$get_growth(), "labels"), c("a", "b", "c"))
  expect_equal(population$get_growth_r(),
               Re((eigen(stage_matrix, only.values = TRUE)$values)[1]))
  expect_is(population, "StagedPopulation")
  expect_s3_class(population, "Population")
  expect_is(population$get_region(), "Region")
  expect_equal(population$get_type(), "stage_structured")
  region <- Region()
  expect_error(population <- StagedPopulation(region, growth = stage_matrix,
                                              capacity = 30),
               paste("Population capacity area is required when capacity is",
                     "specified and the region is spatially implicit"))
  expect_error(population <- StagedPopulation(region, growth = stage_matrix,
                                              capacity = 30,
                                              capacity_area = 0),
               paste("Population capacity area should be a numeric value",
                     "> 0."))
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                               capacity = 30,
                                               capacity_area = 1e+06))
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
  attr(stage_matrix, "labels") <- c("a", "b", "c")
  expect_error(population <- StagedPopulation(region, growth = stage_matrix,
                                 incursion_stages = 2:4),
               "Incursion stages should specify index values between 1 and 3.")
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                               incursion_stages = 2:3,
                                               incursion_mean = 10))
  incursion <- template[region$get_indices()][,1] > 0
  expect_silent(n <- population$make(incursion = incursion))
  expect_equal(dim(n), c(region$get_locations(), 3))
  expect_equal(colnames(n), c("a", "b", "c"))
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
  region <- Region()
  population <- StagedPopulation(region, growth = stage_matrix, capacity = 300,
                                 capacity_area = 1e+06, capacity_stages = 2:3)
  n <- population$make(initial = 350)
  r <- exp(log(population$get_growth_r())*(1 - sum(n[2:3])/300))
  r_mult <- get("r_mult", envir = environment(population$grow))
  mult <- r_mult$mult[which.min(abs(r_mult$r - r))]
  survivals <- attr(stage_matrix, "survivals")
  reproductions <- stage_matrix - survivals
  attr(reproductions, "survivals") <- NULL
  set.seed(1243)
  new_n <- array(0L, c(1, 3))
  for (stage in 1:3) {
    new_n <- new_n +
      stats::rpois(3, (n[1, stage]*t(reproductions[, stage])*mult))
    stage_surv <- stats::rbinom(1, n[1, stage], sum(survivals[, stage])*mult)
    new_n <- new_n + t(stats::rmultinom(1, size = stage_surv,
                                        prob = survivals[, stage]*mult))
  }
  colnames(new_n) <- c("stage 1", "stage 2", "stage 3")
  set.seed(1243)
  expect_equal(population$grow(n), new_n)
  attr(n, "diffusion_rate") <- 2000
  attr(n, "diffusion_radius") <- 1000
  r <- exp(log(population$get_growth_r())*
             (1 - sum(n[2:3])/(300*pi*3000^2/1e+06)))
  mult <- r_mult$mult[which.min(abs(r_mult$r - r))]
  set.seed(1243)
  new_n <- array(0L, c(1, 3))
  for (stage in 1:3) {
    new_n <- new_n +
      stats::rpois(3, (n[1, stage]*t(reproductions[, stage])*mult))
    stage_surv <- stats::rbinom(1, n[1, stage], sum(survivals[, stage])*mult)
    new_n <- new_n + t(stats::rmultinom(1, size = stage_surv,
                                        prob = survivals[, stage]*mult))
  }
  colnames(new_n) <- c("stage 1", "stage 2", "stage 3")
  attr(new_n, "diffusion_rate") <- 2000
  attr(new_n, "diffusion_radius") <- 1000
  set.seed(1243)
  expect_equal(population$grow(n), new_n)
})
