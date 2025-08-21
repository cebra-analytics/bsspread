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
               paste("Population capacity should be a vector or matrix with",
                     "a value or row for each region location."))

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

test_that("makes populations with initial ages and first stages", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  attr(stage_matrix, "labels") <- c("a", "b", "c")
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix))
  initial_n <- +(template[region$get_indices()][,1] > 0.5)
  idx <- which(initial_n > 0)
  initial_n[idx] <- stats::rpois(length(idx), 20)
  initial_age <- initial_n*0
  initial_age[idx[1:150]] <- c(rep(3,50), rep(2,50), rep(1,50))
  attr(initial_n, "age") <- initial_age
  attr(initial_n, "stages") <- 1
  expect_silent(n <- population$make(initial = initial_n))
  expect_equal(dim(n), c(region$get_locations(), 3))
  expect_equal(colnames(n), c("a", "b", "c"))

  expect_equal(rowSums(n[idx,]), initial_n[idx])
  expect_true(all(rowSums(n[-idx,]) == 0))
  expect_equal(colSums(n[idx[1:50],]) > 0, c(a = TRUE, b = TRUE, c = TRUE))
  expect_equal(sum(n[idx[1:50],]), sum(initial_n[idx[1:50]]))
  expect_equal(colSums(n[idx[51:100],]) > 0, c(a = TRUE, b = FALSE, c = TRUE))
  expect_equal(sum(n[idx[51:100],]), sum(initial_n[idx[51:100]]))
  expect_equal(colSums(n[idx[101:150],]),
               c(a = 0, b = sum(initial_n[idx[101:150]]), c = 0))
  expect_equal(colSums(n[idx[151:length(idx)],]),
               c(a = sum(initial_n[idx[151:length(idx)]]), b = 0, c = 0))
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
  expect_equal(round(mean(rowSums(n)[which(incursion)])), 10)
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix))
  expect_silent(population$set_incursion_mean(15))
  expect_silent(population$set_incursion_stages(1:2))
  expect_silent(n <- population$make(incursion = incursion))
  expect_true(all(n[,3] == 0))
  expect_equal(round(mean(rowSums(n)[which(incursion)])), 15)
})

test_that("makes populations with variable growth", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  expect_error(StagedPopulation(region, growth = stage_matrix,
                                growth_mult = as.matrix(1:5)),
               paste("Growth multiplier should be a matrix with a single row",
                     "or a row for each region location."))
  expect_error(StagedPopulation(region, growth = stage_matrix,
                          growth_mult = t(as.matrix(-1:5))),
               "Growth multiplier values should be >= 0 and <= 1.")
  growth_mult_s <- as.matrix(template[region$get_indices()][,1])
  growth_mult_st <- cbind(growth_mult_s, growth_mult_s*0.8, growth_mult_s*0.6)
  growth_mult_t <- growth_mult_st[1,,drop = F]
  attr(growth_mult_s, "apply_to") <- "dummy"
  expect_error(StagedPopulation(region, growth = stage_matrix,
                                growth_mult = growth_mult_s),
               paste("Growth multiplier 'apply to' attribute should be",
                     "'reproductions' or 'survivals'."))
  attr(growth_mult_s, "apply_to") <- "reproductions"
  attr(growth_mult_s, "stages") <- 2:4
  expect_error(StagedPopulation(region, growth = stage_matrix,
                                growth_mult = growth_mult_s),
               paste("Growth multiplier stages should specify index values",
                     "between 1 and 3."))
  attr(growth_mult_s, "stages") <- 2:3
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                         growth_mult = growth_mult_s))
  expect_equal(population$get_growth_mult(cells = 1:10), growth_mult_s[1:10,])
  expect_equal(population$get_growth_mult(tm = 2), as.numeric(growth_mult_s))
  expect_equal(population$get_growth_mult(cells = 1:10, tm = 2),
               growth_mult_s[1:10,])
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                         growth_mult = growth_mult_t))
  expect_equal(population$get_growth_mult(cells = 1:10), growth_mult_t[1])
  expect_equal(population$get_growth_mult(tm = 2), growth_mult_t[2])
  expect_equal(population$get_growth_mult(cells = 1:10, tm = 2),
               growth_mult_t[2])
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                         growth_mult = growth_mult_st))
  expect_equal(population$get_growth_mult(cells = 1:10),
               growth_mult_st[1:10,1])
  expect_equal(population$get_growth_mult(tm = 2), growth_mult_st[,2])
  expect_equal(population$get_growth_mult(cells = 1:10, tm = 2),
               growth_mult_st[1:10,2])
  expect_equal(population$get_growth_mult(cells = 1:10, tm = 5),
               growth_mult_st[1:10,2])
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
  set.seed(1243)
  expect_silent(n <- population$make(initial = initial))
  idx <- which(rowSums(n) > 0)
  expected_r <- Re((eigen(stage_matrix,)$values)[1])
  set.seed(1243)
  expect_silent(n1 <- population$grow(n, 1))
  expect_true(abs(mean(rowSums(n1[idx,]))/20 - expected_r) < 0.05)
  # growth variation
  growth_mult <- rep(1, region$get_locations())
  growth_mult <- cbind(growth_mult, growth_mult*0.8, growth_mult*0.6,
                       growth_mult*0)
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                               growth_mult = growth_mult))
  set.seed(1243)
  expect_equal(mean(rowSums(population$grow(n, 1)[idx,]))/20,
               mean(rowSums(n1[idx,]))/20)
  expect_true(abs(round(mean(rowSums(population$grow(n, 2)[idx,]))/20, 2) -
                    expected_r*0.8) <= 0.02)
  expect_true(abs(round(mean(rowSums(population$grow(n, 3)[idx,]))/20, 2) -
                    expected_r*0.6) <= 0.02)
  expect_true(all(rowSums(population$grow(n, 4)[idx,]) == 0))
  expect_true(abs(round(mean(rowSums(population$grow(n, 6)[idx,]))/20, 2) -
                    expected_r*0.8) <= 0.02)

  attr(growth_mult, "apply_to") <- "reproductions"
  growth_mult[1:10,2] <- 0
  idx1 <- idx[idx <= 10]
  idx2 <- idx[idx > 10]
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                               growth_mult = growth_mult))
  set.seed(1243)
  expect_equal(mean(rowSums(population$grow(n, 1)[idx,]))/20,
               mean(rowSums(n1[idx,]))/20)
  expect_silent(n2 <- population$grow(n, 2))
  expect_equal(sum(n2[idx1,1]), 0)
  expect_equal(round(sum(n2[,2])/sum(n[idx2, 1]), 2), stage_matrix[2,1])
  expect_true(
    abs((sum(n[idx2, 2])*stage_matrix[3,2] +
           sum(n[idx2, 3])*stage_matrix[3,3])/sum(n2[,3]) - 1) < 0.01)
  expect_true(
    abs((sum(n[idx2, 2])*stage_matrix[1,2]*0.8 +
           sum(n[idx2, 3])*stage_matrix[1,3]*0.8)/sum(n2[,1]) - 1) < 0.01)
  expect_silent(n3 <- population$grow(n, 3))
  expect_equal(round(sum(n3[,2])/sum(n[idx2, 1]), 2), stage_matrix[2,1])
  expect_true(
    abs((sum(n[,2])*stage_matrix[3,2] +
           sum(n[,3])*stage_matrix[3,3])/sum(n3[,3]) - 1) < 0.01)
  expect_true(
    abs((sum(n[,2])*stage_matrix[1,2]*0.6 +
           sum(n[,3])*stage_matrix[1,3]*0.6)/sum(n3[,1]) - 1) < 0.01)
  attr(growth_mult, "apply_to") <- "survivals"
  attr(growth_mult, "stages") <- 2:3
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                               growth_mult = growth_mult))
  set.seed(1243)
  expect_equal(mean(rowSums(population$grow(n, 1)[idx,]))/20,
               mean(rowSums(n1[idx,]))/20)
  expect_silent(n2 <- population$grow(n, 2))
  expect_equal(round(sum(n2[,2])/sum(n[,1]), 2), stage_matrix[2,1])
  expect_equal(sum(n2[idx1, 3]), 0)
  expect_equal(round(sum(n2[,2])/sum(n[,1]), 2), stage_matrix[2,1])
  expect_true(
    abs((sum(n[idx2, 2])*stage_matrix[3,2]*0.8 +
           sum(n[idx2, 3])*stage_matrix[3,3]*0.8)/sum(n2[idx2, 3]) - 1) < 0.01)
  expect_true(
    abs((sum(n[idx2, 2])*stage_matrix[1,2] +
           sum(n[idx2, 3])*stage_matrix[1,3])/sum(n2[idx2, 1]) - 1) < 0.01)
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
  capacity <- cbind(capacity, capacity*1.5)
  expect_error(population <- StagedPopulation(region, growth = stage_matrix,
                                              capacity = capacity,
                                              capacity_stages = 2:4),
               "Capacity stages should specify index values between 1 and 3.")
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                               capacity = capacity,
                                               capacity_stages = 2:3))

  idx <- which(template[region$get_indices()] > 0)
  initial <- rep(0, region$get_locations())
  initial[idx] <- stats::rpois(length(idx), capacity[idx, 1]*3)
  set.seed(1243)
  n <- population$make(initial = initial) # silent
  idx <- which(rowSums(n) > 0)
  set.seed(1243)
  expect_silent(n1 <- population$grow(n, 1))
  mean_growth_1 <- mean(rowSums(n1[idx,])/rowSums(n[idx,]))
  expect_equal(round(mean_growth_1, 1), 1.0)
  set.seed(1243)
  expect_silent(n2 <- population$grow(n, 2))
  mean_growth_2 <- mean(rowSums(n2[idx,])/rowSums(n[idx,]))
  expect_true(mean_growth_2 > mean_growth_1)
  set.seed(1243)
  expect_silent(n3 <- population$grow(n, 3))
  expect_equal(n3, n1)
  # growth variation
  expected_r <- Re((eigen(stage_matrix,)$values)[1])
  growth_mult <- rep(1, region$get_locations())
  growth_mult <- cbind(growth_mult, growth_mult*0.8, growth_mult*0.9,
                       growth_mult*0)
  capacity <- template[region$get_indices()][,1]*10
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                               growth_mult = growth_mult,
                                               capacity = capacity,
                                               capacity_stages = 2:3))
  set.seed(1243)
  expect_equal(mean(rowSums(population$grow(n, 1)[idx,])/rowSums(n[idx,])),
               mean_growth_1)
  set.seed(1243)
  expect_silent(mean_growth_2 <- mean(rowSums(population$grow(n, 2)[idx,])/
                                        rowSums(n[idx,])))
  expect_true(mean_growth_2 < expected_r*0.8)
  expect_true(mean(rowSums(population$grow(n, 3)[idx,])/rowSums(n[idx,])) <
                expected_r*0.9)
  expect_true(all(rowSums(population$grow(n, 4)[idx,]) == 0))
  set.seed(1243)
  expect_equal(mean(rowSums(population$grow(n, 6)[idx,])/rowSums(n[idx,])),
               mean_growth_2)
  attr(growth_mult, "apply_to") <- "reproductions"
  growth_mult[1:10,2] <- 0
  idx1 <- idx[idx <= 10]
  idx2 <- idx[idx > 10]
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                               growth_mult = growth_mult,
                                               capacity = capacity,
                                               capacity_stages = 2:3))
  set.seed(1243)
  expect_equal(mean(rowSums(population$grow(n, 1)[idx,])),
               mean(rowSums(n1[idx,])))
  set.seed(1243)
  expect_silent(n2 <- population$grow(n, 2))
  mod_matrix <- stage_matrix
  mod_matrix[1,2:3] <- mod_matrix[1,2:3]*0.8
  mod_r <- Re((eigen(mod_matrix,)$values)[1])
  mult <- exp(log(mod_r)*(1 - sum(n[idx2, 2:3])/sum(capacity[idx2])))/mod_r
  expect_equal(sum(n2[idx1,1]), 0)
  expect_true(abs(sum(n2[idx2,1])/
                    (sum(n[idx2, 2])*stage_matrix[1,2]*mult*0.8 +
                       sum(n[idx2, 3])*stage_matrix[1,3]*mult*0.8) - 1) < 0.1)
  expect_true(abs(sum(n2[idx2,2])/
                    (sum(n[idx2, 1])*stage_matrix[2,1]*mult) - 1) < 0.1)
  expect_true(abs(sum(n2[idx2,3])/
                    (sum(n[idx2, 2])*stage_matrix[3,2]*mult +
                       sum(n[idx2, 3])*stage_matrix[3,3]*mult) - 1) < 0.1)
  attr(growth_mult, "apply_to") <- "survivals"
  attr(growth_mult, "stages") <- 2:3
  expect_silent(population <- StagedPopulation(region, growth = stage_matrix,
                                               growth_mult = growth_mult,
                                               capacity = capacity,
                                               capacity_stages = 2:3))
  set.seed(1243)
  expect_equal(mean(rowSums(population$grow(n, 1)[idx,]))/20,
               mean(rowSums(n1[idx,]))/20)
  set.seed(1243)
  expect_silent(n2 <- population$grow(n, 2))
  mod_matrix <- stage_matrix
  mod_matrix[3,2:3] <- mod_matrix[3,2:3]*0.8
  mod_r <- Re((eigen(mod_matrix,)$values)[1])
  mult <- exp(log(mod_r)*(1 - sum(n[idx2, 2:3])/sum(capacity[idx2])))/mod_r
  expect_equal(sum(n2[idx1, 3]), 0)
  expect_true(abs(sum(n2[idx2,1])/
                    (sum(n[idx2, 2])*stage_matrix[1,2]*mult +
                       sum(n[idx2, 3])*stage_matrix[1,3]*mult) - 1) < 0.1)
  expect_true(abs(sum(n2[idx2,2])/
                    (sum(n[idx2,1])*stage_matrix[2,1]*mult) - 1) < 0.1)
  expect_true(abs(sum(n2[idx2,3])/
                    (sum(n[idx2, 2])*stage_matrix[3,2]*mult*0.8 +
                       sum(n[idx2, 3])*stage_matrix[3,3]*mult*0.8) - 1) < 0.1)
  region <- Region()
  population <- StagedPopulation(region, growth = stage_matrix, capacity = 300,
                                 capacity_area = 1e+06, capacity_stages = 2:3)
  n <- population$make(initial = 350)
  survivals <- attr(stage_matrix, "survivals")
  reproductions <- stage_matrix - survivals
  attr(reproductions, "survivals") <- NULL
  set.seed(1243)
  new_n <- array(0L, c(1, 3))
  for (stage in 1:3) {
    new_n <- new_n + stats::rpois(3, (n[1, stage]*t(reproductions[, stage])))
    stage_surv <- stats::rbinom(1, n[1, stage], sum(survivals[, stage]))
    new_n <- new_n + t(stats::rmultinom(1, size = stage_surv,
                                        prob = survivals[, stage]))
  }
  colnames(new_n) <- c("stage 1", "stage 2", "stage 3")
  set.seed(1243)
  expect_equal(population$grow(n, 1), new_n)
  region$set_max_implicit_area(1e+06)
  r <- exp(log(population$get_growth_r())*(1 - sum(n[2:3])/300))
  mult <- r/population$get_growth_r()
  set.seed(1243)
  new_n <- array(0L, c(1, 3))
  for (stage in 1:3) {
    new_n <-
      new_n + stats::rpois(3, (n[1, stage]*t(reproductions[, stage])*mult))
    stage_surv <- stats::rbinom(1, n[1, stage], sum(survivals[, stage])*mult)
    new_n <- new_n + t(stats::rmultinom(1, size = stage_surv,
                                        prob = survivals[, stage]*mult))
  }
  colnames(new_n) <- c("stage 1", "stage 2", "stage 3")
  set.seed(1243)
  expect_equal(population$grow(n, 1), new_n)
  attr(n, "diffusion_rate") <- 2000
  attr(n, "diffusion_radius") <- 1000
  r <- exp(log(population$get_growth_r())*
             (1 - sum(n[2:3])/(300*pi*3000^2/1e+06)))
  mult <- r/population$get_growth_r()
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
  expect_equal(population$grow(n, 1), new_n)
})
