context("Removals")

test_that("initializes with region, population, and other parameters", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  template_vect <- template[region$get_indices()][,1]
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  expect_error(Removals(region, population_model,
                        removal_pr = (0:10)/10),
               paste("Removal probability should be a vector with a value 0-1",
                     "for each region location."))
  expect_error(Removals(region, population_model,
                        removal_pr = 2),
               paste("Removal probability should be a vector with a value 0-1",
                     "for each region location."))
  expect_error(removals <- Removals(region, population_model,
                                    removal_pr = 0.5,
                                    removal_cost = 1:5),
               paste("The removal cost parameter must be a numeric vector",
                     "with values for each location."))
  expect_error(removals <- Removals(region, population_model,
                                    radius = -1),
               "The radius (m) parameter must be numeric and >= 0.",
               fixed = TRUE)
  expect_message(removals <- Removals(region, population_model,
                                      detected_only = TRUE,
                                      radius = 2000),
                 paste("Radius is not used when only detected individuals are",
                       "removed."))

  expect_silent(removals <- Removals(region, population_model))
  expect_silent(removals <- Removals(region, population_model,
                                     removal_pr = template_vect,
                                     removal_cost = 2,
                                     radius = 2000,
                                     stages = 2:3,
                                     schedule = 4:6))
  expect_is(removals, "Removals")
  expect_s3_class(removals, "Actions")
  expect_named(removals,
               c(c("get_type", "get_id", "set_id", "get_label", "get_stages",
                   "get_schedule", "include_cost", "get_cost_label",
                   "get_cost_unit", "get_attributes", "clear_attributes",
                   "apply", "get_removal_pr_type")))
  expect_equal(removals$get_type(), "removal")
  expect_equal(removals$get_label(), "removed")
  expect_equal(removals$get_removal_pr_type(), "individual")
  expect_equal(removals$get_stages(), 2:3)
  expect_equal(removals$get_schedule(), 4:6)
  expect_true(removals$include_cost())
  expect_equal(removals$get_cost_label(), "removal_cost")
  expect_null(removals$get_cost_unit())
  removal_cost <- 2
  attr(removal_cost, "unit") <- "beans"
  removals <- Removals(region, population_model,
                       removal_pr = template_vect,
                       removal_cost = removal_cost) # silent
  expect_equal(removals$get_cost_unit(), "beans")
  expect_silent(removals$set_id(1))
  expect_equal(removals$get_id(), 1)
  expect_equal(removals$get_label(), "1_removed")
  expect_equal(removals$get_label(include_id = FALSE), "removed")
  expect_equal(removals$get_cost_label(), "1_removal_cost")
  expect_equal(removals$get_cost_label(include_id = FALSE),
               "removal_cost")
})

test_that("applies stochastic removals to invasive population", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  idx1 <- 5916:5922
  region <- Region(template)
  template[region$get_indices()][idx1,] <- c(0.5, 0.75, 1, 0.5, 0.5, 0.75, 1)
  idx <- idx1[5:7]
  template_vect <- template[region$get_indices()][,1]
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  initial_n <- rep(0, region$get_locations())
  initial_n[idx] <- (10:12)*10
  initializer <- Initializer(initial_n, region = region,
                             population_model = population_model)
  # without detected
  set.seed(1234)
  n <- initializer$initialize()
  set.seed(1234)
  expected_removals <- array(c(rep(0, 3),
                               stats::rbinom(6, size = n[idx,2:3],
                                             c(0.5, 0.75, 1))), c(3, 3))
  colnames(expected_removals) <- colnames(n)
  expected_removed <- rep(FALSE, region$get_locations())
  expected_removed[idx] <-
    rowSums((n[idx, 2:3] - expected_removals[,2:3])) == 0
  removal_cost <- 2
  attr(removal_cost, "unit") <- "$"
  expect_silent(removals <- Removals(region, population_model,
                                     removal_pr = template_vect,
                                     detected_only = FALSE,
                                     removal_cost = removal_cost,
                                     radius = NULL,
                                     stages = 2:3,
                                     schedule = 4:6))
  expected_removal_cost <- rep(0, region$get_locations())
  attr(expected_removal_cost, "unit") <- "$"
  expect_silent(new_n <- removals$apply(n, 2)) # not scheduled
  expect_equal(attr(attr(new_n, "removed"), "number")[idx,],
               expected_removals*0)
  expect_equal(as.logical(attr(new_n, "removed")), expected_removed & FALSE)
  expect_equal(new_n[idx,], n[idx,])
  expect_equal(attr(new_n, "removal_cost"), expected_removal_cost)
  set.seed(1234)
  expect_silent(new_n <- removals$apply(n, 4))
  expect_equal(attr(attr(new_n, "removed"), "number")[idx,], expected_removals)
  expect_equal(as.logical(attr(new_n, "removed")), expected_removed)
  expect_equal(new_n[idx,], n[idx,] - expected_removals)
  expected_removal_cost[idx] <- 2
  expect_equal(attr(new_n, "removal_cost"), expected_removal_cost)
  expect_equal(attr(attr(new_n, "removal_cost"), "unit"), "$")
  expect_equal(removals$get_attributes(new_n),
               list(removed = attr(new_n, "removed"),
                    removal_cost = attr(new_n, "removal_cost")))
  expect_silent(new_n <- removals$clear_attributes(new_n))
  expect_null(attr(new_n, "removed"))
  expect_null(attr(new_n, "removal_cost"))
  # single value removal_pr
  set.seed(1234)
  expected_removals <- array(c(rep(0, 3),
                               stats::rbinom(6, size = n[idx,2:3],
                                             0.65)), c(3, 3))
  colnames(expected_removals) <- colnames(n)
  expected_removed <- rep(FALSE, region$get_locations())
  expected_removed[idx] <-
    rowSums((n[idx, 2:3] - expected_removals[,2:3])) == 0
  expected_removal_cost[] <- 0
  expected_removal_cost[idx] <- 2
  expect_silent(removals <- Removals(region, population_model,
                                     removal_pr = 0.65,
                                     detected_only = FALSE,
                                     removal_cost = removal_cost,
                                     radius = NULL,
                                     stages = 2:3, schedule = 4:6))
  set.seed(1234)
  expect_silent(new_n <- removals$apply(n, 4))
  expect_equal(attr(attr(new_n, "removed"), "number")[idx,], expected_removals)
  expect_equal(as.logical(attr(new_n, "removed")), expected_removed)
  expect_equal(attr(new_n, "removal_cost"), expected_removal_cost)
  # detection only
  expect_silent(removals <- Removals(region, population_model,
                                     removal_pr = template_vect,
                                     detected_only = TRUE,
                                     removal_cost = removal_cost,
                                     radius = NULL,
                                     stages = 2:3,
                                     schedule = 4:6))
  set.seed(1234)
  expect_silent(new_n <- removals$apply(n, 4))
  expect_equal(attr(attr(new_n, "removed"), "number")[idx,], n[idx,]*0)
  expect_equal(as.logical(attr(new_n, "removed")), expected_removed & FALSE)
  expect_equal(new_n[idx,], n[idx,])
  expect_equal(attr(new_n, "removal_cost"), expected_removal_cost*0)
  # with detected/undetected
  detected <- n*0
  detected[idx,2:3] <- trunc(n[idx,2:3]*c(0, 0.5, 1))
  undetected <- n - detected
  attr(n, "detected") <- detected
  attr(n, "undetected") <- undetected
  set.seed(1234)
  expected_removals <- array(
    c(rep(0, 3), stats::rbinom(6, size = attr(n, "detected")[idx,2:3],
                               c(0.5, 0.75, 1))), c(3, 3))
  colnames(expected_removals) <- colnames(n)
  expected_removed <- rep(FALSE, region$get_locations())
  expected_removed[idx] <-
    rowSums((n[idx, 2:3] - expected_removals[,2:3])) == 0
  attr(n, "attachment") <- "extra"
  set.seed(1234)
  expect_silent(new_n <- removals$apply(n, 4))
  expect_equal(attr(attr(new_n, "removed"), "number")[idx,], expected_removals)
  expect_equal(as.logical(attr(new_n, "removed")), expected_removed)
  expect_equal(new_n[idx,], n[idx,] - expected_removals)
  expect_equal(attr(new_n, "attachment"), "extra")
  expect_equal(attr(new_n, "undetected"), attr(n, "undetected"))
  expect_equal(attr(new_n, "removal_cost"),
               expected_removal_cost*(rowSums(detected) > 0))
  attr(n, "attachment") <- NULL
  expect_silent(removals <- Removals(region, population_model,
                                     removal_pr = template_vect,
                                     detected_only = FALSE,
                                     removal_cost = removal_cost,
                                     radius = NULL,
                                     stages = 2:3))
  n_apply <- array(as.numeric(n*(detected > 0)), dim(n))
  n_undetected <- array(as.numeric(attr(n, "undetected")*(n_apply > 0)),
                        dim(n))
  n_apply <- list(detected = n_apply - n_undetected, undetected = n_undetected)
  set.seed(1234)
  expected_removals <- lapply(n_apply, function(a) {
    removed <- array(c(rep(0, 3),
                       stats::rbinom(6, size = a[idx,2:3], c(0.5, 0.75, 1))),
                     c(3, 3))
    colnames(removed) <- colnames(n)
    return(removed)
  })
  expected_removed <- rep(FALSE, region$get_locations())
  expected_removed[idx] <-
    rowSums((n[idx, 2:3] - (expected_removals$detected +
                              expected_removals$undetected)[,2:3])) == 0
  set.seed(1234)
  expect_silent(new_n <- removals$apply(n, 4))
  expect_equal(attr(attr(new_n, "removed"), "number")[idx,],
               expected_removals$detected + expected_removals$undetected)
  expect_equal(as.logical(attr(new_n, "removed")), expected_removed)
  expect_equal(new_n[idx,], n[idx,] - (expected_removals$detected +
                                         expected_removals$undetected))
  expect_equal(attr(new_n, "undetected")[idx,],
               attr(n, "undetected")[idx,] - expected_removals$undetected)
  expect_equal(attr(new_n, "removal_cost"),
               expected_removal_cost*(rowSums(detected) > 0))
  expect_silent(removals$set_id(4))
  set.seed(1234)
  expect_silent(new_n <- removals$apply(n, 4))
  expect_equal(attr(attr(new_n, "4_removed"), "number")[idx,],
               expected_removals$detected + expected_removals$undetected)
  expect_equal(as.logical(attr(new_n, "4_removed")), expected_removed)
  expect_equal(new_n[idx,], n[idx,] - (expected_removals$detected +
                                         expected_removals$undetected))
  expect_equal(attr(new_n, "undetected")[idx,],
               attr(n, "undetected")[idx,] - expected_removals$undetected)
  expect_equal(attr(new_n, "4_removal_cost"),
               expected_removal_cost*(rowSums(detected) > 0))
  # remove always (as per without detected)
  expect_silent(removals <-
                  Removals(region, population_model,
                           removal_pr = template_vect,
                           remove_always = TRUE,
                           removal_cost = removal_cost,
                           radius = NULL,
                           stages = 2:3, schedule = 4:6))
  n_apply <- array(as.numeric(n), dim(n))
  n_undetected <- array(as.numeric(attr(n, "undetected")*(n_apply > 0)),
                        dim(n))
  n_apply <- list(detected = n_apply - n_undetected, undetected = n_undetected)
  set.seed(1234)
  expected_removals <- lapply(n_apply, function(a) {
    removed <- array(c(rep(0, 3),
                       stats::rbinom(6, size = a[idx,2:3], c(0.5, 0.75, 1))),
                     c(3, 3))
    colnames(removed) <- colnames(n)
    return(removed)
  })
  expected_removed <- rep(FALSE, region$get_locations())
  expected_removed[idx] <-
    rowSums((n[idx, 2:3] - (expected_removals$detected +
                              expected_removals$undetected)[,2:3])) == 0
  set.seed(1234)
  expect_silent(new_n <- removals$apply(n, 4))
  expect_equal(attr(attr(new_n, "removed"), "number")[idx,],
               expected_removals$detected + expected_removals$undetected)
  expect_equal(as.logical(attr(new_n, "removed")), expected_removed)
  expect_equal(new_n[idx,], n[idx,] - (expected_removals$detected +
                                         expected_removals$undetected))
  expect_equal(attr(new_n, "undetected")[idx,],
               attr(n, "undetected")[idx,] - expected_removals$undetected)
  expect_equal(attr(new_n, "removal_cost"),
               expected_removal_cost*(rowSums(n) > 0))
  # with radius
  expect_silent(removals <- Removals(region, population_model,
                                     removal_pr = template_vect,
                                     detected_only = FALSE,
                                     removal_cost = removal_cost,
                                     radius = 3000,
                                     stages = 2:3))
  n[idx1[1:4],] <- rep(n[idx1[5],] - 6, each = 4)
  attr(n, "undetected")[idx1[1:4],] <- n[idx1[1:4],]
  n_apply <- array(as.numeric(n), dim(n))
  n_undetected <- array(as.numeric(attr(n, "undetected")*(n_apply > 0)),
                        dim(n))
  n_apply <- list(detected = n_apply - n_undetected, undetected = n_undetected)
  set.seed(1234)
  expected_removals <- lapply(n_apply, function(a) {
    removed <- rbind(array(0, c(2, 3)), array(
      c(rep(0, 5), stats::rbinom(10, size = a[idx1[3:7], 2:3],
                                 template_vect[idx1[3:7]])), c(5, 3)))
    colnames(removed) <- colnames(n)
    return(removed)
  })
  expected_removed <- rep(FALSE, region$get_locations())
  expected_removed[idx1] <-
    rowSums((n[idx1, 2:3] - (expected_removals$detected +
                               expected_removals$undetected)[,2:3])) == 0
  set.seed(1234)
  expect_silent(new_n <- removals$apply(n, 4))
  expect_equal(attr(attr(new_n, "removed"), "number")[idx1,],
               expected_removals$detected + expected_removals$undetected)
  expect_equal(as.logical(attr(new_n, "removed")), expected_removed)
  expect_equal(new_n[idx1,], n[idx1,] - (expected_removals$detected +
                                           expected_removals$undetected))
  expect_equal(attr(new_n, "undetected")[idx1,],
               attr(n, "undetected")[idx1,] - expected_removals$undetected)
  expected_removal_cost[] <- 0
  expected_removal_cost[rowSums(attr(attr(new_n, "removed"),
                                     "number")) > 0] <- 2
  expect_equal(attr(new_n, "removal_cost"), expected_removal_cost)
  # population level removal probability
  expect_silent(
    removals <- Removals(region, population_model,
                         removal_pr = template_vect,
                         removal_pr_type = "population",
                         detected_only = FALSE,
                         stages = 2:3))
  idx2 <- which(rowSums(n[,2:3]) > 0)
  attr(n, "undetected") <- NULL
  removed <- zeroed <- rep(0, length(idx2))
  set.seed(1234)
  for (i in 1:1000) {
    new_n <- removals$apply(n, 4)
    removed <- removed + attr(new_n, "removed")[idx2]
    zeroed <- zeroed + (rowSums(new_n[idx2, 2:3]) == 0)
  }
  expect_true(all(abs(removed/1000 - template_vect[idx2]) < 0.05))
  expect_equal(zeroed, removed)
  # unstructured population
  population_model <- UnstructPopulation(region, growth = 1.2)
  set.seed(1234)
  n <- rowSums(initializer$initialize())
  set.seed(1234)
  expected_removals <- stats::rbinom(3, size = n[idx], c(0.5, 0.75, 1))
  expected_removed <- rep(FALSE, region$get_locations())
  expected_removed[idx] <- n[idx] - expected_removals == 0
  removal_cost <- 2
  attr(removal_cost, "unit") <- "$"
  expect_silent(
    removals <- Removals(region, population_model,
                         removal_pr = template_vect,
                         detected_only = FALSE,
                         removal_cost = removal_cost,
                         radius = NULL,
                         schedule = 4:6))
  expected_removal_cost <- rep(0, region$get_locations())
  attr(expected_removal_cost, "unit") <- "$"
  set.seed(1234)
  expect_silent(new_n <- removals$apply(n, 4))
  expect_equal(attr(attr(new_n, "removed"), "number")[idx], expected_removals)
  expect_equal(as.logical(attr(new_n, "removed")), expected_removed)
  expect_equal(new_n[idx], n[idx] - expected_removals)
  expected_removal_cost[idx] <- 2
  expect_equal(attr(new_n, "removal_cost"), expected_removal_cost)
  # population level removal probability
  expect_silent(
    removals <- Removals(region, population_model,
                         removal_pr = template_vect,
                         removal_pr_type = "population",
                         detected_only = FALSE))
  idx2 <- which(n > 0)
  attr(n, "undetected") <- NULL
  removed <- zeroed <- rep(0, length(idx2))
  set.seed(1234)
  for (i in 1:1000) {
    new_n <- removals$apply(n, 4)
    removed <- removed + attr(new_n, "removed")[idx2]
    zeroed <- zeroed + (new_n[idx2] == 0)
  }
  expect_true(all(abs(removed/1000 - template_vect[idx2]) < 0.05))
  expect_equal(zeroed, removed)
  # presence-only population
  population_model <- PresencePopulation(region)
  n <- n > 0
  set.seed(121)
  expected_removals <- stats::rbinom(3, size = n[idx], c(0.5, 0.75, 1))
  removal_cost <- 2
  attr(removal_cost, "unit") <- "$"
  expect_silent(
    removals <- Removals(region, population_model,
                         removal_pr = template_vect,
                         detected_only = FALSE,
                         removal_cost = removal_cost,
                         radius = NULL,
                         schedule = 4:6))
  expected_removal_cost <- rep(0, region$get_locations())
  attr(expected_removal_cost, "unit") <- "$"
  set.seed(121)
  expect_silent(new_n <- removals$apply(n, 4))
  expect_equal(attr(new_n, "removed")[idx], as.logical(expected_removals))
  expect_equal(new_n[idx], n[idx] & !expected_removals)
  expected_removal_cost[idx] <- 2
  expect_equal(attr(new_n, "removal_cost"), expected_removal_cost)
  # spatially-implicit
  region <- Region()
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  initializer <- Initializer(50, region = region,
                             population_model = population_model)
  set.seed(1234)
  n <- initializer$initialize()
  expect_silent(
    removals <- Removals(region, population_model,
                         removal_pr = 0.65,
                         detected_only = FALSE,
                         removal_cost = 2,
                         radius = NULL,
                         stages = 2:3, schedule = 4:6))
  expect_error(new_n <- removals$apply(n, 4),
               paste("Cannot calculate spatially implicit removal cost",
                     "without area occupied."))
  attr(n, "spread_area") <- 300
  set.seed(1234)
  expected_removals <-
    array(c(0, stats::rbinom(2, size = n[,2:3], 0.65)), c(1, 3))
  expected_removed <- rowSums((n - expected_removals)[,2:3, drop = FALSE]) == 0
  colnames(expected_removals) <- colnames(n)
  set.seed(1234)
  expect_silent(new_n <- removals$apply(n, 4))
  expect_equal(new_n[,], n[,] - expected_removals[,])
  expect_equal(attr(attr(new_n, "removed"), "number"), expected_removals)
  expect_equal(as.logical(attr(new_n, "removed")), expected_removed)
  expect_equal(attr(new_n, "removal_cost"), 2*300)
})
