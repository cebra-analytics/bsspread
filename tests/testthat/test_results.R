context("Results")

test_that("initializes with region, population model, and other parameters", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  expect_error(results <- Results(region, population_model = 0),
               paste("Population model must be a 'Population' or inherited",
                     "class object."))
  expect_silent(results <- Results(region, population_model = population))
  expect_is(results, "Results")
  expect_equal(results$get_params(),
               list(time_steps = 1, step_duration = 1, step_units = "years",
                    collation_steps = 1, replicates = 1, stages = NULL,
                    combine_stages = NULL))
})

test_that("collates single replicate results", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  expect_silent(results <- Results(region, population_model = population,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 1,
                                   combine_stages = NULL))
  expect_silent(result_list <- results$get_list())
  expect_equal(names(result_list), c("occupancy", "total_occup", "area"))
  expect_named(result_list$occupancy, as.character(seq(0, 10, 2)))
  expect_true(all(
    unlist(lapply(result_list$occupancy, length)) == region$get_locations()))
  expect_true(all(unlist(result_list$occupancy) == 0))
  expect_named(result_list$total_occup, as.character(0:10))
  expect_true(all(unlist(lapply(result_list$total_occup, length)) == 1))
  expect_true(all(unlist(result_list$total_occup) == 0))
  expect_named(result_list$area, as.character(0:10))
  expect_true(all(unlist(lapply(result_list$area, length)) == 1))
  expect_true(all(unlist(result_list$area) == 0))
  n <- rep(0, region$get_locations())
  for (tm in 0:10) {
    n[tm + 1] <- 1
    results$collate(r = 1, tm, n)
  }
  expect_silent(result_list <- results$get_list())
  expected <- lapply(seq(0, 10, 2),
                     function(tm) c(rep(1, tm + 1), rep(0, 11 - tm)))
  names(expected) <- as.character(seq(0, 10, 2))
  expect_equal(lapply(result_list$occupancy, function(rl) rl[1:12]), expected)
  expect_equal(unname(unlist(result_list$total_occup)), 1:11)
  expect_equal(signif(unname(unlist(result_list$area)), 5), (1:11)*1000^2)
  expect_equal(attr(result_list$area, "units"), "square metres")
  # unstructured population
  population <- UnstructPopulation(region)
  expect_silent(results <- Results(region, population_model = population,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 1,
                                   combine_stages = NULL))
  n <- rep(0, region$get_locations())
  for (tm in 0:10) {
    n[tm + 1] <- 1
    results$collate(r = 1, tm, n)
  }
  expect_silent(result_list <- results$get_list())
  expect_named(result_list,
               c("population", "total_pop", "occupancy", "total_occup",
                 "area"))
  expect_named(result_list$population, as.character(seq(0, 10, 2)))
  expect_true(all(unlist(lapply(result_list$population, length)) ==
                    region$get_locations()))
  expect_named(result_list$total_pop, as.character(0:10))
  expect_true(all(unlist(lapply(result_list$total_pop, length)) == 1))
  expect_named(result_list$occupancy, as.character(seq(0, 10, 2)))
  expect_true(all(unlist(lapply(result_list$occupancy, length)) ==
                    region$get_locations()))
  expect_named(result_list$total_occup, as.character(0:10))
  expect_true(all(unlist(lapply(result_list$total_occup, length)) == 1))
  expect_equal(lapply(result_list$population, function(rl) rl[1:12]), expected)
  expect_equal(unname(unlist(result_list$total_pop)), 1:11)
  expect_equal(lapply(result_list$occupancy, function(rl) rl[1:12]), expected)
  expect_equal(unname(unlist(result_list$total_occup)), 1:11)
  # staged population
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population <- StagedPopulation(region, growth = stage_matrix)
  expect_silent(results <- Results(region, population_model = population,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 1,
                                   combine_stages = NULL))
  n <- rep(0, region$get_locations())
  n[5922] <- 20
  set.seed(1234)
  n <- population$make(initial = n)
  n_t0 <- n
  for (tm in 0:10) {
    n[5922 + tm, ] <- n[5922,]*(tm == 0) + tm
    results$collate(r = 1, tm, n)
  }
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("population", "total_pop", "occupancy",
                              "total_occup", "area"))
  expect_named(result_list$population, as.character(seq(0, 10, 2)))
  expect_true(all(sapply(result_list$population, dim) ==
                    matrix(c(region$get_locations(), 3), nrow = 2, ncol = 6)))
  expect_named(result_list$total_pop, as.character(0:10))
  expect_true(all(sapply(result_list$total_pop, dim) ==
                    matrix(c(1, 3), nrow = 2, ncol = 11)))
  expect_named(result_list$occupancy, as.character(seq(0, 10, 2)))
  expect_true(all(unlist(lapply(result_list$occupancy, length)) ==
                    region$get_locations()))
  expect_named(result_list$total_occup, as.character(0:10))
  expect_true(all(unlist(lapply(result_list$total_occup, length)) == 1))
  n_total <- list('0' = n_t0[5922,, drop = FALSE])
  for (i in 1:10) {
    n_total[[as.character(i)]] <-  n_total[[as.character(i - 1)]] + i
  }
  expect_equal(result_list$total_pop, n_total)
  expect_equal(lapply(result_list$occupancy, function(rl) rl[5922:5933]),
               expected)
  unname(unlist(result_list$total_occup)) ; 1:11
  expect_silent(results <- Results(region, population_model = population,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 1,
                                   combine_stages = 2:3))
  n <- n_t0
  for (tm in 0:10) {
    n[5922 + tm, ] <- n[5922,]*(tm == 0) + tm
    results$collate(r = 1, tm, n)
  }
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("population", "total_pop", "occupancy",
                              "total_occup", "area"))
  expect_named(result_list$population, as.character(seq(0, 10, 2)))
  expect_true(all(sapply(result_list$population, dim) ==
                    matrix(c(region$get_locations(), 1), nrow = 2, ncol = 6)))
  expect_named(result_list$total_pop, as.character(0:10))
  expect_true(all(sapply(result_list$total_pop, dim) ==
                    matrix(c(1, 1), nrow = 2, ncol = 11)))
  n_total <- list('0' = as.matrix(sum(n_t0[5922, 2:3])))
  colnames(n_total[[1]]) <- "combined"
  for (i in 1:10) {
    n_total[[as.character(i)]] <-  n_total[[as.character(i - 1)]] + 2*i
  }
  expect_equal(result_list$total_pop, n_total)
})

test_that("collates and finalises multiple replicate results", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- UnstructPopulation(region)
  expect_silent(results <- Results(region, population_model = population,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 5,
                                   combine_stages = NULL))
  expect_silent(result_list <- results$get_list())
  expect_named(result_list,
               c("population", "total_pop", "occupancy", "total_occup",
                 "area"))
  expect_equal(unname(unlist(lapply(result_list$population, names))),
               rep(c("mean", "sd"), 6))
  expect_equal(unname(unlist(lapply(result_list$total_pop, names))),
               rep(c("mean", "sd"), 11))
  expect_equal(unname(unlist(lapply(result_list$area, names))),
               rep(c("mean", "sd"), 11))
  expect_equal(unname(unlist(lapply(result_list$occupancy, names))),
               rep("mean", 6))
  expect_equal(unname(unlist(lapply(result_list$total_occup, names))),
               rep(c("mean", "sd"), 11))
  expect_true(all(unlist(lapply(
    result_list$population,
    function(rl) lapply(rl, length))) == region$get_locations()))
  expect_true(all(unlist(result_list$population) == 0))
  expect_length(unlist(result_list$total_pop), 11*2)
  expect_true(all(unlist(result_list$total_pop) == 0))
  expect_length(unlist(result_list$area), 11*2)
  expect_true(all(unlist(result_list$area) == 0))
  expect_true(all(unlist(lapply(
    result_list$occupancy,
    function(rl) lapply(rl, length))) == region$get_locations()))
  expect_length(unlist(result_list$total_occup), 11*2)
  expected <- array(0, c(12, 5))
  for (r in 1:5) {
    n <- rep(0, region$get_locations())
    for (tm in 0:10) {
      n[tm + 1] <- tm %% r + (r %% 2)
      results$collate(r, tm, n)
    }
    expected[, r] <- n[1:12]
  }
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expected_mean <- rowMeans(expected)
  expected_sd <- apply(expected, 1, stats::sd)
  expected_mask <-
    lapply(seq(0, 10, 2), function(tm) c(rep(1, tm + 1), rep(0, 11 - tm)))
  names(expected_mask) <- as.character(seq(0, 10, 2))
  expect_equal(lapply(result_list$population, function(rl) rl$mean[1:12]),
               lapply(expected_mask, function(e) e*expected_mean))
  expect_equal(lapply(result_list$population, function(rl) rl$sd[1:12]),
               lapply(expected_mask, function(e) e*expected_sd))
  expect_equal(unname(unlist(lapply(result_list$total_pop,
                                    function(rl) rl$mean))),
               sapply(1:11, function(a)
                 mean(colSums(expected[1:a,,drop = FALSE]))))
  expect_equal(unname(unlist(lapply(result_list$total_pop,
                                    function(rl) rl$sd))),
               sapply(1:11, function(a)
                 stats::sd(colSums(expected[1:a,,drop = FALSE]))))

  expect_equal(unname(signif(unlist(lapply(result_list$area,
                                           function(rl) rl$mean)), 5)),
               sapply(1:11, function(a) mean(colSums(
                 expected[1:a,,drop = FALSE] > 0)))*1000^2)
  expect_equal(unname(signif(unlist(lapply(result_list$area,
                                           function(rl) rl$sd)), 5)),
               signif(sapply(1:11, function(a) stats::sd(colSums(
                 expected[1:a,,drop = FALSE] > 0)))*1000^2, 5))
  expect_equal(attr(result_list$area, "units"), "square metres")
  expected_mean <- rowMeans(expected > 0)
  expect_equal(lapply(result_list$occupancy, function(rl) rl$mean[1:12]),
               lapply(expected_mask, function(e) e*expected_mean))
  expect_equal(unname(unlist(lapply(result_list$total_occup,
                                    function(rl) rl$mean))),
               sapply(1:11, function(a)
                 mean(colSums((expected > 0)[1:a,,drop = FALSE]))))
  expect_equal(unname(unlist(lapply(result_list$total_occup,
                                    function(rl) rl$sd))),
               sapply(1:11, function(a)
                 stats::sd(colSums((expected > 0)[1:a,,drop = FALSE]))))
  # staged population
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population <- StagedPopulation(region, growth = stage_matrix)
  expect_silent(results <- Results(region, population_model = population,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 5,
                                   combine_stages = NULL))
  for (r in 1:5) {
    n <- rep(0, region$get_locations())
    n[5922] <- 20
    set.seed(1234)
    n <- population$make(initial = n)
    for (tm in 0:10) {
      n[5922 + tm, ] <- n[5922,]*(tm == 0) + tm %% r + (r %% 2)
      results$collate(r, tm, n)
    }
  }
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("population", "total_pop", "occupancy",
                              "total_occup", "area"))
  expect_named(result_list$population, as.character(seq(0, 10, 2)))
  expect_true(all(sapply(result_list$population, function(x) dim(x$mean)) ==
                    matrix(c(region$get_locations(), 3), nrow = 2, ncol = 6)))
  expect_true(all(sapply(result_list$population, function(x) dim(x$sd)) ==
                    matrix(c(region$get_locations(), 3), nrow = 2, ncol = 6)))
  expect_named(result_list$total_pop, as.character(0:10))
  expect_true(all(sapply(result_list$total_pop, function(x) dim(x$mean)) ==
                    matrix(c(1, 3), nrow = 2, ncol = 11)))
  expect_true(all(sapply(result_list$total_pop, function(x) dim(x$sd)) ==
                    matrix(c(1, 3), nrow = 2, ncol = 11)))
  expect_named(result_list$occupancy, as.character(seq(0, 10, 2)))
  expect_true(all(unlist(lapply(result_list$occupancy,
                                function(x) names(x))) == "mean"))
  expect_true(all(unlist(lapply(result_list$occupancy,
                                function(x) length(x$mean))) ==
                    region$get_locations()))
  expect_named(result_list$total_occup, as.character(0:10))
  expect_true(all(unlist(lapply(result_list$total_occup,
                                function(x) length(x$mean))) == 1))
  expect_true(all(unlist(lapply(result_list$total_occup,
                                function(x) length(x$sd))) == 1))
  expect_silent(results <- Results(region, population_model = population,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 5,
                                   combine_stages = 2:3))
  for (r in 1:5) {
    n <- rep(0, region$get_locations())
    n[5922] <- 20
    set.seed(1234)
    n <- population$make(initial = n)
    for (tm in 0:10) {
      n[5922 + tm, ] <- n[5922,]*(tm == 0) + tm %% r + (r %% 2)
      results$collate(r, tm, n)
    }
  }
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("population", "total_pop", "occupancy",
                              "total_occup", "area"))
  expect_named(result_list$population, as.character(seq(0, 10, 2)))
  expect_true(all(sapply(result_list$population,
                         function(x) colnames(x$mean)) == "combined"))
  expect_true(all(sapply(result_list$population, function(x) dim(x$mean)) ==
                    matrix(c(region$get_locations(), 1), nrow = 2, ncol = 6)))
  expect_true(all(sapply(result_list$population, function(x) dim(x$sd)) ==
                    matrix(c(region$get_locations(), 1), nrow = 2, ncol = 6)))
  expect_named(result_list$total_pop, as.character(0:10))
  expect_true(all(sapply(result_list$total_pop,
                         function(x) colnames(x$mean)) == "combined"))
  expect_true(all(sapply(result_list$total_pop, function(x) dim(x$mean)) ==
                    matrix(c(1, 1), nrow = 2, ncol = 11)))
  expect_true(all(sapply(result_list$total_pop, function(x) dim(x$sd)) ==
                    matrix(c(1, 1), nrow = 2, ncol = 11)))
})

test_that("collates spatially implicit area results", {
  region <- Region()
  region$set_max_implicit_area(1e8)
  population_model <- UnstructPopulation(region,
                                         growth = 2,
                                         capacity = 100,
                                         capacity_area = 1e6)
  expect_silent(results <- Results(region,
                                   population_model = population_model,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 1,
                                   combine_stages = NULL))
  n <- 10
  for (tm in 0:10) {
    attr(n, "diffusion_radius") <- 2000*tm
    results$collate(r = 1, tm, n)
  }
  result_list <- results$get_list()
  expect_named(result_list, c("population", "occupancy", "area"))
  expect_equal(unname(unlist(result_list$area)), pi*((0:10)*2000)^2)
  expect_equal(attr(result_list$area, "units"), "square metres")
  expect_silent(results <- Results(region,
                                   population_model = population_model,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 1,
                                   combine_stages = NULL))
  n <- 100
  for (tm in 0:10) {
    attr(n, "spread_area") <- as.numeric(n)/100*1e6
    results$collate(r = 1, tm, n)
    n <- n*population_model$get_growth()
  }
  result_list <- results$get_list()
  expect_equal(unname(unlist(result_list$area)),
               as.numeric(result_list$population)*1e6/100)
  expect_equal(attr(result_list$area, "units"), "square metres")
})
