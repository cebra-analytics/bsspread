context("Results")

test_that("initializes with region, population model, and other parameters", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  expect_error(results <- Results(region, population_model = 0),
               paste("Population model must be a 'Population' or inherited",
                     "class object."))
  expect_error(results <- Results(region, population_model = population,
                                  impacts = as.list(1:2)),
               paste("Impacts must be a list of 'Impacts' or inherited",
                     "class objects."))
  expect_error(results <- Results(region, population_model = population,
                                  actions = as.list(1:2)),
               paste("Actions must be a list of 'Actions' or inherited",
                     "class objects."))
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
  expect_named(result_list, c("occupancy", "total_occup", "area"))
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
  expect_equal(unname(unlist(result_list$total_occup)), 1:11)
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
  population <- PresencePopulation(region)
  expect_silent(results <- Results(region, population_model = population,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 5,
                                   combine_stages = NULL))
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("occupancy", "total_occup", "area"))
  expect_equal(unname(unlist(lapply(result_list$area, names))),
               rep(c("mean", "sd"), 11))
  expect_equal(unname(unlist(lapply(result_list$occupancy, names))),
               rep("mean", 6))
  expect_equal(unname(unlist(lapply(result_list$total_occup, names))),
               rep(c("mean", "sd"), 11))
  expect_equal(length(unlist(result_list$area)), 11*2)
  expect_true(all(unlist(result_list$area) == 0))
  expect_true(all(unlist(lapply(result_list$occupancy,
                                function(rl) lapply(rl, length))) ==
                    region$get_locations()))
  expect_equal(length(unlist(result_list$total_occup)), 11*2)
  expected <- array(FALSE, c(12, 5))
  for (r in 1:5) {
    n <- rep(FALSE, region$get_locations())
    for (tm in 0:10) {
      n[tm + 1] <- ((tm %% r) > 0) | ((r %% 2) > 0)
      results$collate(r, tm, n)
    }
    expected[, r] <- n[1:12]
  }
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expected_mean <- rowMeans(expected)
  expected_sd <- apply(expected, 1, stats::sd)
  expected_mask <- lapply(seq(0, 10, 2),
                          function(tm) c(rep(1, tm + 1), rep(0, 11 - tm)))
  names(expected_mask) <- as.character(seq(0, 10, 2))
  expect_equal(unname(signif(unlist(lapply(result_list$area,
                                           function(rl) rl$mean)), 5)),
               sapply(1:11, function(a) mean(
                 colSums(expected[1:a,,drop = FALSE] > 0)))*1000^2)
  expect_equal(unname(signif(unlist(lapply(result_list$area,
                                           function(rl) rl$sd)), 5)),
               signif(sapply(1:11, function(a) stats::sd(
                 colSums(expected[1:a,,drop = FALSE] > 0)))*1000^2, 5))
  expect_equal(attr(result_list$area, "units"), "square metres")
  expected_mean <- rowMeans(expected > 0)
  expect_equal(lapply(result_list$occupancy, function(rl) rl$mean[1:12]),
               lapply(expected_mask, function(e) e*expected_mean))
  expect_equal(unname(unlist(lapply(result_list$total_occup,
                                    function(rl) rl$mean))),
               sapply(1:11, function(a) mean(
                 colSums((expected > 0)[1:a,,drop = FALSE]))))
  expect_equal(unname(unlist(lapply(result_list$total_occup,
                                    function(rl) rl$sd))),
               sapply(1:11, function(a) stats::sd(
                 colSums((expected > 0)[1:a,,drop = FALSE]))))
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
  population_model <- PresencePopulation(region)
  expect_silent(results <- Results(region, population_model = population_model,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 1,
                                   combine_stages = NULL))
  n <- TRUE
  for (tm in 0:10) {
    attr(n, "diffusion_radius") <- 2000*tm
    results$collate(r = 1, tm, n)
  }
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("occupancy", "area"))
  expect_named(result_list$occupancy, as.character(seq(0, 10, 1)))
  expect_equal(unname(unlist(result_list$occupancy)), rep(1, 11))
  expect_named(result_list$area, as.character(seq(0, 10, 1)))
  expect_equal(unname(unlist(result_list$area)), pi*((0:10)*2000)^2)
  expect_equal(attr(result_list$area, "units"), "square metres")
  expect_silent(results <- Results(region, population_model = population_model,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 5,
                                   combine_stages = NULL))
  for (r in 1:5) {
    n <- (r < 5)
    for (tm in 0:10) {
      attr(n, "diffusion_radius") <- 2000*tm + 100*r
      results$collate(r, tm, n)
    }
  }
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("occupancy", "area"))
  expect_named(result_list$occupancy, as.character(seq(0, 10, 1)))
  expect_named(unlist(unname(result_list$occupancy)), rep("mean", 11))
  expect_equal(unname(unlist(unname(result_list$occupancy))), rep(0.8, 11))
  expect_named(result_list$area, as.character(seq(0, 10, 1)))
  expect_named(unlist(unname(result_list$area)), rep(c("mean", "sd"), 11))

  population_model <- UnstructPopulation(region,
                                         growth = 2,
                                         capacity = 100,
                                         capacity_area = 1e6)
  expect_silent(results <- Results(region, population_model = population_model,
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
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("population", "occupancy", "area"))
  expect_equal(unname(unlist(result_list$area)),
               as.numeric(result_list$population)*1e6/100)
  expect_equal(attr(result_list$area, "units"), "square metres")
  expect_silent(results <- Results(region, population_model = population_model,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 5,
                                   combine_stages = NULL))
  for (r in 1:5) {
    n <- 100 + r
    for (tm in 0:10) {
      attr(n, "spread_area") <- as.numeric(n)/100*1e6
      results$collate(r, tm, n)
      n <- n*population_model$get_growth()
    }
  }
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("population", "occupancy", "area"))
  expect_named(unlist(unname(result_list$population)),
               rep(c("mean", "sd"), 11))
  expect_named(unlist(unname(result_list$occupancy)), rep("mean", 11))
  expect_named(unlist(unname(result_list$area)), rep(c("mean", "sd"), 11))
  # staged population
  region <- Region()
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, growth = stage_matrix)
  expect_silent(results <- Results(region, population_model = population_model,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 1,
                                   combine_stages = NULL))
  n <- 20
  set.seed(1234)
  n <- population_model$make(initial = n)
  for (tm in 0:10) {
    n <- n + tm
    results$collate(r = 1, tm, n)
  }
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("population", "occupancy", "area"))
  expect_equal(unname(unlist(lapply(result_list$population, colnames))),
               rep(c("stage 1", "stage 2", "stage 3"), 11))
  expect_equal(unname(unlist(result_list$occupancy)), rep(1, 11))
  expect_equal(unname(unlist(result_list$area)), rep(1, 11))
  expect_silent(results <- Results(region, population_model = population_model,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 1,
                                   combine_stages = 2:3))
  n <- 20
  set.seed(1234)
  n <- population_model$make(initial = n)
  for (tm in 0:10) {
    n <- n + tm
    results$collate(r = 1, tm, n)
  }
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("population", "occupancy", "area"))
  expect_equal(unname(unlist(lapply(result_list$population, colnames))),
               rep("combined", 11))
  expect_equal(unname(unlist(result_list$occupancy)), rep(1, 11))
  expect_equal(unname(unlist(result_list$area)), rep(1, 11))
  expect_silent(results <- Results(region, population_model = population_model,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 5,
                                   combine_stages = NULL))
  for (r in 1:5) {
    n <- 20
    set.seed(1234)
    n <- population_model$make(initial = n)
    for (tm in 0:10) {
      n <- n + r + tm
      results$collate(r, tm, n)
    }
  }
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("population", "occupancy", "area"))
  expect_equal(unname(unlist(lapply(result_list$population, names))),
               rep(c("mean", "sd"), 11))
  expect_equal(unname(unlist(lapply(result_list$population,
                                    function(x) lapply(x, colnames)))),
               rep(c("stage 1", "stage 2", "stage 3"), 22))
  expect_named(unlist(unname(result_list$occupancy)), rep("mean", 11))
  expect_named(unlist(unname(result_list$area)), rep(c("mean", "sd"), 11))
  expect_silent(results <- Results(region, population_model = population_model,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 5,
                                   combine_stages = 2:3))
  for (r in 1:5) {
    n <- 20
    set.seed(1234)
    n <- population_model$make(initial = n)
    for (tm in 0:10) {
      n <- n + r + tm
      results$collate(r, tm, n)
    }
  }
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("population", "occupancy", "area"))
  expect_equal(unname(unlist(lapply(result_list$population, names))),
               rep(c("mean", "sd"), 11))
  expect_equal(unname(unlist(lapply(result_list$population,
                                    function(x) lapply(x, colnames)))),
               rep("combined", 22))
  expect_named(unlist(unname(result_list$occupancy)), rep("mean", 11))
  expect_named(unlist(unname(result_list$area)), rep(c("mean", "sd"), 11))
})

test_that("initializes results with impacts", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  pops <- region$get_locations()
  population_model <- UnstructPopulation(region, growth = 1.2)
  asset_value_1 <- 100*(template > 0.1 & template < 0.3)
  asset_value_2 <- 200*(template > 0.2 & template < 0.4)
  asset_value_3 <- 300*(template > 0.1 & template < 0.3)
  asset_value_4 <- 400*(template > 0.2 & template < 0.4)
  impacts_1 <- Impacts(region, population_model, # monetary
                       asset_name = "impact1",
                       asset_value = asset_value_1,
                       loss_rate = 0.3)
  impacts_1$set_id(1)
  impacts_2 <- Impacts(region, population_model, # monetary
                       asset_name = "impact2",
                       asset_value = asset_value_2,
                       loss_rate = 0.4)
  impacts_2$set_id(2)
  impacts_3 <- Impacts(region, population_model,
                       valuation_type = "non-monetary",
                       asset_name = "impact3",
                       asset_value = asset_value_3,
                       loss_rate = 0.3)
  impacts_3$set_id(3)
  impacts_4 <- Impacts(region, population_model,
                       valuation_type = "non-monetary",
                       asset_name = "impact4",
                       asset_value = asset_value_4,
                       loss_rate = 0.4)
  impacts_4$set_id(4)
  # single replicate
  expect_silent(results <- Results(region, population_model = population_model,
                     impacts = list(impacts_1, impacts_2, impacts_3,
                                    impacts_4),
                     time_steps = 10, collation_steps = 2,
                     replicates = 1))
  expect_silent(result_list <- results$get_list())
  expect_named(result_list, c("population", "total_pop", "occupancy",
                              "total_occup", "area", "impacts"))
  expect_named(result_list$impacts, c("idx", "monetary", "non-monetary"))
  expect_equal(result_list$impacts$idx,
               list(monetary = c("impact1" = 1, "impact2" = 2),
                    `non-monetary` = c("impact3" = 3, "impact4" = 4)))
  expect_named(result_list$impacts$monetary,
               c("total", "impact1", "impact2", "combined", "cumulative"))
  expect_named(result_list$impacts[["non-monetary"]],
               c("total", "impact3", "impact4"))
  collated <- lapply(1:6, function(i) rep(0, pops))
  names(collated) <- as.character(seq(0, 10, 2))
  expect_equal(result_list$impacts[["non-monetary"]][c("impact3", "impact4")],
               list(impact3 = collated, impact4 = collated))
  attr(collated, "unit") <- "$"
  expect_equal(result_list$impacts$monetary[c("impact1", "impact2",
                                              "combined")],
               list(impact1 = collated, impact2 = collated,
                    combined = collated))
  totals <- lapply(1:11, function(i) 0)
  names(totals) <- as.character(0:10)
  expect_equal(result_list$impacts[["non-monetary"]]$total,
               list(impact3 = totals, impact4 = totals))
  attr(totals, "unit") <- "$"
  expect_equal(result_list$impacts$monetary$total,
               list(impact1 = totals, impact2 = totals, combined = totals))
  expect_named(result_list$impacts$monetary$cumulative,
               c("total", "impact1", "impact2", "combined"))
  expect_equal(result_list$impacts$monetary$cumulative[c("impact1", "impact2",
                                                         "combined")],
               list(impact1 = collated, impact2 = collated,
                    combined = collated))
  expect_equal(result_list$impacts$monetary$cumulative$total,
               list(impact1 = totals, impact2 = totals, combined = totals))
  # multiple replicates
  expect_silent(results <- Results(region, population_model = population_model,
                                   impacts = list(impacts_1, impacts_2,
                                                  impacts_3, impacts_4),
                                   time_steps = 4, collation_steps = 2,
                                   replicates = 3))
  expect_silent(result_list <- results$get_list())
  collated <- lapply(1:3, function(i)
    list(mean = rep(0, pops), sd = rep(0, pops)))
  names(collated) <- as.character(seq(0, 4, 2))
  expect_equal(result_list$impacts[["non-monetary"]][c("impact3", "impact4")],
               list(impact3 = collated, impact4 = collated))
  attr(collated, "unit") <- "$"
  expect_equal(result_list$impacts$monetary[c("impact1", "impact2",
                                              "combined")],
               list(impact1 = collated, impact2 = collated,
                    combined = collated))
  totals <- lapply(1:5, function(i) list(mean = 0, sd = 0))
  names(totals) <- as.character(0:4)
  expect_equal(result_list$impacts[["non-monetary"]]$total,
               list(impact3 = totals, impact4 = totals))
  attr(totals, "unit") <- "$"
  expect_equal(result_list$impacts$monetary$total,
               list(impact1 = totals, impact2 = totals, combined = totals))
  expect_equal(result_list$impacts$monetary$cumulative[c("impact1", "impact2",
                                                         "combined")],
               list(impact1 = collated, impact2 = collated,
                    combined = collated))
  expect_equal(result_list$impacts$monetary$cumulative$total,
               list(impact1 = totals, impact2 = totals, combined = totals))
})

test_that("collates and finalizes impact results", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  idx <- 5918:5922
  region <- Region(template)
  template[region$get_indices()][idx[4:5],] <- 0.25
  population_model <-
    UnstructPopulation(region, growth = 1.2,
                       capacity = rep(10, region$get_locations()))
  n <- rep(0, region$get_locations())
  n[idx] <- 7:11
  asset_value_1 <- 100*(template > 0.1 & template < 0.3)
  asset_value_2 <- 200*(template > 0.2 & template < 0.4)
  asset_value_3 <- 300*(template > 0.1 & template < 0.3)
  asset_value_4 <- 400*(template > 0.2 & template < 0.4)
  impacts_1 <- Impacts(region, population_model, # monetary
                       impact_type = "density",
                       asset_name = "impact1",
                       asset_value = asset_value_1,
                       loss_rate = 0.3)
  impacts_1$set_id(1)
  impacts_2 <- Impacts(region, population_model, # monetary
                       impact_type = "density",
                       asset_name = "impact2",
                       asset_value = asset_value_2,
                       loss_rate = 0.4)
  impacts_2$set_id(2)
  impacts_3 <- Impacts(region, population_model,
                       impact_type = "density",
                       valuation_type = "non-monetary",
                       asset_name = "impact3",
                       asset_value = asset_value_3,
                       loss_rate = 0.3)
  impacts_3$set_id(3)
  impacts_4 <- Impacts(region, population_model,
                       impact_type = "density",
                       valuation_type = "non-monetary",
                       asset_name = "impact4",
                       asset_value = asset_value_4,
                       loss_rate = 0.4)
  impacts_4$set_id(4)
  impacts <- list(impacts_1, impacts_2, impacts_3, impacts_4)
  calc_impacts <- lapply(impacts, function(impacts_i)
    attr(impacts_i$calculate(n, 0), "impacts")[[impacts_i$get_id()]])
  attr(n, "impacts") <- calc_impacts
  # single replicate
  expect_silent(results <- Results(region, population_model = population_model,
                                   impacts = impacts, time_steps = 4,
                                   collation_steps = 2, replicates = 1))
  expect_silent(results$collate(r = 1, tm = 0, n = n))
  for (tm in 1:4) {
    results$collate(r = 1, tm = tm, n = n)
  }
  expect_silent(results$finalize())
  result_list <- results$get_list()
  calc_impacts$combined <- calc_impacts[[1]] + calc_impacts[[2]]
  collated <- lapply(1:5, function(i) {
    collated_tm <- lapply(1:3, function(tm) calc_impacts[[i]])
    names(collated_tm) <- as.character(seq(0, 4, 2))
    collated_tm
  })
  totals <- lapply(1:5, function(i) {
    totals_tm <- lapply(1:5, function(tm) sum(calc_impacts[[i]]))
    names(totals_tm) <- as.character(seq(0, 4, 1))
    totals_tm
  })
  mult <- c(1, 3, 5)
  cumulative <- lapply(1:5, function(i) {
    cumul_tm <- lapply(1:3, function(tm) collated[[i]][[tm]]*mult[[tm]])
    names(cumul_tm) <- as.character(seq(0, 4, 2))
    cumul_tm
  })
  cumul_totals <- lapply(1:5, function(i) {
    totals_tm <- lapply(1:5, function(tm) totals[[i]][[tm]]*tm)
    names(totals_tm) <- as.character(seq(0, 4, 1))
    totals_tm
  })
  for (i in c(1:2, 5)) {
    attr(collated[[i]], "unit") <- "$"
    attr(totals[[i]], "unit") <- "$"
    attr(cumulative[[i]], "unit") <- "$"
    attr(cumul_totals[[i]], "unit") <- "$"
  }
  expect_equal(result_list$impacts$monetary[c("impact1", "impact2",
                                              "combined")],
               list(impact1 = collated[[1]], impact2 = collated[[2]],
                    combined = collated[[5]]))
  expect_equal(result_list$impacts[["non-monetary"]][c("impact3", "impact4")],
               list(impact3 = collated[[3]], impact4 = collated[[4]]))
  expect_equal(result_list$impacts$monetary$total,
               list(impact1 = totals[[1]], impact2 = totals[[2]],
                    combined = totals[[5]]))
  expect_equal(result_list$impacts[["non-monetary"]]$total,
               list(impact3 = totals[[3]], impact4 = totals[[4]]))
  expect_equal(result_list$impacts$monetary$cumulative[c("impact1", "impact2",
                                                         "combined")],
               list(impact1 = cumulative[[1]], impact2 = cumulative[[2]],
                    combined = cumulative[[5]]))
  expect_equal(result_list$impacts$monetary$cumulative$total,
               list(impact1 = cumul_totals[[1]], impact2 = cumul_totals[[2]],
                    combined = cumul_totals[[5]]))
  # multiple replicates
  expect_silent(results <- Results(region, population_model = population_model,
                                   impacts = impacts, time_steps = 4,
                                   collation_steps = 2, replicates = 3))
  attr(n, "impacts") <- NULL
  n_r <- lapply(-1:1, function(i) {
    n[idx] <- n[idx] + i
    calc_impacts <- lapply(impacts, function(impacts_i)
      attr(impacts_i$calculate(n, 0), "impacts")[[impacts_i$get_id()]])
    attr(n, "impacts") <- calc_impacts
    n
  })
  for (i in 1:3) {
    for (tm in 0:4) {
      results$collate(r = i, tm = tm, n = n_r[[i]])
    }
  }
  expect_silent(results$finalize())
  result_list <- results$get_list()
  calc_impacts <- lapply(n_r, function(n) {
    calc_impacts_i <- attr(n, "impacts")
    calc_impacts_i$combined <- calc_impacts_i[[1]] + calc_impacts_i[[2]]
    calc_impacts_i
  })
  calc_impacts <- lapply(1:5, function(i) {
    sapply(calc_impacts, function(calc_impacts_j) calc_impacts_j[[i]])
  })
  collated <- lapply(1:5, function(i) {
    collated_tm <- lapply(1:3, function(tm) {
      list(mean = rowMeans(calc_impacts[[i]]),
           sd = apply(calc_impacts[[i]], 1, sd))
    })
    names(collated_tm) <- as.character(seq(0, 4, 2))
    collated_tm
  })
  totals <- lapply(1:5, function(i) {
    totals_tm <- lapply(1:5, function(tm) {
      total_i_tm <- colSums(calc_impacts[[i]])
      list(mean = mean(total_i_tm), sd = sd(total_i_tm))
    })
    names(totals_tm) <- as.character(seq(0, 4, 1))
    totals_tm
  })
  cumulative <- lapply(1:5, function(i) {
    cumul_tm <- lapply(1:3, function(tm) {
      cumul_i_tm <- calc_impacts[[i]]*mult[[tm]]
      list(mean = rowMeans(cumul_i_tm), sd = apply(cumul_i_tm, 1, sd))
    })
    names(cumul_tm) <- as.character(seq(0, 4, 2))
    cumul_tm
  })
  cumul_totals <- lapply(1:5, function(i) {
    totals_tm <- lapply(1:5, function(tm) {
      cumul_i_tm <- colSums(calc_impacts[[i]])*tm
      list(mean = mean(cumul_i_tm), sd = sd(cumul_i_tm))
    })
    names(totals_tm) <- as.character(seq(0, 4, 1))
    totals_tm
  })
  for (i in c(1:2, 5)) {
    attr(collated[[i]], "unit") <- "$"
    attr(totals[[i]], "unit") <- "$"
    attr(cumulative[[i]], "unit") <- "$"
    attr(cumul_totals[[i]], "unit") <- "$"
  }
  expect_equal(result_list$impacts$monetary[c("impact1", "impact2",
                                              "combined")],
               list(impact1 = collated[[1]], impact2 = collated[[2]],
                    combined = collated[[5]]))
  expect_equal(result_list$impacts[["non-monetary"]][c("impact3", "impact4")],
               list(impact3 = collated[[3]], impact4 = collated[[4]]))
  expect_equal(result_list$impacts$monetary$total,
               list(impact1 = totals[[1]], impact2 = totals[[2]],
                    combined = totals[[5]]))
  expect_equal(result_list$impacts[["non-monetary"]]$total,
               list(impact3 = totals[[3]], impact4 = totals[[4]]))
  expect_equal(result_list$impacts$monetary$cumulative[c("impact1", "impact2",
                                                         "combined")],
               list(impact1 = cumulative[[1]], impact2 = cumulative[[2]],
                    combined = cumulative[[5]]))
  expect_equal(result_list$impacts$monetary$cumulative$total,
               list(impact1 = cumul_totals[[1]], impact2 = cumul_totals[[2]],
                    combined = cumul_totals[[5]]))
})

test_that("collates and finalizes spatially implicit impact results", {
  region <- Region()
  population_model <- UnstructPopulation(region, growth = 1.2)
  impacts_1 <- Impacts(region, population_model, # monetary
                       impact_type = "area",
                       asset_name = "impact1",
                       asset_value = 0.01,
                       loss_rate = 0.3)
  impacts_1$set_id(1)
  impacts_2 <- Impacts(region, population_model, # monetary
                       impact_type = "area",
                       asset_name = "impact2",
                       asset_value = 0.02,
                       loss_rate = 0.4)
  impacts_2$set_id(2)
  impacts_3 <- Impacts(region, population_model,
                       impact_type = "area",
                       valuation_type = "non-monetary",
                       asset_name = "impact3",
                       asset_value = 0.03,
                       loss_rate = 0.3)
  impacts_3$set_id(3)
  impacts_4 <- Impacts(region, population_model,
                       impact_type = "area",
                       valuation_type = "non-monetary",
                       asset_name = "impact4",
                       asset_value = 0.04,
                       loss_rate = 0.4)
  impacts_4$set_id(4)
  impacts <- list(impacts_1, impacts_2, impacts_3, impacts_4)
  n <- 8
  attr(n, "spread_area") <- 10000
  calc_impacts <-
    lapply(impacts, function(impacts_i) attr(impacts_i$calculate(n, 0),
                                             "impacts")[[impacts_i$get_id()]])
  attr(n, "impacts") <- calc_impacts
  # single replicate
  expect_silent(results <- Results(region, population_model = population_model,
                                   impacts = impacts, time_steps = 4,
                                   collation_steps = 2, replicates = 1))
  expect_silent(results$collate(r = 1, tm = 0, n = n))
  for (tm in 1:4) {
    results$collate(r = 1, tm = tm, n = n)
  }
  expect_silent(results$finalize())
  result_list <- results$get_list()
  calc_impacts$combined <- calc_impacts[[1]] + calc_impacts[[2]]
  collated <- lapply(1:5, function(i) {
    collated_tm <- lapply(1:5, function(tm) calc_impacts[[i]])
    names(collated_tm) <- as.character(seq(0, 4, 1))
    collated_tm
  })
  cumulative <- lapply(1:5, function(i) {
    cumul_tm <- lapply(1:5, function(tm) calc_impacts[[i]]*tm)
    names(cumul_tm) <- as.character(seq(0, 4, 1))
    cumul_tm
  })
  for (i in c(1:2, 5)) {
    attr(collated[[i]], "unit") <- "$"
    attr(cumulative[[i]], "unit") <- "$"
  }
  expect_named(result_list$impacts, c("idx", "monetary", "non-monetary"))
  expect_equal(result_list$impacts$idx,
               list(monetary = c("impact1" = 1, "impact2" = 2),
                    `non-monetary` = c("impact3" = 3, "impact4" = 4)))
  expect_named(result_list$impacts$monetary,
               c("impact1", "impact2", "combined", "cumulative"))
  expect_named(result_list$impacts[["non-monetary"]], c("impact3", "impact4"))
  expect_named(result_list$impacts$monetary$cumulative,
               c("impact1", "impact2", "combined"))
  expect_equal(result_list$impacts$monetary[c("impact1", "impact2",
                                              "combined")],
               list(impact1 = collated[[1]], impact2 = collated[[2]],
                    combined = collated[[5]]))
  expect_equal(result_list$impacts[["non-monetary"]][c("impact3", "impact4")],
               list(impact3 = collated[[3]], impact4 = collated[[4]]))
  expect_equal(result_list$impacts$monetary$cumulative[c("impact1", "impact2",
                                                         "combined")],
               list(impact1 = cumulative[[1]], impact2 = cumulative[[2]],
                    combined = cumulative[[5]]))
  # multiple replicates
  expect_silent(results <- Results(region, population_model = population_model,
                     impacts = impacts, time_steps = 4,
                     collation_steps = 2, replicates = 3))
  attr(n, "impacts") <- NULL
  n_r <- lapply(-1:1, function(i) {
    n <- n + 2*i
    attr(n ,"spread_area") <- attr(n,"spread_area") + 2000*i
    calc_impacts <- lapply(impacts, function(impacts_i)
      attr(impacts_i$calculate(n, 0), "impacts")[[impacts_i$get_id()]])
    attr(n, "impacts") <- calc_impacts
    n
  })
  for (i in 1:3) {
    for (tm in 0:4) {
      results$collate(r = i, tm = tm, n = n_r[[i]])
    }
  }
  expect_silent(results$finalize())
  result_list <- results$get_list()
  calc_impacts <- lapply(n_r, function(n) {
    calc_impacts_i <- attr(n, "impacts")
    calc_impacts_i$combined <- calc_impacts_i[[1]] + calc_impacts_i[[2]]
    calc_impacts_i
  })
  calc_impacts <- lapply(1:5, function(i) {
    sapply(calc_impacts, function(calc_impacts_j) calc_impacts_j[[i]])
  })
  collated <- lapply(1:5, function(i) {
    collated_tm <- lapply(1:5, function(tm) {
      list(mean = mean(calc_impacts[[i]]),
           sd = sd(calc_impacts[[i]]))
    })
    names(collated_tm) <- as.character(seq(0, 4, 1))
    collated_tm
  })
  cumulative <- lapply(1:5, function(i) {
    cumul_tm <- lapply(1:5, function(tm) {
      cumul_i_tm <- calc_impacts[[i]]*tm
      list(mean = mean(cumul_i_tm), sd = sd(cumul_i_tm))
    })
    names(cumul_tm) <- as.character(seq(0, 4, 1))
    cumul_tm
  })
  for (i in c(1:2, 5)) {
    attr(collated[[i]], "unit") <- "$"
    attr(cumulative[[i]], "unit") <- "$"
  }
  expect_equal(result_list$impacts$monetary[c("impact1", "impact2",
                                              "combined")],
               list(impact1 = collated[[1]], impact2 = collated[[2]],
                    combined = collated[[5]]))
  expect_equal(result_list$impacts[["non-monetary"]][c("impact3", "impact4")],
               list(impact3 = collated[[3]], impact4 = collated[[4]]))
  expect_equal(result_list$impacts$monetary$cumulative[c("impact1", "impact2",
                                                         "combined")],
               list(impact1 = cumulative[[1]], impact2 = cumulative[[2]],
                    combined = cumulative[[5]]))
})

test_that("initializes results with actions (and monetary impacts)", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  sensitivity <- manage_pr <- template[region$get_indices()][,1]
  pops <- region$get_locations()
  population_model <- UnstructPopulation(region, growth = 1.2)
  actions_1 <- Detection(region, population_model,
                         sensitivity = sensitivity,
                         sensitivity_type = "presence")
  actions_1$set_id(1)
  actions_2 <- Controls(region, population_model,
                        control_type = "growth",
                        control_mult = 0.4)
  actions_2$set_id(2)
  actions_3 <- Controls(region, population_model,
                        control_type = "spread",
                        control_mult = 0.5)
  actions_3$set_id(3)
  actions_4 <- Removals(region, population_model,
                        removal_pr_type = "population")
  actions_4$set_id(4)
  # single replicate
  expect_silent(results <- Results(region,
                                   population_model = population_model,
                                   actions = list(actions_1, actions_2,
                                                  actions_3, actions_4),
                                   time_steps = 10, collation_steps = 2,
                                   replicates = 1))
  result_list <- results$get_list()
  expect_named(result_list, c("population", "total_pop", "occupancy",
                              "total_occup", "area", "actions"))
  collated <- lapply(1:6, function(i) rep(FALSE, pops))
  names(collated) <- as.character(seq(0, 10, 2))
  totals <- lapply(1:11, function(i) 0)
  names(totals) <- as.character(0:10)
  expect_equal(result_list$actions[[1]],
               list(detected = collated, total = totals))
  expect_equal(result_list$actions[[2]],
               list(control_growth = collated, total = totals))
  expect_equal(result_list$actions[[3]],
               list(control_spread = collated, total = totals))
  expect_equal(result_list$actions[[4]],
               list(removed = collated, total = totals))
  # multiple replicates
  expect_silent(results <- Results(region,
                                   population_model = population_model,
                                   actions = list(actions_1, actions_2,
                                                  actions_3, actions_4),
                                   time_steps = 4, collation_steps = 2,
                                   replicates = 3))
  result_list <- results$get_list()
  collated_no_sd <- lapply(1:3, function(i) list(mean = rep(0, pops)))
  collated_with_sd <- lapply(1:3, function(i)
    list(mean = rep(0, pops), sd = rep(0, pops)))
  names(collated_no_sd) <- as.character(seq(0, 4, 2))
  names(collated_with_sd) <- as.character(seq(0, 4, 2))
  totals <- lapply(1:5, function(i) list(mean = 0, sd = 0))
  names(totals) <- as.character(0:4)
  expect_equal(result_list$actions[[1]],
               list(detected = collated_no_sd, total = totals))
  expect_equal(result_list$actions[[2]],
               list(control_growth = collated_no_sd, total = totals))
  expect_equal(result_list$actions[[3]],
               list(control_spread = collated_no_sd, total = totals))
  expect_equal(result_list$actions[[4]],
               list(removed = collated_no_sd, total = totals))
  # include individual numbers
  actions_1 <- Detection(region, population_model,
                         sensitivity = sensitivity,
                         sensitivity_type = "individual")
  actions_1$set_id(1)
  actions_2 <- Controls(region, population_model,
                        control_type = "search_destroy",
                        manage_pr = manage_pr,
                        manage_pr_type = "individual")
  actions_2$set_id(2)
  actions_3 <- Removals(region, population_model,
                        removal_pr_type = "individual")
  actions_3$set_id(3)
  expect_silent(results <- Results(region,
                                   population_model = population_model,
                                   actions = list(actions_1, actions_2,
                                                  actions_3),
                                   time_steps = 4, collation_steps = 2,
                                   replicates = 3))
  result_list <- results$get_list()
  expect_equal(result_list$actions[[1]],
               list(detected = collated_no_sd, total = totals,
                    number = list(detected = collated_with_sd,
                                  total = totals)))
  expect_equal(result_list$actions[[2]],
               list(control_search_destroy = collated_no_sd, total = totals,
                    number = list(control_search_destroy = collated_with_sd,
                                  total = totals)))
  expect_equal(result_list$actions[[3]],
               list(removed = collated_no_sd, total = totals,
                    number = list(removed = collated_with_sd,
                                  total = totals)))
  # with action costs and monetary-only impacts
  asset_value_1 <- 100*(template > 0.1 & template < 0.3)
  asset_value_2 <- 200*(template > 0.2 & template < 0.4)
  impacts_1 <- Impacts(region, population_model, # monetary
                       asset_name = "impact1",
                       asset_value = asset_value_1,
                       loss_rate = 0.3)
  impacts_1$set_id(1)
  impacts_2 <- Impacts(region, population_model, # monetary
                       asset_name = "impact2",
                       asset_value = asset_value_2,
                       loss_rate = 0.4)
  impacts_2$set_id(2)
  surv_cost <- 3; attr(surv_cost, "unit") <- "$"
  actions_1 <- Detection(region, population_model,
                         sensitivity = sensitivity,
                         sensitivity_type = "presence",
                         surv_cost = surv_cost)
  actions_1$set_id(1)
  control_cost <- 4; attr(control_cost, "unit") <- "$"
  actions_2 <- Controls(region, population_model,
                        control_type = "growth",
                        control_mult = 0.4,
                        control_cost = control_cost)
  actions_2$set_id(2)
  removal_cost <- 6; attr(removal_cost, "unit") <- "$"
  actions_3 <- Removals(region, population_model,
                        removal_pr_type = "population",
                        removal_cost = removal_cost)
  actions_3$set_id(3)
  expect_silent(results <- Results(region,
                                   population_model = population_model,
                                   impacts = list(impacts_1, impacts_2),
                                   actions = list(actions_1, actions_2,
                                                  actions_3),
                                   time_steps = 4, collation_steps = 2,
                                   replicates = 3))
  result_list <- results$get_list()
  expected_costs <- list(detected = collated_with_sd, total = totals,
                         cumulative = list(detected = collated_with_sd,
                                           total = totals))
  attr(expected_costs, "unit") <- "$"
  expect_equal(result_list$actions[[1]],
               list(detected = collated_no_sd, total = totals,
                    cost = expected_costs))
  names(expected_costs)[1] <- "control_growth"
  names(expected_costs$cumulative)[1] <- "control_growth"
  expect_equal(result_list$actions[[2]],
               list(control_growth = collated_no_sd, total = totals,
                    cost = expected_costs))
  names(expected_costs)[1] <- "removed"
  names(expected_costs$cumulative)[1] <- "removed"
  expect_equal(result_list$actions[[3]],
               list(removed = collated_no_sd, total = totals,
                    cost = expected_costs))
  names(expected_costs)[1] <- "combined"
  names(expected_costs$cumulative)[1] <- "combined"
  expect_equal(result_list$actions$cost, expected_costs)
  expect_equal(result_list$cost, expected_costs)
})

test_that("collates and finalizes action results with costs", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  idx <- 5916:5922
  region <- bsspread::Region(template)
  template[region$get_indices()][idx,] <- c(rep(0.5, 4), 0.5, 0.75, 1)
  sensitivity <- removal_pr <- template[region$get_indices()][,1]
  population_model <- bsspread::UnstructPopulation(region, growth = 1.2)
  n <- rep(0, region$get_locations())
  idx <- idx[5:7]
  n_list <- list(n, n, n)
  n_list[[1]][idx] <- (10:12)*5
  n_list[[2]][idx] <- c(9, 0, 11)*5
  n_list[[3]][idx] <- c(0, 12, 13)*5
  asset_value_1 <- 10*(template < 0.9)
  asset_value_2 <- 20*(template > 0.6)
  impacts_1 <- Impacts(region, population_model, # monetary
                       asset_name = "impact1",
                       asset_value = asset_value_1,
                       loss_rate = 0.1)
  impacts_1$set_id(1)
  impacts_2 <- Impacts(region, population_model, # monetary
                       asset_name = "impact2",
                       asset_value = asset_value_2,
                       loss_rate = 0.2)
  impacts_2$set_id(2)
  impacts <- list(impacts_1, impacts_2)
  surv_cost <- 3; attr(surv_cost, "unit") <- "$"
  actions_1 <- Detection(region, population_model,
                         sensitivity = sensitivity,
                         surv_cost = surv_cost,
                         schedule = 2:3)
  actions_1$set_id(1)
  removal_cost <- 4; attr(removal_cost, "unit") <- "$"
  actions_2 <- Removals(region, population_model,
                        removal_pr = removal_pr,
                        removal_cost = removal_cost,
                        schedule = 2:3)
  actions_2$set_id(2)
  control_cost <- 5; attr(control_cost, "unit") <- "$"
  actions_3 <- Controls(region, population_model,
                        control_type = "growth",
                        control_mult = 0.7,
                        control_cost = control_cost,
                        schedule = 2:3)
  actions_3$set_id(3)
  actions <- list(actions_1, actions_2, actions_3)
  # single replicate
  expect_silent(
    results <- Results(region, population_model = population_model,
                       impacts = impacts, actions = actions,
                       time_steps = 4, collation_steps = 2,
                       replicates = 1))
  expect_silent(action_results_0 <- results$get_list()$actions)
  action_results_0$combined_cost_plus_impacts <- action_results_0$cost$combined
  action_results_0$combined_cum_cost_plus_impacts <-
    action_results_0$cost$combined
  action_results_0$total_cost_plus_impacts <- action_results_0$cost$total
  action_results_0$total_cum_cost_plus_impacts <- action_results_0$cost$total
  expect_list <- list(action_results_0, action_results_0, action_results_0)
  set.seed(1234)
  for (r in 1:3) { # collect expected for single and multiple replicates
    n <- n_list[[r]]
    for (tm in 0:4) {
      tmc <- as.character(tm)
      tmc_prev <- as.character(max(tm - 1, 0))
      for (i in 1:2) {
        n <- impacts[[i]]$calculate(n)
      }
      for (i in 1:3) {
        n <- actions[[i]]$clear_attributes(n) # clear
      }
      for (i in 1:3) {
        n <- actions[[i]]$apply(n, tm)
      }
      if (tm %in% c(0, 2, 4)) {
        expect_list[[r]][[1]]$detected[[tmc]] <-
          as.logical(attr(n, "1_detected"))
        expect_list[[r]][[1]]$number$detected[[tmc]] <-
          attr(attr(n, "1_detected"), "number")
        expect_list[[r]][[2]]$removed[[tmc]] <-
          as.logical(attr(n, "2_removed"))
        expect_list[[r]][[2]]$number$removed[[tmc]] <-
          attr(attr(n, "2_removed"), "number")
        expect_list[[r]][[3]]$control_growth[[tmc]] <-
          attr(n, "3_control_growth") < 1
        expect_list[[r]][[1]]$cost$detected[[tmc]] <-
          as.numeric(attr(n, "1_surv_cost"))
        expect_list[[r]][[2]]$cost$removed[[tmc]] <-
          as.numeric(attr(n, "2_removal_cost"))
        expect_list[[r]][[3]]$cost$control_growth[[tmc]] <-
          as.numeric(attr(n, "3_control_growth_cost"))
      }
      expect_list[[r]][[1]]$total[[tmc]] <- sum(attr(n, "1_detected"))
      expect_list[[r]][[1]]$number$total[[tmc]] <-
        sum(attr(attr(n, "1_detected"), "number"))
      expect_list[[r]][[2]]$total[[tmc]] <- sum(attr(n, "2_removed"))
      expect_list[[r]][[2]]$number$total[[tmc]] <-
        sum(attr(attr(n, "2_removed"), "number"))
      expect_list[[r]][[3]]$total[[tmc]] <-
        sum(attr(n, "3_control_growth") < 1)
      expect_list[[r]][[1]]$cost$total[[tmc]] <- sum(attr(n, "1_surv_cost"))
      expect_list[[r]][[2]]$cost$total[[tmc]] <- sum(attr(n, "2_removal_cost"))
      expect_list[[r]][[3]]$cost$total[[tmc]] <-
        sum(attr(n, "3_control_growth_cost"))
      expect_list[[r]][[1]]$cost$cumulative$total[[tmc]] <-
        (expect_list[[r]][[1]]$cost$cumulative$total[[tmc_prev]] +
           sum(attr(n, "1_surv_cost")))
      expect_list[[r]][[2]]$cost$cumulative$total[[tmc]] <-
        (expect_list[[r]][[2]]$cost$cumulative$total[[tmc_prev]] +
           sum(attr(n, "2_removal_cost")))
      expect_list[[r]][[3]]$cost$cumulative$total[[tmc]] <-
        (expect_list[[r]][[3]]$cost$cumulative$total[[tmc_prev]] +
           sum(attr(n, "3_control_growth_cost")))
      combined_cost <- as.numeric(
        attr(n, "1_surv_cost") + attr(n, "2_removal_cost") +
          attr(n, "3_control_growth_cost"))
      combined_cost_plus_impacts <-
        (combined_cost +
           rowSums(sapply(1:2, function(i) attr(n, "impacts")[[i]])))
      if (tm %in% c(0, 2, 4)) {
        expect_list[[r]]$cost$combined[[tmc]] <- combined_cost
        expect_list[[r]]$combined_cost_plus_impacts[[tmc]] <-
          combined_cost_plus_impacts
      }
      expect_list[[r]]$cost$total[[tmc]] <- sum(combined_cost)
      expect_list[[r]]$cost$cumulative$total[[tmc]] <-
        (expect_list[[r]]$cost$cumulative$total[[tmc_prev]] +
           sum(combined_cost))
      expect_list[[r]]$total_cost_plus_impacts[[tmc]] <-
        sum(combined_cost_plus_impacts)
      expect_list[[r]]$total_cum_cost_plus_impacts[[tmc]] <-
        (expect_list[[r]]$total_cum_cost_plus_impacts[[tmc_prev]] +
           sum(combined_cost_plus_impacts))
      if (tm %in% c(0, 2, 4)) {
        tmc_prev <- as.character(tm)
      } else {
        tmc <- as.character(tm + 1)
      }
      expect_list[[r]][[1]]$cost$cumulative$detected[[tmc]] <-
        (expect_list[[r]][[1]]$cost$cumulative$detected[[tmc_prev]] +
           as.numeric(attr(n, "1_surv_cost")))
      expect_list[[r]][[2]]$cost$cumulative$removed[[tmc]] <-
        (expect_list[[r]][[2]]$cost$cumulative$removed[[tmc_prev]] +
           as.numeric(attr(n, "2_removal_cost")))
      expect_list[[r]][[3]]$cost$cumulative$control_growth[[tmc]] <-
        (expect_list[[r]][[3]]$cost$cumulative$control_growth[[tmc_prev]] +
           as.numeric(attr(n, "3_control_growth_cost")))
      expect_list[[r]]$cost$cumulative$combined[[tmc]] <-
        expect_list[[r]]$cost$cumulative$combined[[tmc_prev]] + combined_cost
      expect_list[[r]]$combined_cum_cost_plus_impacts[[tmc]] <-
        (expect_list[[r]]$combined_cum_cost_plus_impacts[[tmc_prev]] +
           combined_cost_plus_impacts)
      if (r == 1) {
        results$collate(r = 1, tm = tm, n = n)
      }
    }
  }
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expect_equal(lapply(result_list$population, function(i) attributes(i)),
               lapply(result_list$population, function(i) NULL))
  expect_equal(result_list$actions[[1]][c("detected", "total")],
               expect_list[[1]][[1]][c("detected", "total")])
  expect_equal(result_list$actions[[1]]$number[c("detected", "total")],
               expect_list[[1]][[1]]$number[c("detected", "total")])
  expect_equal(result_list$actions[[1]]$cost[c("detected", "total")],
               expect_list[[1]][[1]]$cost[c("detected", "total")])
  expect_equal(
    result_list$actions[[1]]$cost$cumulative[c("detected", "total")],
    expect_list[[1]][[1]]$cost$cumulative[c("detected", "total")])
  expect_equal(result_list$actions[[2]][c("removed", "total")],
               expect_list[[1]][[2]][c("removed", "total")])
  expect_equal(result_list$actions[[2]]$number[c("removed", "total")],
               expect_list[[1]][[2]]$number[c("removed", "total")])
  expect_equal(result_list$actions[[2]]$cost[c("removed", "total")],
               expect_list[[1]][[2]]$cost[c("removed", "total")])
  expect_equal(result_list$actions[[2]]$cost$cumulative[c("removed", "total")],
               expect_list[[1]][[2]]$cost$cumulative[c("removed", "total")])
  expect_equal(result_list$actions[[3]][c("control_growth", "total")],
               expect_list[[1]][[3]][c("control_growth", "total")])
  expect_equal(result_list$actions[[3]]$cost[c("control_growth", "total")],
               expect_list[[1]][[3]]$cost[c("control_growth", "total")])
  expect_equal(
    result_list$actions[[3]]$cost$cumulative[c("control_growth", "total")],
    expect_list[[1]][[3]]$cost$cumulative[c("control_growth", "total")])
  expect_equal(result_list$actions$cost$combined,
               expect_list[[1]]$cost$combined)
  expect_equal(result_list$actions$cost$cumulative$combined,
               expect_list[[1]]$cost$cumulative$combined)
  expect_equal(result_list$actions$cost$total,
               expect_list[[1]]$cost$total)
  expect_equal(result_list$actions$cost$cumulative$total,
               expect_list[[1]]$cost$cumulative$total)
  expect_equal(attr(result_list$actions$cost, "unit"), "$")
  expect_equal(result_list$cost$combined,
               expect_list[[1]]$combined_cost_plus_impacts)
  expect_equal(result_list$cost$cumulative$combined,
               expect_list[[1]]$combined_cum_cost_plus_impacts)
  expect_equal(result_list$cost$total,
               expect_list[[1]]$total_cost_plus_impacts)
  expect_equal(result_list$cost$cumulative$total,
               expect_list[[1]]$total_cum_cost_plus_impacts)
  expect_equal(attr(result_list$cost, "unit"), "$")
  # multiple replicates
  expect_silent(
    results <- Results(region, population_model = population_model,
                       impacts = impacts, actions = actions,
                       time_steps = 4, collation_steps = 2,
                       replicates = 3))
  expect_silent(result_list_0 <- results$get_list())
  set.seed(1234)
  for (r in 1:3) {
    n <- n_list[[r]]
    for (tm in 0:4) {
      for (i in 1:2) {
        n <- impacts[[i]]$calculate(n)
      }
      for (i in 1:3) {
        n <- actions[[i]]$clear_attributes(n) # clear
      }
      for (i in 1:3) {
        n <- actions[[i]]$apply(n, tm)
      }
      results$collate(r = r, tm = tm, n = n)
    }
  }
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expect_collated <- result_list_0$actions[[1]]$number[[1]]
  expect_collated_no_sd <- lapply(expect_collated, function(ec) ec["mean"])
  expect_totals <- result_list_0$actions[[1]]$total
  generate_summary <- function(rep_list, summary_list) {
    for (i in 1:length(summary_list)) {
      nrow <- length(summary_list[[i]][[1]])
      rep_values <- matrix(sapply(rep_list, function(l) l[[i]]), nrow = nrow)
      summary_list[[i]]$mean <- rowMeans(rep_values)
      if (any(sapply(summary_list, function(sl) "sd" %in% names(sl)))) {
        summary_list[[i]]$sd <- apply(rep_values, 1, sd)
      }
    }
    return(summary_list)
  }
  expect_equal(
    lapply(result_list$population,
           function(i) lapply(i, function(j) attributes(j))),
    lapply(expect_collated, function(i) lapply(i, function(j) NULL)))
  expect_equal(
    result_list$actions[[1]]$detected,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[1]]$detected),
                     expect_collated_no_sd))
  expect_equal(
    result_list$actions[[1]]$total,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[1]]$total),
                     expect_totals))
  expect_equal(
    result_list$actions[[1]]$number$detected,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[1]]$number$detected),
                     expect_collated))
  expect_equal(
    result_list$actions[[1]]$number$total,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[1]]$number$total),
                     expect_totals))
  expect_equal(
    result_list$actions[[1]]$cost$detected,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[1]]$cost$detected),
                     expect_collated))
  expect_equal(
    result_list$actions[[1]]$cost$total,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[1]]$cost$total),
                     expect_totals))
  expect_equal(
    result_list$actions[[1]]$cost$cumulative$detected,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]][[1]]$cost$cumulative$detected),
      expect_collated))
  expect_equal(
    result_list$actions[[1]]$cost$cumulative$total,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]][[1]]$cost$cumulative$total),
      expect_totals))
  expect_equal(
    result_list$actions[[2]]$removed,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[2]]$removed),
                     expect_collated_no_sd))
  expect_equal(
    result_list$actions[[2]]$total,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[2]]$total),
                     expect_totals))
  expect_equal(
    result_list$actions[[2]]$number$removed,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[2]]$number$removed),
                     expect_collated))
  expect_equal(
    result_list$actions[[2]]$number$total,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[2]]$number$total),
                     expect_totals))
  expect_equal(
    result_list$actions[[2]]$cost$removed,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[2]]$cost$removed),
                     expect_collated))
  expect_equal(
    result_list$actions[[2]]$cost$total,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[2]]$cost$total),
                     expect_totals))
  expect_equal(
    result_list$actions[[2]]$cost$cumulative$removed,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]][[2]]$cost$cumulative$removed),
      expect_collated))
  expect_equal(
    result_list$actions[[2]]$cost$cumulative$total,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]][[2]]$cost$cumulative$total),
      expect_totals))
  expect_equal(
    result_list$actions[[3]]$control_growth,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[3]]$control_growth),
                     expect_collated_no_sd))
  expect_equal(
    result_list$actions[[3]]$total,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[3]]$total),
                     expect_totals))
  expect_equal(
    result_list$actions[[3]]$cost$control_growth,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]][[3]]$cost$control_growth),
      expect_collated))
  expect_equal(
    result_list$actions[[3]]$cost$total,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[3]]$cost$total),
                     expect_totals))
  expect_equal(
    result_list$actions[[3]]$cost$cumulative$control_growth,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]][[3]]$cost$cumulative$control_growth),
      expect_collated))
  expect_equal(
    result_list$actions[[3]]$cost$cumulative$total,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]][[3]]$cost$cumulative$total),
      expect_totals))
  expect_equal(
    result_list$actions$cost$combined,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]]$cost$combined),
                     expect_collated))
  expect_equal(
    result_list$actions$cost$cumulative$combined,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]]$cost$cumulative$combined),
      expect_collated))
  expect_equal(
    result_list$actions$cost$total,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]]$cost$total),
                     expect_totals))
  expect_equal(
    result_list$actions$cost$cumulative$total,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]]$cost$cumulative$total),
      expect_totals))
  expect_equal(attr(result_list$actions$cost, "unit"), "$")
  expect_equal(
    result_list$cost$combined,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]]$combined_cost_plus_impacts),
      expect_collated))
  expect_equal(
    result_list$cost$cumulative$combined,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]]$combined_cum_cost_plus_impacts),
      expect_collated))
  expect_equal(
    result_list$cost$total,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]]$total_cost_plus_impacts),
      expect_totals))
  expect_equal(
    result_list$cost$cumulative$total,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]]$total_cum_cost_plus_impacts),
      expect_totals))
  expect_equal(attr(result_list$cost, "unit"), "$")
})

test_that("collates and finalizes spatially implicit results with costs", {
  region <- Region()
  population_model <- UnstructPopulation(region, growth = 1.2)
  n_list <- list()
  n_list[[1]] <- 60; attr(n_list[[1]], "spread_area") <- 8000
  n_list[[2]] <- 80; attr(n_list[[2]], "spread_area") <- 10000
  n_list[[3]] <- 100; attr(n_list[[3]],"spread_area") <- 12000
  impacts_1 <- Impacts(region, population_model, # monetary
                       impact_type = "area",
                       asset_name = "impact1",
                       asset_value = 0.1,
                       loss_rate = 0.1)
  impacts_1$set_id(1)
  impacts_2 <- Impacts(region, population_model, # monetary
                       impact_type = "area",
                       asset_name = "impact2",
                       asset_value = 0.2,
                       loss_rate = 0.2)
  impacts_2$set_id(2)
  impacts <- list(impacts_1, impacts_2)
  surv_cost <- 300; attr(surv_cost, "unit") <- "$"
  actions_1 <- Detection(region, population_model,
                         sensitivity = 0.6,
                         surv_cost = surv_cost,
                         schedule = 2:3)
  actions_1$set_id(1)
  removal_cost <- 0.04; attr(removal_cost, "unit") <- "$"
  actions_2 <- Removals(region, population_model,
                        removal_pr = 0.3,
                        removal_cost = removal_cost,
                        schedule = 2:3)
  actions_2$set_id(2)
  control_cost <- 0.05; attr(control_cost, "unit") <- "$"
  actions_3 <- Controls(region, population_model,
                        control_type = "growth",
                        control_mult = 0.9,
                        control_cost = control_cost,
                        schedule = 2:3)
  actions_3$set_id(3)
  actions <- list(actions_1, actions_2, actions_3)
  # single replicate
  expect_silent(results <- Results(region, population_model = population_model,
                                   impacts = impacts, actions = actions,
                                   time_steps = 4, collation_steps = 2,
                                   replicates = 1))
  action_results_0 <- results$get_list()$actions
  action_results_0$combined_cost_plus_impacts <-
    action_results_0$cost$combined
  action_results_0$combined_cum_cost_plus_impacts <-
    action_results_0$cost$combined
  expect_list <- list(action_results_0, action_results_0, action_results_0)
  set.seed(1234)
  for (r in 1:3) { # collect expected for single and multiple replicates
    n <- n_list[[r]]
    for (tm in 0:4) {
      tmc <- as.character(tm)
      tmc_prev <- as.character(max(tm - 1, 0))
      for (i in 1:2) {
        n <- impacts[[i]]$calculate(n)
      }
      for (i in 1:3) {
        n <- actions[[i]]$clear_attributes(n) # clear
      }
      for (i in 1:3) {
        n <- actions[[i]]$apply(n, tm)
      }
      expect_list[[r]][[1]]$detected[[tmc]] <-
        as.logical(attr(n, "1_detected"))
      expect_list[[r]][[1]]$number$detected[[tmc]] <-
        attr(attr(n, "1_detected"), "number")
      expect_list[[r]][[2]]$removed[[tmc]] <- as.logical(attr(n, "2_removed"))
      expect_list[[r]][[2]]$number$removed[[tmc]] <-
        attr(attr(n, "2_removed"), "number")
      expect_list[[r]][[3]]$control_growth[[tmc]] <-
        attr(n, "3_control_growth") < 1
      expect_list[[r]][[2]]$cost$removed[[tmc]] <-
        as.numeric(attr(n, "2_removal_cost"))
      expect_list[[r]][[1]]$cost$detected[[tmc]] <-
        as.numeric(attr(n, "1_surv_cost"))
      expect_list[[r]][[3]]$cost$control_growth[[tmc]] <-
        as.numeric(attr(n, "3_control_growth_cost"))
      combined_cost <-
        (as.numeric(attr(n, "1_surv_cost") + attr(n, "2_removal_cost") +
                      attr(n, "3_control_growth_cost")))
      combined_cost_plus_impacts <-
        combined_cost + sum(unlist(attr(n, "impacts")))
      expect_list[[r]]$cost$combined[[tmc]] <- combined_cost
      expect_list[[r]]$combined_cost_plus_impacts[[tmc]] <-
        combined_cost_plus_impacts
      expect_list[[r]][[1]]$cost$cumulative$detected[[tmc]] <-
        (expect_list[[r]][[1]]$cost$cumulative$detected[[tmc_prev]] +
           as.numeric(attr(n, "1_surv_cost")))
      expect_list[[r]][[2]]$cost$cumulative$removed[[tmc]] <-
        (expect_list[[r]][[2]]$cost$cumulative$removed[[tmc_prev]] +
           as.numeric(attr(n, "2_removal_cost")))
      expect_list[[r]][[3]]$cost$cumulative$control_growth[[tmc]] <-
        (expect_list[[r]][[3]]$cost$cumulative$control_growth[[tmc_prev]] +
           as.numeric(attr(n, "3_control_growth_cost")))
      expect_list[[r]]$cost$cumulative$combined[[tmc]] <-
        expect_list[[r]]$cost$cumulative$combined[[tmc_prev]] + combined_cost
      expect_list[[r]]$combined_cum_cost_plus_impacts[[tmc]] <-
        (expect_list[[r]]$combined_cum_cost_plus_impacts[[tmc_prev]] +
           combined_cost_plus_impacts)
      if (r == 1) {
        results$collate(r = 1, tm = tm, n = n)
      }
    }
  }
  expect_silent(results$finalize())
  result_list <- results$get_list()
  expect_equal(result_list$actions[[1]], expect_list[[1]][[1]])
  expect_equal(result_list$actions[[2]], expect_list[[1]][[2]])
  expect_equal(result_list$actions[[3]], expect_list[[1]][[3]])
  expect_equal(result_list$actions$cost, expect_list[[1]]$cost)
  expect_equal(result_list$cost$combined,
               expect_list[[1]]$combined_cost_plus_impacts)
  expect_equal(result_list$cost$cumulative$combined,
               expect_list[[1]]$combined_cum_cost_plus_impacts)
  # multiple replicates
  expect_silent(
    results <- Results(region, population_model = population_model,
                       impacts = impacts, actions = actions,
                       time_steps = 4, collation_steps = 2,
                       replicates = 3))
  expect_silent(result_list_0 <- results$get_list())
  set.seed(1234)
  for (r in 1:3) {
    n <- n_list[[r]]
    for (tm in 0:4) {
      for (i in 1:2) {
        n <- impacts[[i]]$calculate(n)
      }
      for (i in 1:3) {
        n <- actions[[i]]$clear_attributes(n) # clear
      }
      for (i in 1:3) {
        n <- actions[[i]]$apply(n, tm)
      }
      results$collate(r = r, tm = tm, n = n)
    }
  }
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expect_collated <- result_list_0$actions[[1]]$number[[1]]
  expect_collated_no_sd <- lapply(expect_collated, function(ec) ec["mean"])
  generate_summary <- function(rep_list, summary_list) {
    for (i in 1:length(summary_list)) {
      nrow <- length(summary_list[[i]][[1]])
      rep_values <- matrix(sapply(rep_list, function(l) l[[i]]), nrow = nrow)
      summary_list[[i]]$mean <- rowMeans(rep_values)
      if (any(sapply(summary_list, function(sl) "sd" %in% names(sl)))) {
        summary_list[[i]]$sd <- apply(rep_values, 1, sd)
      }
    }
    return(summary_list)
  }
  expect_equal(
    result_list$actions[[1]]$detected,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[1]]$detected),
                     expect_collated_no_sd))
  expect_equal(
    result_list$actions[[1]]$number$detected,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[1]]$number$detected),
                     expect_collated))
  expect_equal(
    result_list$actions[[1]]$cost$detected,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[1]]$cost$detected),
                     expect_collated))
  expect_equal(
    result_list$actions[[1]]$cost$cumulative$detected,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]][[1]]$cost$cumulative$detected),
      expect_collated))
  expect_equal(
    result_list$actions[[2]]$removed,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[2]]$removed),
                     expect_collated_no_sd))
  expect_equal(
    result_list$actions[[2]]$number$removed,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[2]]$number$removed),
                     expect_collated))
  expect_equal(
    result_list$actions[[2]]$cost$removed,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[2]]$cost$removed),
                     expect_collated))
  expect_equal(
    result_list$actions[[2]]$cost$cumulative$removed,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]][[2]]$cost$cumulative$removed),
      expect_collated))
  expect_equal(
    result_list$actions[[3]]$control_growth,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]][[3]]$control_growth),
                     expect_collated_no_sd))
  expect_equal(
    result_list$actions[[3]]$cost$control_growth,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]][[3]]$cost$control_growth),
      expect_collated))
  expect_equal(
    result_list$actions[[3]]$cost$cumulative$control_growth,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]][[3]]$cost$cumulative$control_growth),
      expect_collated))
  expect_equal(
    result_list$actions$cost$combined,
    generate_summary(lapply(as.list(1:3),
                            function(i) expect_list[[i]]$cost$combined),
                     expect_collated))
  expect_equal(
    result_list$actions$cost$cumulative$combined,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]]$cost$cumulative$combined),
      expect_collated))
  expect_equal(attr(result_list$actions$cost, "unit"), "$")
  expect_equal(
    result_list$cost$combined,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]]$combined_cost_plus_impacts),
      expect_collated))
  expect_equal(
    result_list$cost$cumulative$combined,
    generate_summary(
      lapply(as.list(1:3),
             function(i) expect_list[[i]]$combined_cum_cost_plus_impacts),
      expect_collated))
  expect_equal(attr(result_list$cost, "unit"), "$")
})

test_that("collates and finalizes staged action results", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  idx <- 5916:5922
  region <- bsspread::Region(template)
  template[region$get_indices()][idx,] <- c(rep(0.5, 4), 0.5, 0.75, 1)
  idx <- idx[5:7]
  sensitivity <- removal_pr <- template[region$get_indices()][,1]
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  initial_n <- rep(0, region$get_locations())
  n <- rep(0, region$get_locations())
  n[idx] <- (10:12)*5
  n <- population_model$make(initial = n)
  actions_1 <- Detection(region, population_model,
                         sensitivity = sensitivity,
                         stages = 2:3, schedule = 2:3)
  actions_1$set_id(1)
  actions_2 <- Removals(region, population_model,
                        removal_pr = removal_pr,
                        stages = 2:3, schedule = 2:3)
  actions_2$set_id(2)
  actions_3 <- Controls(region, population_model,
                        control_type = "growth",
                        control_mult = 0.7,
                        stages = 2:3, schedule = 2:3)
  actions_3$set_id(3)
  actions <- list(actions_1, actions_2, actions_3)
  # staged single replicate
  expect_silent(
    results <- Results(region,
                       population_model = population_model,
                       actions = actions,
                       time_steps = 4, collation_steps = 2,
                       replicates = 1))
  set.seed(1234)
  n <- actions[[1]]$apply(n, 2)
  n <- actions[[2]]$apply(n, 2)
  n <- actions[[3]]$apply(n, 2)
  expect_silent(results$collate(r = 1, tm = 2, n = n))
  expect_silent(result_list <- results$get_list())
  action_results <- lapply(result_list$actions,
                           function(i) lapply(i[names(i)[1:2]],
                                              function(j) j[["2"]]))
  action_num_results <- lapply(result_list$actions[1:2],
                               function(i) lapply(i$number[names(i)[1:2]],
                                                  function(j) j[["2"]]))
  expect_equal(lapply(action_results[[1]], function(a) length(a)),
               list(detected = region$get_locations(), total = 1))
  expect_equal(lapply(action_num_results[[1]], function(a) dim(a)),
               list(detected = c(region$get_locations(), 3), total = c(1, 3)))
  expect_equal(lapply(action_results[[2]], function(a) length(a)),
               list(removed = region$get_locations(), total = 1))
  expect_equal(lapply(action_num_results[[2]], function(a) dim(a)),
               list(removed = c(region$get_locations(), 3), total = c(1, 3)))
  expect_equal(lapply(action_results[[3]], function(a) length(a)),
               list(control_growth = region$get_locations(), total = 1))
  expect_true(action_results[[1]]$total > 0)
  expect_equal(as.logical(action_num_results[[1]]$total > 0),
               c(FALSE, TRUE, TRUE))
  expect_true(action_results[[2]]$total > 0)
  expect_equal(as.logical(action_num_results[[2]]$total > 0),
               c(FALSE, TRUE, TRUE))
  expect_equal(action_results[[3]]$total,
               sum(action_results[[3]]$control_growth))
  # staged multiple replicates
  set.seed(1234)
  n_r <- list(); n_r2 <- list()
  n <- rep(0, region$get_locations())
  n[idx] <- (9:11)*5
  n <- population_model$make(initial = n)
  n <- actions[[1]]$apply(n, 2); n <- actions[[2]]$apply(n, 2)
  n <- actions[[3]]$apply(n, 2); n_r[[1]] <- n
  attributes(n)[c("1_detected", "undetected", "2_removed",
                  "3_control_growth")] <- NULL
  n <- actions[[1]]$apply(n, 4); n <- actions[[2]]$apply(n, 4)
  n <- actions[[3]]$apply(n, 4); n_r2[[1]] <- n
  n <- rep(0, region$get_locations())
  n[idx] <- (10:12)*5
  n <- population_model$make(initial = n)
  n <- actions[[1]]$apply(n, 2); n <- actions[[2]]$apply(n, 2)
  n <- actions[[3]]$apply(n, 2); n_r[[2]] <- n
  attributes(n)[c("1_detected", "undetected", "2_removed",
                  "3_control_growth")] <- NULL
  n <- actions[[1]]$apply(n, 4); n <- actions[[2]]$apply(n, 4)
  n <- actions[[3]]$apply(n, 4); n_r2[[2]] <- n
  n <- rep(0, region$get_locations())
  n[idx] <- (11:13)*5
  n <- population_model$make(initial = n)
  n <- actions[[1]]$apply(n, 2); n <- actions[[2]]$apply(n, 2)
  n <- actions[[3]]$apply(n, 2); n_r[[3]] <- n
  attributes(n)[c("1_detected", "undetected", "2_removed",
                  "3_control_growth")] <- NULL
  n <- actions[[1]]$apply(n, 4); n <- actions[[2]]$apply(n, 4)
  n <- actions[[3]]$apply(n, 4); n_r2[[3]] <- n
  collated_r <- lapply(n_r, function(n) lapply(actions, function(a) {
    if (a$get_label(include_id = FALSE) == "control_growth") {
      n_a <- +(attr(n, a$get_label()) < 1)
    } else {
      n_a <- as.numeric(attr(n, a$get_label()))
    }
    collated <- list(n_a, total = sum(n_a))
    names(collated)[1] <- a$get_label(include_id = FALSE)
    collated
  }))
  collated_arrays <- collated_r[[1]]
  for (i in 2:3) {
    for (j in names(actions)) {
      for (a in names(collated_r[[i]][[j]])) {
        collated_arrays[[j]][[a]] <- rbind(collated_arrays[[j]][[a]],
                                           collated_r[[i]][[j]][[a]])
      }
    }
  }
  expect_silent(
    results <- Results(region,
                       population_model = population_model,
                       actions = actions,
                       time_steps = 4, collation_steps = 2,
                       replicates = 3))
  expect_silent(results$collate(r = 1, tm = 2, n = n_r[[1]]))
  expect_silent(results$collate(r = 1, tm = 4, n = n_r2[[1]]))
  expect_silent(results$collate(r = 2, tm = 2, n = n_r[[2]]))
  expect_silent(results$collate(r = 2, tm = 4, n = n_r2[[2]]))
  expect_silent(results$collate(r = 3, tm = 2, n = n_r[[3]]))
  expect_silent(results$collate(r = 3, tm = 4, n = n_r2[[3]]))
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  action_results <- lapply(result_list$actions,
                           function(i) lapply(i[names(i)[1:2]],
                                              function(j) j[["2"]]))
  action_num_results <- lapply(result_list$actions[1:2],
                               function(i) lapply(i$number[names(i)[1:2]],
                                                  function(j) j[["2"]]))
  locs <- region$get_locations()
  expect_equal(lapply(action_results[[1]], function(a) lapply(a, length)),
               list(detected = list(mean = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(action_num_results[[1]], function(a) lapply(a, dim)),
               list(detected = list(mean = c(locs, 3), sd = c(locs, 3)),
                    total = list(mean = c(1, 3), sd = c(1, 3))))
  expect_equal(lapply(action_results[[2]], function(a) lapply(a, length)),
               list(removed = list(mean = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(action_num_results[[2]], function(a) lapply(a, dim)),
               list(removed = list(mean = c(locs, 3), sd = c(locs, 3)),
                    total = list(mean = c(1, 3), sd = c(1, 3))))
  expect_equal(lapply(action_results[[3]], function(a) lapply(a, length)),
               list(control_growth = list(mean = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(action_results[[1]]$total,
                      function(tot) as.logical(tot > 0)),
               list(mean = TRUE, sd = FALSE))
  expect_equal(lapply(action_num_results[[1]]$total,
                      function(tot) as.logical(tot > 0)),
               list(mean = c(FALSE, TRUE, TRUE), sd = c(FALSE, TRUE, TRUE)))
  expect_equal(lapply(action_results[[2]]$total,
                      function(tot) as.logical(tot > 0)),
               list(mean = TRUE, sd = FALSE))
  expect_equal(lapply(action_num_results[[2]]$total,
                      function(tot) as.logical(tot > 0)),
               list(mean = c(FALSE, TRUE, TRUE), sd = c(FALSE, TRUE, TRUE)))
  expect_equal(action_results[[3]]$total,
               list(mean = sum(action_results[[3]]$control_growth$mean),
                    sd = 0))
  # Population level probabilities
  set.seed(1234)
  manage_pr <- runif(region$get_locations())
  surv_cost <- 3; attr(surv_cost, "unit") <- "$"
  actions_1 <- Detection(region, population_model,
                         sensitivity = sensitivity,
                         sensitivity_type = "presence",
                         surv_cost = surv_cost,
                         stages = 2:3, schedule = 2:3)
  actions_1$set_id(1)
  removal_cost <- 4; attr(removal_cost, "unit") <- "$"
  actions_2 <- Removals(region, population_model,
                        removal_pr = removal_pr,
                        removal_pr_type = "population",
                        removal_cost = removal_cost,
                        stages = 2:3, schedule = 2:3)
  actions_2$set_id(2)
  control_cost <- 5; attr(control_cost, "unit") <- "$"
  actions_3 <- Controls(region, population_model,
                        control_type = "search_destroy",
                        manage_pr = manage_pr,
                        manage_pr_type = "population",
                        control_cost = control_cost,
                        sstages = 2:3, chedule = 2:3)
  actions_3$set_id(3)
  actions <- list(actions_1, actions_2, actions_3)
  expect_silent(results <- Results(region,
                     population_model = population_model,
                     actions = actions,
                     time_steps = 4, collation_steps = 2,
                     replicates = 3))
  zero_results <- results$get_list()$actions
  locs <- region$get_locations()
  expect_equal(lapply(zero_results[[1]][1:2],
                      function(a) lapply(a[["2"]], length)),
               list(detected = list(mean = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(zero_results[[1]]$cost[1:2],
                      function(a) lapply(a[["2"]], length)),
               list(detected = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(zero_results[[1]]$cost$cumulative,
                      function(a) lapply(a[["2"]], length)),
               list(detected = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(zero_results[[2]][1:2],
                      function(a) lapply(a[["2"]], length)),
               list(removed = list(mean = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(zero_results[[2]]$cost[1:2],
                      function(a) lapply(a[["2"]], length)),
               list(removed = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(zero_results[[2]]$cost$cumulative,
                      function(a) lapply(a[["2"]], length)),
               list(removed = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(zero_results[[3]][1:2],
                      function(a) lapply(a[["2"]], length)),
               list(control_search_destroy = list(mean = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(zero_results[[3]]$cost[1:2],
                      function(a) lapply(a[["2"]], length)),
               list(control_search_destroy = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(zero_results[[3]]$cost$cumulative,
                      function(a) lapply(a[["2"]], length)),
               list(control_search_destroy = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  set.seed(1234)
  n_r <- list(); n_r2 <- list()
  n <- rep(0, region$get_locations())
  n[idx] <- (9:11)*5
  n <- population_model$make(initial = n)
  n <- actions[[1]]$apply(n, 0); n <- actions[[2]]$apply(n, 0);
  n <- actions[[3]]$apply(n, 0); n_0 <- n
  attributes(n)[c("1_detected", "undetected", "2_removed",
                  "3_control_search_destroy")] <- NULL
  n <- actions[[1]]$apply(n, 2); n <- actions[[2]]$apply(n, 2)
  n <- actions[[3]]$apply(n, 2); n_r[[1]] <- n
  attributes(n)[c("1_detected", "undetected", "2_removed",
                  "3_control_search_destroy")] <- NULL
  n <- actions[[1]]$apply(n, 4); n <- actions[[2]]$apply(n, 4)
  n <- actions[[3]]$apply(n, 4); n_r2[[1]] <- n
  n <- rep(0, region$get_locations())
  n[idx] <- (10:12)*5
  n <- population_model$make(initial = n)
  n <- actions[[1]]$apply(n, 2); n <- actions[[2]]$apply(n, 2)
  n <- actions[[3]]$apply(n, 2); n_r[[2]] <- n
  attributes(n)[c("1_detected", "undetected", "2_removed",
                  "3_control_search_destroy")] <- NULL
  n <- actions[[1]]$apply(n, 4); n <- actions[[2]]$apply(n, 4)
  n <- actions[[3]]$apply(n, 4); n_r2[[2]] <- n
  n <- rep(0, region$get_locations())
  n[idx] <- (11:13)*5
  n <- population_model$make(initial = n)
  n <- actions[[1]]$apply(n, 2); n <- actions[[2]]$apply(n, 2)
  n <- actions[[3]]$apply(n, 2); n_r[[3]] <- n
  attributes(n)[c("1_detected", "undetected", "2_removed",
                  "3_control_search_destroy")] <- NULL
  n <- actions[[1]]$apply(n, 4); n <- actions[[2]]$apply(n, 4)
  n <- actions[[3]]$apply(n, 4); n_r2[[3]] <- n
  expect_silent(results$collate(r = 1, tm = 0, n = n_0))
  expect_silent(results$collate(r = 1, tm = 2, n = n_r[[1]]))
  expect_silent(results$collate(r = 1, tm = 4, n = n_r2[[1]]))
  expect_silent(results$collate(r = 2, tm = 0, n = n_0))
  expect_silent(results$collate(r = 2, tm = 2, n = n_r[[2]]))
  expect_silent(results$collate(r = 2, tm = 4, n = n_r2[[2]]))
  expect_silent(results$collate(r = 3, tm = 0, n = n_0))
  expect_silent(results$collate(r = 3, tm = 2, n = n_r[[3]]))
  expect_silent(results$collate(r = 3, tm = 4, n = n_r2[[3]]))
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expect_equal(lapply(result_list$actions[[1]][1:2],
                      function(a) lapply(a[["2"]], length)),
               list(detected = list(mean = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(result_list$actions[[1]]$cost[1:2],
                      function(a) lapply(a[["2"]], length)),
               list(detected = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(result_list$actions[[1]]$cost$cumulative,
                      function(a) lapply(a[["2"]], length)),
               list(detected = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(result_list$actions[[2]][1:2],
                      function(a) lapply(a[["2"]], length)),
               list(removed = list(mean = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(result_list$actions[[2]]$cost[1:2],
                      function(a) lapply(a[["2"]], length)),
               list(removed = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(result_list$actions[[2]]$cost$cumulative,
                      function(a) lapply(a[["2"]], length)),
               list(removed = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(result_list$actions[[3]][1:2],
                      function(a) lapply(a[["2"]], length)),
               list(control_search_destroy = list(mean = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(result_list$actions[[3]]$cost[1:2],
                      function(a) lapply(a[["2"]], length)),
               list(control_search_destroy = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(lapply(result_list$actions[[3]]$cost$cumulative,
                      function(a) lapply(a[["2"]], length)),
               list(control_search_destroy = list(mean = locs, sd = locs),
                    total = list(mean = 1, sd = 1)))
  expect_equal(result_list$actions[[1]]$total[["2"]]$mean,
               sum(result_list$actions[[1]]$detected[["2"]]$mean))
  expect_equal(result_list$actions[[1]]$cost$total[["2"]]$mean,
               sum(result_list$actions[[1]]$cost$detected[["2"]]$mean))
  expect_equal(result_list$actions[[2]]$total[["2"]]$mean,
               sum(result_list$actions[[2]]$removed[["2"]]$mean))
  expect_equal(result_list$actions[[2]]$cost$total[["2"]]$mean,
               sum(result_list$actions[[2]]$cost$removed[["2"]]$mean))
  expect_equal(
    result_list$actions[[3]]$total[["2"]]$mean,
    sum(result_list$actions[[3]]$control_search_destroy[["2"]]$mean))
  expect_equal(
    result_list$actions[[3]]$cost$total[["2"]]$mean,
    sum(result_list$actions[[3]]$cost$control_search_destroy[["2"]]$mean))
})

test_that("collates and finalizes spatially implicit staged action results", {
  region <- Region()
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  actions_1 <- Detection(region, population_model,
                         sensitivity = 0.7,
                         stages = 2:3, schedule = 2:3)
  actions_1$set_id(1)
  actions_2 <- Removals(region, population_model,
                        removal_pr = 0.8,
                        stages = 2:3, schedule = 2:3)
  actions_2$set_id(2)
  actions_3 <- Controls(region, population_model,
                        control_type = "growth",
                        control_mult = 0.7,
                        stages = 2:3, schedule = 2:3)
  actions_3$set_id(3)
  actions <- list(actions_1, actions_2, actions_3)
  # staged single replicate
  expect_silent(
    results <- Results(region, population_model = population_model,
                       actions = actions,
                       time_steps = 4, collation_steps = 2,
                       replicates = 1))
  n <- population_model$make(initial = 50)
  n <- actions[[1]]$apply(n, 2)
  n <- actions[[2]]$apply(n, 2)
  n <- actions[[3]]$apply(n, 2)
  expect_silent(results$collate(r = 1, tm = 2, n = n))
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  action_results <- lapply(result_list$actions,
                           function(i) lapply(i[names(i)[1]],
                                              function(j) j[["2"]]))
  action_num_results <- lapply(result_list$actions[1:2],
                               function(i) lapply(i$number[names(i)[1]],
                                                  function(j) j[["2"]]))
  expect_equal(lapply(action_results[[1]], function(a) length(a)),
               list(detected = 1))
  expect_equal(lapply(action_num_results[[1]], function(a) dim(a)),
               list(detected = c(1, 3)))
  expect_equal(lapply(action_results[[2]], function(a) length(a)),
               list(removed = 1))
  expect_equal(lapply(action_num_results[[2]], function(a) dim(a)),
               list(removed = c(1, 3)))
  expect_equal(lapply(action_results[[3]], function(a) length(a)),
               list(control_growth = 1))
  # staged multiple replicates
  set.seed(4321)
  n_r <- list()
  n <- population_model$make(initial = 40)
  n <- actions[[1]]$apply(n, 2)
  n <- actions[[2]]$apply(n, 2)
  n <- actions[[3]]$apply(n, 2)
  n_r[[1]] <- n
  n <- population_model$make(initial = 50)
  n <- actions[[1]]$apply(n, 2)
  n <- actions[[2]]$apply(n, 2)
  n <- actions[[3]]$apply(n, 2)
  n_r[[2]] <- n
  n <- population_model$make(initial = 60)
  n <- actions[[1]]$apply(n, 2)
  n <- actions[[2]]$apply(n, 2)
  n <- actions[[3]]$apply(n, 2)
  n_r[[3]] <- n
  collated_arrays <- list(
    list(detected = as.matrix(sapply(n_r, function(m) attr(m, "1_detected"))),
         number = t(sapply(n_r, function(m)
           attr(attr(m, "1_detected"), "number")))),
    list(removed = as.matrix(sapply(n_r, function(m) attr(m, "2_removed"))),
         number = t(sapply(n_r, function(m)
           attr(attr(m, "2_removed"), "number")))),
    list(control_growth = as.matrix(sapply(n_r, function(m)
      attr(m, "3_control_growth") < 1))))
  collated_summaries <- lapply(collated_arrays, function(a) {
    a_summary <- lapply(a, function(a_i) {
      list(mean = colMeans(a_i))
    })
    if ("number" %in% names(a)) {
      a_summary$number$mean <- t(as.matrix(a_summary$number$mean))
      a_summary$number$sd <- t(as.matrix(apply(a$number, 2, sd)))
      a_summary$number <- lapply(a_summary$number, function(m) {
        colnames(m) <- attr(population_model$get_growth(), "labels")
        m
      })
    }
    a_summary
  })
  expect_silent(results <- Results(region, population_model = population_model,
                                   actions = actions, time_steps = 4,
                                   collation_steps = 2, replicates = 3))
  expect_silent(results$collate(r = 1, tm = 2, n = n_r[[1]]))
  expect_silent(results$collate(r = 2, tm = 2, n = n_r[[2]]))
  expect_silent(results$collate(r = 3, tm = 2, n = n_r[[3]]))
  expect_silent(results$finalize())
  expect_silent(result_list <- results$get_list())
  expect_equal(result_list$actions[[1]]$detected[["2"]], collated_summaries[[1]]$detected)
  expect_equal(result_list$actions[[1]]$number$detected[["2"]], collated_summaries[[1]]$number)
  expect_equal(result_list$actions[[2]]$removed[["2"]], collated_summaries[[2]]$removed)
  expect_equal(result_list$actions[[2]]$number$removed[["2"]], collated_summaries[[2]]$number)
  expect_equal(result_list$actions[[3]]$control_growth[["2"]], collated_summaries[[3]]$control_growth)
})
