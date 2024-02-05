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
  expect_equal(names(result_list), c("collated", "total", "area"))
  expect_named(result_list$collated, as.character(seq(0, 10, 2)))
  expect_true(all(
    unlist(lapply(result_list$collated, length)) == region$get_locations()))
  expect_true(all(unlist(result_list$collated) == 0))
  expect_named(result_list$total, as.character(0:10))
  expect_true(all(unlist(lapply(result_list$total, length)) == 1))
  expect_true(all(unlist(result_list$total) == 0))
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
  expect_equal(lapply(result_list$collated, function(rl) rl[1:12]), expected)
  expect_equal(unname(unlist(result_list$total)), 1:11)
  expect_equal(unname(unlist(result_list$area)), (1:11)*1000^2)
  expect_equal(attr(result_list$area, "units"), "square metres")
})

test_that("collates and finalises multiple replicate results", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  expect_silent(results <- Results(region, population_model = population,
                                   time_steps = 10,
                                   step_duration = 1,
                                   step_units = "years",
                                   collation_steps = 2,
                                   replicates = 5,
                                   combine_stages = NULL))
  expect_silent(result_list <- results$get_list())
  expect_equal(unname(unlist(lapply(result_list$collated, names))),
               rep(c("mean", "sd"), 6))
  expect_equal(unname(unlist(lapply(result_list$total, names))),
               rep(c("mean", "sd"), 11))
  expect_equal(unname(unlist(lapply(result_list$area, names))),
               rep(c("mean", "sd"), 11))
  expect_true(all(unlist(lapply(
    result_list$collated,
    function(rl) lapply(rl, length))) == region$get_locations()))
  expect_true(all(unlist(result_list$collated) == 0))
  expect_length(unlist(result_list$total), 11*2)
  expect_true(all(unlist(result_list$total) == 0))
  expect_length(unlist(result_list$area), 11*2)
  expect_true(all(unlist(result_list$area) == 0))
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
  expect_equal(lapply(result_list$collated, function(rl) rl$mean[1:12]),
               lapply(expected_mask, function(e) e*expected_mean))
  expect_equal(lapply(result_list$collated, function(rl) rl$sd[1:12]),
               lapply(expected_mask, function(e) e*expected_sd))
  expect_equal(unname(unlist(lapply(result_list$total, function(rl) rl$mean))),
               sapply(1:11, function(a)
                 mean(colSums(expected[1:a,,drop = FALSE]))))
  expect_equal(unname(unlist(lapply(result_list$total, function(rl) rl$sd))),
               sapply(1:11, function(a)
                 stats::sd(colSums(expected[1:a,,drop = FALSE]))))
  expect_equal(unname(unlist(lapply(result_list$area, function(rl) rl$mean))),
               sapply(1:11, function(a)
                 mean(colSums(expected[1:a,,drop = FALSE] > 0)))*1000^2)
  expect_equal(unname(unlist(lapply(result_list$area, function(rl) rl$sd))),
               sapply(1:11, function(a)
                 stats::sd(colSums(expected[1:a,,drop = FALSE] > 0)))*1000^2)
  expect_equal(attr(result_list$area, "units"), "square metres")
})

test_that("collates spatially implicit area results", {
  region <- Region()
  population <- Population(region)
  expect_silent(results <- Results(region, population_model = population,
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
  expect_equal(unname(unlist(result_list$area)), pi*((0:10)*2000)^2)
  expect_equal(attr(result_list$area, "units"), "square metres")
})
