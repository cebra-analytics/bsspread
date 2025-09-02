context("AreaSpread")

test_that("initializes with region and population model", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  expect_error(area_spread <- AreaSpread(Region(template), Population(region)),
               "The region must be spatially implicit for area spread.")
  region <- Region()
  expect_error(area_spread <- AreaSpread(region, Population(region)),
               paste("The population model must be unstructured or stage",
                     "structured for area spread."))
  population_model <- UnstructPopulation(region)
  expect_error(area_spread <- AreaSpread(region, population_model),
               "The population capacity must be specified for area spread.")
  population_model <- UnstructPopulation(region,
                                         growth = 2,
                                         capacity = 100,
                                         capacity_area = 1e6)
  expect_error(AreaSpread(region, population_model, allow_contraction = 0),
               paste("The allow contraction indicator must be logical TRUE or",
                     "FALSE."))
  expect_silent(area_spread <- AreaSpread(region, population_model,
                                          allow_contraction = TRUE))
  expect_is(area_spread, "AreaSpread")
  expect_s3_class(area_spread, "Dispersal")
})

test_that("performs implicit area spread for an unstructured population", {
  TEST_DIRECTORY <- test_path("test_inputs")
  region <- Region()
  region$set_max_implicit_area(1e8)
  population_model <- UnstructPopulation(region,
                                         growth = 2,
                                         capacity = 100,
                                         capacity_area = 1e6)
  expect_silent(area_spread <- AreaSpread(region, population_model))
  n <- 100
  n_t <- list()
  set.seed(1234)
  expect_silent({
    for (i in 1:12) {
      n <- population_model$grow(n, tm = i)
      n <- area_spread$pack(n) # silent
      n <- area_spread$disperse(n, tm = i)  # silent
      n <- area_spread$unpack(n) # silent
      n_t[[i]] <- n
    }
  })
  pop <- unlist(n_t)
  delta_pop <- pop[1:12] - c(100, pop[1:11])
  expect_true(all(delta_pop[2:6] > delta_pop[1:5]))
  expect_true(all(delta_pop[8:12] < delta_pop[7:11]))
  expect_equal(sapply(n_t, function(n) attr(n, "spread_area")),
               pmin(pop*1e6/100, 1e8))
  n_decline <- n_t[[10]] - 1000
  expect_silent(n <- area_spread$pack(n_decline))
  expect_silent(n <- area_spread$disperse(n, tm = 11))
  expect_silent(n <- area_spread$unpack(n))
  expect_equal(attr(n, "spread_area"), as.numeric(n)*1e6/100)
  expect_silent(area_spread <- AreaSpread(region, population_model,
                                          allow_contraction = FALSE))
  expect_silent(n <- area_spread$pack(n_decline))
  expect_silent(n <- area_spread$disperse(n, tm = 11))
  expect_silent(n <- area_spread$unpack(n))
  expect_equal(attr(n, "spread_area"), attr(n_decline, "spread_area"))
})

test_that("performs implicit area spread for a staged population", {
  TEST_DIRECTORY <- test_path("test_inputs")
  region <- Region()
  region$set_max_implicit_area(1e8)
  stage_matrix <- matrix(c(0.0, 4.0, 10.0,
                           0.4, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region,
                                       growth = stage_matrix,
                                       capacity = 100,
                                       capacity_area = 1e6,
                                       capacity_stages = 2:3)
  expect_silent(area_spread <- AreaSpread(region, population_model))
  set.seed(1234)
  n <- population_model$make(initial = 500)
  n_t <- list()
  expect_silent({
    for (i in 1:12) {
      n <- population_model$grow(n, tm = i)
      n <- area_spread$pack(n) # silent
      n <- area_spread$disperse(n, tm = i)  # silent
      n <- area_spread$unpack(n) # silent
      n_t[[i]] <- n
    }
  })
  pop <- colSums(sapply(n_t, function(x) x)[2:3,])
  delta_pop <- pop[2:12] - pop[1:11]
  expect_true(all(delta_pop[2:6] > delta_pop[1:5]))
  expect_true(all(delta_pop[7:10] < delta_pop[6:9]))
  expect_equal(sapply(n_t, function(n) attr(n, "spread_area")),
               pmin(pop*1e6/100, 1e8))
})
