context("Simulator")

test_that("initializes with components and parameters", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  expect_error(simulator <- Simulator(region, time_steps = 0),
               paste("Time step and replicate parameters should be numeric",
                     "values > 0."))
  expect_error(simulator <- Simulator(region, step_duration = 0),
               paste("Time step and replicate parameters should be numeric",
                     "values > 0."))
  expect_error(simulator <- Simulator(region, collation_steps = 0),
               paste("Time step and replicate parameters should be numeric",
                     "values > 0."))
  expect_error(simulator <- Simulator(region, replicates = 0),
               paste("Time step and replicate parameters should be numeric",
                     "values > 0."))
  expect_error(simulator <- Simulator(region, parallel_cores = 0),
               "The number of parallel cores should be a numeric value > 0.")
  expect_error(simulator <- Simulator(region, initializer = 0),
               paste("Initializer must be a 'Initializer' or inherited class",
                     "object."))
  expect_error(simulator <- Simulator(region, population_model = 0),
               paste("Population model must be a 'Population' or inherited",
                     "class object."))
  expect_error(simulator <- Simulator(region, dispersal_models = 0),
               paste("Dispersal models must be 'Dispersal' or inherited class",
                     "objects."))
  expect_error(simulator <- Simulator(region, user_function = 0),
               paste("User-defined function must define a function for",
                     "transforming the population at each simulation time",
                     "step."))
  staged_population <- StagedPopulation(region, growth = array(1, c(3, 3)))
  expect_error(simulator <- Simulator(region,
                                      population_model = staged_population,
                                      result_stages = 1:4),
               paste("Result stages must be a numeric vector of stage indices",
                     "consistent with that defined in the population model."))
  expect_silent(simulator <- Simulator(region,
                                       population_model = staged_population,
                                       result_stages = 2:3))
  population <- Population(region)
  initializer <- Initializer(region$get_template(),
                             population_model = population)
  dispersal <- Dispersal(region, population)
  expect_silent(simulator <- Simulator(
    region,
    time_steps = 10,
    step_duration = 1,
    step_units = "years",
    collation_steps = 2,
    replicates = 5,
    initializer = initializer,
    population_model = population,
    dispersal_models = list(dispersal),
    user_function = function(n, r, tm) n + 1))
})

test_that("runs simulator with correct configuration", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  simulator <- Simulator(region)
  expect_error(simulator$run(),
               paste("The simulator requires the initializer and population",
                     "model to be set."))
  population <- UnstructPopulation(region, growth = 2)
  initial_n <- rep(0, region$get_locations())
  initial_n[5922] <- 10
  initializer <- Initializer(initial_n, region = region,
                             population_model = population)
  dispersal <- Dispersal(region, population, proportion = 1,
                         max_distance = 1000)
  expect_silent(simulator <- Simulator(region,
                                       initializer = initializer,
                                       population_model = population,
                                       dispersal_models = list(dispersal),
                                       user_function = function(n, r, tm) {
                                         idx <- which(n > 0)
                                         n[idx] <- n[idx] + 1
                                         return(n)
                                       }))
  set.seed(2143); n_grown <- stats::rpois(1, 2*10); set.seed(2143)
  expect_silent(results <- simulator$run())
  expect_is(results, "Results")
  results_list <- results$get_list()
  expect_named(results_list$collated, as.character(0:1))
  expect_equal(results_list$collated[["0"]][5922], 10)
  expect_true(all(results_list$collated[["0"]][-5922] == 0))
  expect_equal(results_list$collated[["1"]][5922], 0)
  idx <- which(results_list$collated[["1"]] > 0)
  expect_true(length(idx) <= 4)
  expect_equal(sum(results_list$collated[["1"]][idx]), n_grown + length(idx))
})

test_that("attaches attributes for spatially implicit diffusion", {
  region <- Region()
  population <- UnstructPopulation(region, growth = 2)
  initial_n <- 10
  initializer <- Initializer(initial_n, region = region,
                             population_model = population)
  diffusion <- Diffusion(region, population_model = population,
                         diffusion_rate = 2000, proportion = 1)
  expect_silent(simulator <- Simulator(region,
                                       initializer = initializer,
                                       population_model = population,
                                       dispersal_models = list(diffusion)))
  set.seed(1248); expect_n_1 <- stats::rpois(1, 2*initial_n); set.seed(1248)
  expect_diff_radius_1 <- sqrt(4*2000^2/(4*log(2))*1*log(expect_n_1/10))
  expect_silent(results <- simulator$run())
  results_list <- results$get_list()
  expect_named(results_list$collated, as.character(0:1))
  expect_equal(attributes(results_list$collated[["0"]]),
               list(initial_n = 10, diffusion_rate = 2000,
                    diffusion_radius = 0))
  expect_equal(attributes(results_list$collated[["1"]]),
               list(initial_n = 10, diffusion_rate = 2000,
                    diffusion_radius = expect_diff_radius_1, tm = 1))
  expect_equal(results_list$area[["0"]], 0)
  expect_equal(results_list$area[["1"]], pi*expect_diff_radius_1^2)
})

test_that("attaches attributes for spatially implicit area spread", {
  region <- Region()
  region$set_max_implicit_area(1e8)
  population_model <- UnstructPopulation(region,
                                         growth = 2,
                                         capacity = 100,
                                         capacity_area = 1e6)
  initial_n <- 100
  initializer <- Initializer(initial_n, region = region,
                             population_model = population_model)
  area_spread <- AreaSpread(region, population_model)
  area_capacity <- 100*1e8/1e6
  r <- exp(log(2)*(1 - initial_n/area_capacity))
  set.seed(1248)
  expect_n_1 <- stats::rpois(1, r*initial_n)
  expect_silent(simulator <- Simulator(region,
                                       initializer = initializer,
                                       population_model = population_model,
                                       dispersal_models = list(area_spread)))
  set.seed(1248)
  expect_silent(results <- simulator$run())
  results_list <- results$get_list()
  expect_named(results_list$collated, as.character(0:1))
  expect_equal(attributes(results_list$collated[["0"]]),
               list(spread_area = 1e6))
  expect_equal(attributes(results_list$collated[["1"]]),
               list(spread_area = expect_n_1*1e6/100))
  expect_equal(results_list$area[["0"]], 1e6)
  expect_equal(results_list$area[["1"]], expect_n_1*1e6/100)
  stage_matrix <- matrix(c(0.0, 4.0, 10.0,
                           0.4, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region,
                                       growth = stage_matrix,
                                       capacity = 100,
                                       capacity_area = 1e6,
                                       capacity_stages = 2:3)
  initializer <- Initializer(initial_n, region = region,
                             population_model = population_model)
  area_spread <- AreaSpread(region, population_model)
  expect_silent(simulator <- Simulator(region,
                                       initializer = initializer,
                                       population_model = population_model,
                                       dispersal_models = list(area_spread)))
  set.seed(1248)
  expect_silent(results <- simulator$run())
  results_list <- results$get_list()
  expect_equal(attr(results_list$collated[["0"]], "spread_area"),
               sum(results_list$total[["0"]][,2:3])*1e6/100)
  expect_equal(attr(results_list$collated[["1"]], "spread_area"),
               sum(results_list$total[["1"]][,2:3])*1e6/100)
  expect_equal(results_list$area[["0"]],
               sum(results_list$total[["0"]][,2:3])*1e6/100)
  expect_equal(results_list$area[["1"]],
               sum(results_list$total[["1"]][,2:3])*1e6/100)
})
