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
  expect_error(simulator <- Simulator(region, time_steps = 1,
                                      impacts = list(1:2)),
               "Impacts must be a list of 'Impacts' objects.")
  expect_error(simulator <- Simulator(region, time_steps = 1,
                                      actions = list(1:2)),
               "Actions must be a list of 'Actions' objects.")
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
  population_model <- UnstructPopulation(region)
  initializer <- Initializer(region$get_template(),
                             population_model = population_model)
  dispersal <- Dispersal(region, population_model)
  impacts <- list(Impacts(region, population_model, # monetary
                          asset_name = "impact1",
                          asset_value = 10*(template < 0.9),
                          loss_rate = 0.1),
                  Impacts(region, population_model, # monetary
                          asset_name = "impact2",
                          asset_value = 20*(template > 0.6),
                          loss_rate = 0.2))
  sensitivity <- removal_pr <- template[region$get_indices()][,1]
  actions <- list(Detection(region, population_model,
                            sensitivity = sensitivity,
                            schedule = 2:3),
                  Removals(region, population_model,
                           removal_pr = removal_pr,
                           schedule = 2:3),
                  Controls(region, population_model,
                           control_type = "growth",
                           control_mult = 0.7,
                           schedule = 2:3))
  expect_silent(simulator <- Simulator(
    region,
    time_steps = 10,
    step_duration = 1,
    step_units = "years",
    collation_steps = 2,
    replicates = 5,
    initializer = initializer,
    population_model = population_model,
    dispersal_models = list(dispersal),
    impacts = impacts,
    actions = actions,
    user_function = function(n, r, tm) n + 1))
  expect_is(simulator, "Simulator")
  expect_named(simulator, c("set_initializer", "set_population_model",
                            "set_dispersal_models", "run"))
  expect_equal(sapply(impacts, function(i) i$get_id()), 1:2)
  expect_equal(sapply(actions, function(i) i$get_id()), 1:3)
})

test_that("runs simulator with correct configuration - growth & spread", {
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
  expect_named(results_list$population, as.character(0:1))
  expect_equal(results_list$population[["0"]][5922], 10)
  expect_true(all(results_list$population[["0"]][-5922] == 0))
  expect_equal(results_list$population[["1"]][5922], 0)
  idx <- which(results_list$population[["1"]] > 0)
  expect_true(length(idx) <= 4)
  expect_equal(sum(results_list$population[["1"]][idx]), n_grown + length(idx))
})

test_that("runs simulator with correct configuration - impacts & actions", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  idx <- 5920:5922
  region <- Region(template)
  simulator <- Simulator(region)
  expect_error(simulator$run(),
               paste("The simulator requires the initializer and population",
                     "model to be set."))
  template[region$get_indices()][idx,] <- c(0.25, 0.5, 1)
  template_vect <- template[region$get_indices()][,1]
  population_model <- UnstructPopulation(region, growth = 2)
  population_model_grow <- population_model$grow
  population_model$grow <- function(n, tm) {
    attr(n, "growth") <- c(attr(n, "growth"), max(c(attr(n, "growth"), 0) + 1))
    population_model_grow(n, tm)
  }
  initial_n <- rep(0, region$get_locations())
  initial_n[idx] <- (10:12)*5
  initializer <- Initializer(initial_n, region = region,
                             population_model = population_model)
  dispersal <- Dispersal(region, population_model)
  dispersal_disperse <- dispersal$disperse
  dispersal$disperse <- function(n, tm) {
    attr(n$relocated, "dispersal") <-
      c(attr(n$relocated, "dispersal"),
        max(c(attr(n$relocated, "dispersal"), 0) + 1))
    dispersal_disperse(n, tm)
  }
  impacts <- list(Impacts(region, population_model, # monetary
                          asset_name = "impact1",
                          asset_value = 100*template,
                          loss_rate = 0.5),
                  Impacts(region, population_model, # monetary
                          asset_name = "impact2",
                          asset_value = 200*template,
                          loss_rate = 0.7,
                          recovery_delay = 2))
  sensitivity <- removal_pr <- template[region$get_indices()][,1]
  actions <- list(Detection(region, population_model,
                            sensitivity = sensitivity),
                  Removals(region, population_model,
                           removal_pr = removal_pr),
                  Controls(region, population_model,
                           control_type = "growth",
                           control_mult = 0.7))
  expect_silent(
    simulator <- Simulator(
      region,
      time_steps = 4,
      collation_steps = 2,
      replicates = 1,
      initializer = initializer,
      population_model = population_model,
      dispersal_models = list(dispersal),
      impacts = impacts,
      actions = actions,
      user_function = function(n, r, tm) {
        attr(n, "user") <- c(attr(n, "user"), max(c(attr(n, "user"), 0) + 1))
        attr(n, "tm") <- c(attr(n, "tm"), tm)
        n
      }))
  set.seed(4321)
  expect_silent(results <- simulator$run())
  expect_silent(results_list <- results$get_list())
  expect_named(results_list, c("population", "total_pop", "occupancy",
                               "total_occup", "area", "impacts", "actions"))
  expect_named(results_list$population, as.character(seq(0, 4, 2)))
  expect_equal(unname(sapply(results_list$population, length)),
               rep(region$get_locations(), 3))
  expect_equal(attributes(results_list$population[["2"]]),
               list(growth = 1:2, dispersal = 1:2, user = 1:2, tm = 1:2))
  expect_equal(attributes(results_list$population[["4"]]),
               list(growth = 1:4, dispersal = 1:4, user = 1:4, tm = 1:4))
  expect_named(results_list$total_pop, as.character(0:4))
  expect_named(results_list$area, as.character(0:4))
  expect_named(results_list$occupancy, as.character(seq(0, 4, 2)))
  expect_named(results_list$total_occup, as.character(0:4))
  expect_named(results_list$impacts, c("idx", "monetary"))
  expect_equal(results_list$impacts$idx,
               list(monetary = c(impact1 = 1, impact2 = 2)))
  expect_named(results_list$impacts$monetary,
               c("total", "impact1", "impact2", "combined", "cumulative"))
  expect_named(results_list$impacts$monetary$total,
               c("impact1", "impact2", "combined"))
  expect_named(results_list$impacts$monetary$cumulative,
               c("total", "impact1", "impact2", "combined"))
  expect_named(results_list$impacts$monetary$cumulative$total,
               c("impact1", "impact2", "combined"))
  tmc_collated <- as.character(seq(0, 4, 2))
  tmc_all <- as.character(0:4)
  locs <- region$get_locations()
  expect_equal(lapply(results_list$impacts$monetary[2:4], names),
               list(impact1 = tmc_collated, impact2 = tmc_collated,
                    combined = tmc_collated))
  expect_equal(lapply(results_list$impacts$monetary[2:4],
                      function(i) length(i[["2"]])),
               list(impact1 = locs, impact2 = locs, combined = locs))
  expect_equal(lapply(results_list$impacts$monetary$total[1:3], names),
               list(impact1 = tmc_all, impact2 = tmc_all, combined = tmc_all))
  expect_equal(lapply(results_list$impacts$monetary$total[1:3],
                      function(i) length(i[["2"]])),
               list(impact1 = 1, impact2 = 1, combined = 1))
  expect_equal(lapply(results_list$impacts$monetary$cumulative[2:4], names),
               list(impact1 = tmc_collated, impact2 = tmc_collated,
                    combined = tmc_collated))
  expect_equal(lapply(results_list$impacts$monetary$cumulative[2:4],
                      function(i) length(i[["2"]])),
               list(impact1 = locs, impact2 = locs, combined = locs))
  expect_equal(lapply(results_list$impacts$monetary$cumulative$total[1:3],
                      names),
               list(impact1 = tmc_all, impact2 = tmc_all, combined = tmc_all))
  expect_equal(lapply(results_list$impacts$monetary$cumulative$total[1:3],
                      function(i) length(i[["2"]])),
               list(impact1 = 1, impact2 = 1, combined = 1))
  impact_mask <- unname(results_list$occupancy)
  impact_mask[[1]] <- +(initial_n > 0) # impacts before removal
  expect_equal(sapply(results_list$impacts$monetary[2:4],
                      function(i) attr(i, "unit")),
               c(impact1 = "$", impact2 = "$", combined = "$"))
  results_list$impacts$monetary[2:4] <-
    lapply(results_list$impacts$monetary[2:4], function(i) {
      attr(i, "unit") <- NULL
      i
    })
  expect_equal(unname(results_list$impacts$monetary$impact1),
               lapply(1:3, function (i)
                 (initial_n > 0)*template_vect*100*0.5*impact_mask[[i]]))
  impact_mask[[2]] <- impact_mask[[1]] # delayed impacts
  expect_equal(unname(results_list$impacts$monetary$impact2),
               lapply(1:3, function (i)
                 (initial_n > 0)*template_vect*200*0.7*impact_mask[[i]]))
  expect_equal(lapply(results_list$actions, names),
               list(c("detected", "total", "number"),
                    c("removed", "total", "number"),
                    c("control_growth", "total")))
  expect_equal(lapply(results_list$actions, function(a) names(a[[1]])),
               list(tmc_collated, tmc_collated, tmc_collated))
  expect_equal(lapply(results_list$actions, function(a) length(a[[1]][["2"]])),
               list(locs, locs, locs))
  expect_equal(lapply(results_list$actions, function(a) names(a$total)),
               list(tmc_all, tmc_all, tmc_all))
  expect_equal(lapply(results_list$actions,
                      function(a) length(a$total[["2"]])),
               list(1, 1, 1))
  expect_equal(lapply(results_list$actions[1:2], function(a) names(a$number)),
               list(c("detected", "total"), c("removed", "total")))
  expect_equal(lapply(results_list$actions[1:2],
                      function(a) names(a$number[[1]])),
               list(tmc_collated, tmc_collated))
  expect_equal(lapply(results_list$actions[1:2],
                      function(a) length(a$number[[1]][["2"]])),
               list(locs, locs))
  expect_equal(lapply(results_list$actions[1:2],
                      function(a) names(a$number$total)),
               list(tmc_all, tmc_all))
  expect_equal(lapply(results_list$actions[1:2],
                      function(a) length(a$number$total[["2"]])),
               list(1, 1))
  action_mask <- unname(lapply(results_list$occupancy, function(n) n > 0))
  action_mask[[1]] <- initial_n > 0 # detection before removal
  expect_equal(unname(results_list$actions[[1]]$detected), action_mask)
  expect_equal(unname(lapply(results_list$actions[[1]]$number$detected,
                             function(n) n > 0)), action_mask)
  expect_equal(unname(lapply(results_list$actions[[2]]$number$removed,
                             function(n) n > 0)), action_mask)
  expect_equal(results_list$actions[[2]]$removed[["0"]],
               initial_n > 0 & !results_list$population[["0"]])
  expect_equal(results_list$actions[[2]]$number$removed[["0"]],
               initial_n - results_list$population[["0"]])
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
  expect_named(results_list$population, as.character(0:1))
  expect_equal(attributes(results_list$population[["0"]]),
               list(initial_n = 10, diffusion_rate = 2000,
                    diffusion_radius = 0))
  expect_equal(attributes(results_list$population[["1"]]),
               list(initial_n = 10, diffusion_rate = 2000,
                    diffusion_radius = expect_diff_radius_1))
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
  expect_named(results_list$population, as.character(0:1))
  expect_equal(attributes(results_list$population[["0"]]),
               list(spread_area = 1e6))
  expect_equal(attributes(results_list$population[["1"]]),
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
  expect_equal(attr(results_list$population[["0"]], "spread_area"),
               sum(results_list$population[["0"]][,2:3])*1e6/100)
  expect_equal(attr(results_list$population[["1"]], "spread_area"),
               sum(results_list$population[["1"]][,2:3])*1e6/100)
  expect_equal(results_list$area[["0"]],
               sum(results_list$population[["0"]][,2:3])*1e6/100)
  expect_equal(results_list$area[["1"]],
               sum(results_list$population[["1"]][,2:3])*1e6/100)
})
