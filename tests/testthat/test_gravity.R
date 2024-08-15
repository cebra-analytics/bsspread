context("Gravity")

test_that("initializes with region, population model, and other parameters", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  expect_silent(gravity <- Gravity(region, population_model = population,
                                   attractors = list()))
  expect_is(gravity, "Gravity")
  expect_s3_class(gravity, "Dispersal")
  expect_error(gravity <- Gravity(region, population_model = population,
                                  attractors = list(a = 0)),
               paste("Attractors must be a list containing zero or more",
                     "'Attractor' or inherited class objects."))
  attractor_vect <- rep(0, region$get_locations())
  attractor_vect[1:5922] <- 1
  attractor <- Attractor(attractor_vect, region)
  expect_error(gravity <- Gravity(region, population_model = population,
                                  attractors = list(attractor),
                                  attractor_function = 0),
               "The attractor function should be defined as a function.")
  expect_error(gravity <- Gravity(region, population_model = population,
                                  attractors = list(attractor),
                                  attractor_function = function(a) 0),
               paste("The attractor function should transform a list of",
                     "'Attractor' or inherited class objects and return",
                     "another list of class objects."))
  expect_error(gravity <- Gravity(region, population_model = population,
                                  attractors = list(attractor),
                                  beta = 0),
               "The beta parameter must be numeric and > 0.")
  expect_error(gravity <- Gravity(region, population_model = population,
                                  attractors = list(attractor),
                                  distance_scale = 0),
               "The distance scale parameter must be numeric and > 0.")
  expect_silent(gravity <- Gravity(region, population_model = population,
                                   attractors = list(attractor),
                                   attractor_function = function(attractors){
                                     a_vector <- rep(0, region$get_locations())
                                     for (a in attractors) {
                                       a_vector <- a_vector + a$get_values()
                                     }
                                     list(Attractor(a_vector, region))
                                   },
                                   beta = 2, distance_scale = 1000,
                                   proportion = 1, events = 5,
                                   density_dependent = TRUE))
  expect_silent(distance_function <-
                  get("distance_function",
                      envir = environment(gravity$disperse)))
  expect_equal(distance_function(c(5000, 10000)), 1/c(5, 10)^2)
})

test_that("disperses via gravity function", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[1] <- TRUE
  attractor_vect <- rep(0, region$get_locations())
  attractor_vect[1:8] <- 5000
  attractor <- Attractor(attractor_vect, region, type = "destination")
  expect_silent(gravity <- Gravity(
    region, population_model = population, attractors = list(attractor),
    attractor_function = function(attractors){
      a_vector <- rep(0, region$get_locations())
      for (a in attractors) {
        a_vector <- a_vector + a$get_values()
      }
      list(Attractor(a_vector, region, type = "destination"))
    }, beta = 2, distance_scale = 1000, proportion = 0.8))
  expect_silent(n <- gravity$pack(n))
  region$calculate_paths(1)
  paths <- region$get_paths(1)
  idx <- paths$idx[["1"]]
  distances <-paths$distances[["1"]]
  set.seed(1259)
  new_locs <- stats::rbinom(13, size = 1,
                            prob = 0.8*attractor_vect[idx]/(distances/1000)^2)
  set.seed(1259)
  expect_silent(new_n <- gravity$unpack(gravity$disperse(n)))
  expect_equal(+(new_n), c(1, new_locs))
})
