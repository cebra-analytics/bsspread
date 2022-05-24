context("Dispersal")

test_that("initializes with region, population model, and other parameters", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  expect_silent(dispersal <- Dispersal(region, population_model = population))
  expect_is(dispersal, "Dispersal")
  expect_error(dispersal <- Dispersal(region, population_model = 1:10),
               paste("Population model must be a 'Population' or inherited",
                     "class object."))
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      proportion = 0),
               "The proportion parameter must be numeric and > 0.")
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      events = 0),
               "The events parameter must be numeric and > 0.")
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      distance_function = 0),
               "The distance function must be a function.")
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      distance_adjust = "y"),
               "The distance adjust must be logical.")
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      direction_function = 0),
               "The direction function must be a function.")
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      attractors = 0),
               paste("Attractors must be a list containing zero or more",
                     "'Attractor' or inherited class objects, and/or a",
                     "numeric attractor (> 0) named 'source_density'."),
               fixed = TRUE)
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      permeability = 0),
               paste("Permeability parameter must be a 'Permeability' or",
                     "inherited class object."))
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      max_distance = 0),
               "The maximum distance parameter must be numeric and > 0.")
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       proportion = 1, events = 5,
                                       distance_function = function(x) 1/x,
                                       distance_adjust = FALSE,
                                       direction_function = function(x) x/360,
                                       attractors = list(source_density = 1),
                                       permeability = Permeability(template,
                                                                   region),
                                       max_distance = 300000))
})

test_that("packs population array", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  population <- Population(region)
  expect_silent(dispersal <- Dispersal(region, population_model = population))
  n <- c(rep(TRUE, 4), rep(FALSE, 10))
  idx <- which(n)
  expected <- list(indices = idx,
                   original = as.matrix(n[idx]),
                   remaining = as.matrix(n[idx]),
                   relocated = rep(FALSE, length(n)))
  expect_equal(dispersal$pack(n), expected)
})

test_that("unpacks population array", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  population <- Population(region, type = "unstructured")
  expect_silent(dispersal <- Dispersal(region, population_model = population))
  n <- c((4:1)*10, rep(0, 10))
  expect_silent(n <- dispersal$pack(n))
  n$remaining[1:2,] <- n$remaining[1:2,] - c(8, 4)
  n$relocated[5:6] <- c(8, 4)
  expect_equal(dispersal$unpack(n), c(n$remaining[,1], n$relocated[5:14]))
})

