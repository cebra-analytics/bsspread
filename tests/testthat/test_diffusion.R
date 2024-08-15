context("Diffusion")

test_that("initializes with region, population model, and other parameters", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  expect_silent(diffusion <- Diffusion(region, population_model = population))
  expect_is(diffusion, "Diffusion")
  expect_s3_class(diffusion, "Dispersal")
  expect_error(diffusion <- Diffusion(region, population_model = population,
                                      diffusion_rate = -1),
               "The diffusion rate must be numeric and >= 0.")
  expect_error(diffusion <- Diffusion(region, population_model = population,
                                      diffusion_coeff = 0),
               "The diffusion coefficient must be numeric and > 0.")
  expect_silent(diffusion <- Diffusion(region, population_model = population,
                                       diffusion_rate = 2000,
                                       diffusion_coeff = 1000000,
                                       proportion = 1,
                                       density_dependent = TRUE,
                                       direction_function = function(x) x/360,
                                       attractors = list(Attractor(template,
                                                                   region)),
                                       permeability = Permeability(template,
                                                                   region)))
})

test_that("sets dispersal parameters and combined function appropriately", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  diffusion <- Diffusion(region, population_model = population)
  expect_silent(combined_function <-
                  get("combined_function",
                      envir = environment(diffusion$disperse)))
  expect_is(combined_function, "function")
  expect_equal(combined_function(list(c(1000, 1414, 2000)), c(0, 45, 90)),
               c(0, 0, 0))
  expect_silent(diffusion_rate <- get("diffusion_rate",
                                      envir = environment(combined_function)))
  expect_equal(diffusion_rate, 0)
  expect_silent(max_distance <- get("max_distance",
                                    envir = environment(diffusion$disperse)))
  expect_equal(max_distance, region$get_res()*2)
  expect_null(get("events", envir = environment(diffusion$disperse)))
  expect_false(get("distance_adjust", envir = environment(diffusion$disperse)))
  diffusion <- Diffusion(region, population_model = population,
                         diffusion_rate = 500)
  combined_function <- get("combined_function",
                           envir = environment(diffusion$disperse))
  expect_equal(get("diffusion_rate", envir = environment(combined_function)),
               500)
  expect_equal(get("max_distance", envir = environment(diffusion$disperse)),
               region$get_res()*2)
  expect_equal(combined_function(list(c(1000, 1414, 2000)), c(0, 45, 90)),
               c((500/1000)^2,
                 (1000 + 1000*cos(pi*45/180) - 1414)/1414*(500/1000)^2, 0))
  expect_silent(diffusion <- Diffusion(region, population_model = population,
                                       diffusion_rate = 2000))
  combined_function <- get("combined_function",
                           envir = environment(diffusion$disperse))
  expect_equal(get("diffusion_rate", envir = environment(combined_function)),
               2000)
  expect_equal(get("max_distance", envir = environment(diffusion$disperse)),
               2000 + region$get_res())
  expect_equal(combined_function(list(c(1414, 2000, 2828)), c(45, 90, 135)),
               c(1, 1, 0))
})

test_that("diffuses population in a raster grid region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[4320] <- TRUE
  diffusion <- Diffusion(region, population_model = population)
  expect_silent(n <- diffusion$pack(n))
  expect_silent(n <- diffusion$disperse(n))
  expect_silent(n <- diffusion$unpack(n))
  expect_equal(which(n), 4320)
  diffusion <- Diffusion(region, population_model = population,
                         diffusion_rate = 2000, proportion = 1)
  idx <- which(diffusion$unpack(diffusion$disperse(diffusion$pack(n))))
  coords <- terra::xyFromCell(template, region$get_indices()[c(4320, idx)])
  distances <- terra::distance(coords[1,,drop = FALSE], coords[-1,], lonlat = FALSE)
  expect_true(all(distances < 2000 + 1000))
  expect_true(max(distances) >= 2000)
})

test_that("diffuses population in a patch region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  region$calculate_paths(1)
  paths <- region$get_paths(1)
  region <- Region(locations)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[1] <- TRUE
  diffusion <- Diffusion(region, population_model = population,
                         diffusion_rate = 200000, proportion = 1)
  expect_silent(n <- diffusion$pack(n))
  expect_silent(n <- diffusion$disperse(n))
  expect_silent(n <- diffusion$unpack(n))
  expect_equal(which(n),
               c(1, paths$idx[["1"]][which(paths$distances[["1"]] <= 200000)]))
})

test_that("spatially implicit radial diffusion with presence-only", {
  TEST_DIRECTORY <- test_path("test_inputs")
  region <- Region()
  region$set_max_implicit_area(1e8)
  population <- PresencePopulation(region)
  n <- TRUE
  diffusion <- Diffusion(region, population_model = population,
                         diffusion_rate = 2000, proportion = 1)
  expect_silent(n <- diffusion$pack(n))
  expect_silent(n <- diffusion$disperse(n))
  expect_silent(n <- diffusion$unpack(n))
  expect_true(n)
  expect_equal(attr(n, "diffusion_radius"), 2000)
  expect_silent(n <- diffusion$pack(n))
  expect_silent(n <- diffusion$disperse(n))
  expect_silent(n <- diffusion$unpack(n))
  expect_equal(attr(n, "diffusion_radius"), 4000)
  expect_silent(n <- diffusion$pack(n))
  expect_silent(n <- diffusion$disperse(n))
  expect_silent(n <- diffusion$unpack(n))
  expect_equal(attr(n, "diffusion_radius"), sqrt(1e8/pi))
})

test_that("spatially implicit reaction diffusion unstructured", {
  TEST_DIRECTORY <- test_path("test_inputs")
  region <- Region()
  region$set_max_implicit_area(1e8)
  population <- UnstructPopulation(region, growth = 1.2)
  n <- 120
  diffusion <- Diffusion(region, population_model = population,
                         diffusion_rate = 2000, proportion = 1)
  expect_silent(n <- diffusion$pack(n))
  expect_error(n <- diffusion$disperse(n),
               paste("The initial population needs to be set as an attribute",
                     "for reaction-diffusion calculations."))
  attr(n$relocated, "initial_n") <- 100
  expect_error(n <- diffusion$disperse(n),
               paste("The current time step needs to be set as an attribute",
               "for reaction-diffusion calculations."))
  attr(n$relocated, "tm") <- 1
  expected_n <- n$relocated + 120
  expect_silent(n <- diffusion$disperse(n))
  expect_silent(n <- diffusion$unpack(n))
  attr(expected_n, "diffusion_radius") <- 2000
  expect_equal(n, expected_n)
  attr(n, "diffusion_radius") <- NULL
  diffusion_coeff <- 2000^2/(4*log(1.2))
  diffusion <- Diffusion(region, population_model = population,
                         diffusion_coeff = diffusion_coeff, proportion = 1)
  expect_silent(n <- diffusion$pack(n))
  expect_silent(n <- diffusion$disperse(n))
  expect_silent(n <- diffusion$unpack(n))
  expect_equal(attr(n, "diffusion_radius"),
               sqrt(4*diffusion_coeff*1*log(120/100)))
  attr(n, "tm") <- 2
  n <- n + 24
  expect_silent(n <- diffusion$pack(n))
  expect_silent(n <- diffusion$disperse(n))
  expect_silent(n <- diffusion$unpack(n))
  expect_equal(attr(n, "diffusion_radius"),
               sqrt(4*diffusion_coeff*2*log(144/100)))
  attr(n, "tm") <- 3
  n <- n + 26
  expect_silent(n <- diffusion$pack(n))
  expect_silent(n <- diffusion$disperse(n))
  expect_silent(n <- diffusion$unpack(n))
  expect_equal(attr(n, "diffusion_radius"), sqrt(1e8/pi))
})
