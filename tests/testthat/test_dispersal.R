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

test_that("disperses population in a raster grid region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[5922] <- TRUE
  dispersal <- Dispersal(region, population_model = population)
  n <- dispersal$pack(n)
  dispersal$unpack(dispersal$disperse(n)) # n
  dispersal <- Dispersal(region, population_model = population,
                         proportion = 1)
  expect_true(all(dispersal$unpack(dispersal$disperse(n))))
  dispersal <- Dispersal(region, population_model = population,
                         proportion = 0.5)
  expect_equal(round(sum(
    dispersal$unpack(dispersal$disperse(n)))/region$get_locations(), 1), 0.5)
  dispersal <- Dispersal(region, population_model = population,
                         events = 100)
  set.seed(1234); events <- stats::rpois(1, 100); set.seed(1234)
  expect_equal(sum(dispersal$unpack(dispersal$disperse(n))), events)
  dispersal <- Dispersal(region, population_model = population,
                         proportion = 1, max_distance = 10000)
  idx <- region$get_paths(5922, max_distance = 10000)$idx[["5922"]]$cell
  expect_equal(sum(dispersal$unpack(dispersal$disperse(n))), length(idx) + 1)
  expect_true(all(dispersal$unpack(dispersal$disperse(n))[c(5922, idx)]))
})

test_that("disperses grid population with distance and direction functions", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[5922] <- TRUE
  expect_silent(dispersal <-
                  Dispersal(region, population_model = population,
                            proportion = 1, distance_adjust = FALSE,
                            distance_function = function(d) +(d <= 10000),
                            direction_function = function(d) +(d <= 180)))
  n <- dispersal$pack(n)
  region$calculate_paths(5922)
  idx <- region$get_paths(5922, max_distance = 10000)$idx[["5922"]]$cell
  directions <-
    region$get_paths(5922, directions = TRUE,
                     max_distance = 10000)$directions[["5922"]]$cell
  expect_equal(sum(dispersal$unpack(dispersal$disperse(n))),
               length(which(directions <= 180)) + 1)
  expect_true(all(
    dispersal$unpack(dispersal$disperse(n))[idx[which(directions <= 180)]]))
})

test_that("disperses grid population with attractors", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[5922] <- TRUE
  attractor_vect <- rep(0, region$get_locations())
  attractor_vect[1:5922] <- 1
  attractor <- Attractor(attractor_vect, region, type = "destination")
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       proportion = 1, max_distance = 10000,
                                       attractors = list(attractor)))
  n <- dispersal$pack(n)
  region$calculate_paths(5922)
  idx <- region$get_paths(5922, max_distance = 10000)$idx[["5922"]]$cell
  expect_equal(sum(dispersal$unpack(dispersal$disperse(n))),
               length(which(idx < 5922)) + 1)
  expect_true(all(
    dispersal$unpack(dispersal$disperse(n))[idx[which(idx < 5922)]]))
})

test_that("disperses grid population with permeability", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  region$configure_paths(max_distance = 30000)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[5922] <- TRUE
  perm_rast <- region$get_template()
  perm_rast[region$get_indices()[1:5922]] <- 0.5
  permeability <- Permeability(perm_rast, region)
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       proportion = 1, max_distance = 20000,
                                       permeability = permeability))
  n <- dispersal$pack(n)
  region$calculate_paths(5922)
  paths <- region$get_paths(5922, perm_id = 1)
  perm_dist <- paths$perm_dist[["5922"]]$cell
  idx <- paths$idx[["5922"]]$cell[which(perm_dist <= 20000)]
  distances <- paths$distances[["5922"]]$cell[which(perm_dist <= 20000)]
  expect_true(all(distances <= 20000*0.5))
  expect_true(all(which(dispersal$unpack(dispersal$disperse(n))) <= 5922))
  expect_equal(sum(dispersal$unpack(dispersal$disperse(n))), length(idx) + 1)
  expect_true(all(dispersal$unpack(dispersal$disperse(n))[idx]))
})

test_that("next", {
  TEST_DIRECTORY <- test_path("test_inputs")
  expect_true(TRUE)
})
