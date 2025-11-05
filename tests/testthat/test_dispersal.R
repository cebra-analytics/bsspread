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
                                      proportion = -1),
               "The proportion parameter must be numeric >= 0.")
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      events = 0:2),
               paste("The events parameter should be a single value or a",
                     "matrix with a single row or a row for each region",
                     "location."))
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      density_dependent = 0),
               "The density dependent parameter must be logical.")
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      distance_function = 0),
               "The distance function must be a function.")
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      direction_function = 0),
               "The direction function must be a function.")
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      combined_function = 0),
               "The combined function must be a function.")
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      distance_adjust = "y"),
               "The distance adjust parameter must be logical.")
  expect_error(dispersal <- Dispersal(region, population_model = population,
                                      attractors = 0),
               paste("Attractors must be a list containing zero or more",
                     "'Attractor' or inherited class objects."),
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
                                       density_dependent = TRUE,
                                       distance_function = function(x) 1/x,
                                       direction_function = function(x) x/360,
                                       combined_function = function(x, y)
                                         1/x[[1]]*y/360,
                                       distance_adjust = FALSE,
                                       attractors = list(Attractor(template,
                                                                   region)),
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
  population <- UnstructPopulation(region)
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
  dispersal <- Dispersal(region, population_model = population)
  n <- rep(FALSE, region$get_locations())
  expect_equal(dispersal$unpack(
    dispersal$disperse(dispersal$pack(n), tm = 1)), n)
  n[5922] <- TRUE
  n <- dispersal$pack(n)
  expect_equal(dispersal$unpack(dispersal$disperse(n, tm = 1)),
               dispersal$unpack(n))
  dispersal <- Dispersal(region, population_model = population,
                         proportion = 1)
  expect_true(all(dispersal$unpack(dispersal$disperse(n, tm = 1))))
  dispersal <- Dispersal(region, population_model = population,
                         proportion = 0.5)
  expect_equal(round(sum(
    dispersal$unpack(dispersal$disperse(n, tm = 1)))/
      region$get_locations(), 1), 0.5)
  dispersal <- Dispersal(region, population_model = population,
                         events = 100)
  set.seed(1234); new_locs <- stats::rpois(1, 100); set.seed(1234)
  expect_true(sum(
    dispersal$unpack(dispersal$disperse(n, tm = 1))) <= new_locs + 1)
  # spread and establishment control
  n <- dispersal$unpack(n)
  attr(n, "control_spread") <- 0.5
  n <- dispersal$pack(n)
  dispersal <- Dispersal(region, population_model = population,
                         proportion = 0.5)
  expect_equal(round(sum(dispersal$unpack(dispersal$disperse(n, tm = 1)))/
                       region$get_locations(), 2), 0.25)
  dispersal <- Dispersal(region, population_model = population,
                         events = 100)
  set.seed(1234); new_locs <- stats::rpois(1, 50); set.seed(1234)
  expect_true(sum(
    dispersal$unpack(dispersal$disperse(n, tm = 1))) <= new_locs + 1)
  attr(n$relocated, "control_spread") <- NULL
  establish_pr <- template[region$get_indices()][,1]
  population <- Population(region, establish_pr = establish_pr)
  dispersal <- Dispersal(region, population_model = population,
                         events = 1000)
  set.seed(1234)
  no_estab_control <- sum(dispersal$unpack(dispersal$disperse(n, tm = 1)))
  attr(n$relocated, "control_establishment") <- 0.5
  set.seed(1234)
  expect_silent(n2 <- dispersal$unpack(dispersal$disperse(n, tm = 2)))
  expect_equal(round(sum(n2)/no_estab_control, 1), 0.5)
  expect_equal(attr(n2, "control_establishment"), 0.5)
  attr(n$relocated, "control_establishment") <- NULL
  population <- Population(region)
  # spatio-temporal
  prop_matrix <- matrix(rep(0, region$get_locations()),
                        nrow = region$get_locations(), ncol = 3)
  prop_matrix[5922,] <- c(1, 0.5, 0.2)
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       proportion = prop_matrix))
  expect_equal(round(sum(dispersal$unpack(dispersal$disperse(n, tm = 1)))/
                       region$get_locations(), 1), 1)
  expect_equal(round(sum(dispersal$unpack(dispersal$disperse(n, tm = 2)))/
                       region$get_locations(), 1), 0.5)
  expect_equal(round(sum(dispersal$unpack(dispersal$disperse(n, tm = 3)))/
                       region$get_locations(), 1), 0.2)
  events_matrix <- matrix(rep(0, region$get_locations()),
                          nrow = region$get_locations(), ncol = 3)
  events_matrix[5922,] <- c(100, 50, 20)
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       events = events_matrix))
  set.seed(1234); new_locs <- stats::rpois(1, 100); set.seed(1234)
  expect_true(sum(
    dispersal$unpack(dispersal$disperse(n, tm = 1))) <= new_locs + 1)
  set.seed(1234); new_locs <- stats::rpois(1, 50); set.seed(1234)
  expect_true(sum(
    dispersal$unpack(dispersal$disperse(n, tm = 2))) <= new_locs + 1)
  set.seed(1234); new_locs <- stats::rpois(1, 20); set.seed(1234)
  expect_true(sum(
    dispersal$unpack(dispersal$disperse(n, tm = 3))) <= new_locs + 1)
  dispersal <- Dispersal(region, population_model = population,
                         proportion = 1, max_distance = 10000)
  idx <- region$get_paths(5922, max_distance = 10000)$idx[["5922"]]$cell
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(sum(new_n), length(idx) + 1)
  expect_true(all(new_n[c(5922, idx)]))
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
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(sum(new_n), length(which(directions <= 180)) + 1)
  expect_true(all(new_n[idx[which(directions <= 180)]]))
})

test_that("disperses grid population with combined function", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[5922] <- TRUE
  expect_silent(dispersal <-
                  Dispersal(region, population_model = population,
                            proportion = 1, distance_adjust = FALSE,
                            combined_function = function(x, d) {
                              +(x[[1]] <= 10000 & d <= 180)}))
  n <- dispersal$pack(n)
  region$calculate_paths(5922)
  idx <- region$get_paths(5922, max_distance = 10000)$idx[["5922"]]$cell
  directions <-
    region$get_paths(5922, directions = TRUE,
                     max_distance = 10000)$directions[["5922"]]$cell
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(sum(new_n), length(which(directions <= 180)) + 1)
  expect_true(all(new_n[idx[which(directions <= 180)]]))
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
  attractor <- Attractor(attractor_vect, region, is_dynamic = TRUE)
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       proportion = 1, max_distance = 10000,
                                       attractors = list(attractor)))
  n <- dispersal$pack(n)
  region$calculate_paths(5922)
  idx <- region$get_paths(5922, max_distance = 10000)$idx[["5922"]]$cell
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(sum(new_n), length(which(idx < 5922)) + 1)
  expect_true(all(new_n[idx[which(idx < 5922)]]))
  # dynamically linked to impacts
  mult <- rep(0, region$get_locations())
  idx <- region$get_paths(5922, max_distance = 5000)$idx[["5922"]]$cell
  mult[idx] <- 1
  n <- dispersal$unpack(n)
  attr(n, "dynamic_mult") <- list(NULL, list(mult))
  attr(attr(n, "dynamic_mult")[[2]], "links") <- "attractors"
  n <- dispersal$pack(n)
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(sum(new_n), length(which(idx < 5922)) + 1)
  expect_true(all(new_n[idx[which(idx < 5922)]]))
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
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_true(all(which(new_n) <= 5922))
  expect_equal(sum(new_n), length(idx) + 1)
  expect_true(all(new_n[idx]))
})

test_that("disperses population in a two-tier raster grid region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  region$configure_paths(max_distance = 80000)
  region$set_aggr(aggr_factor = 5, inner_radius = 10000)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[5922] <- TRUE
  perm_rast <- region$get_template()
  perm_rast[region$get_indices()[1:5922]] <- 0.5
  permeability <- Permeability(perm_rast, region)
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       proportion = 1, max_distance = 60000,
                                       permeability = permeability))
  n <- dispersal$pack(n)
  region$calculate_paths(5922)
  paths <- region$get_paths(5922, perm_id = 1)
  perm_dist <- paths$perm_dist[["5922"]]$aggr
  distances <- paths$distances[["5922"]]$aggr[which(perm_dist <= 60000)]
  expect_true(all(distances <= 60000*0.5))
  idx <- c(paths$idx$`5922`$cell,
           region$get_aggr()$get_cells(
             paths$idx$`5922`$aggr[which(perm_dist <= 60000)]))
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(sort(which(new_n)), sort(c(5922, idx)))
  expect_true(all(which(new_n) <= 5922))
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                         events = 400, max_distance = 60000,
                         permeability = permeability))
  expect_true(
    all(which(dispersal$unpack(
      dispersal$disperse(n, tm = 1))) %in% c(5922, idx)))
})

test_that("disperses unstructured population in a two-tier grid region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  region$set_aggr(aggr_factor = 5, inner_radius = 10000)
  establish_pr = template[region$get_indices()][,1]
  capacity_vect <- establish_pr*0 + 2000
  population <- UnstructPopulation(region, establish_pr = establish_pr,
                                   capacity = capacity_vect)
  n <- rep(0, region$get_locations())
  n[5922] <- 2000
  expect_silent(dispersal <- Dispersal(
    region, population_model = population, proportion = 1,
    density_dependent = TRUE, max_distance = 30000))
  n <- dispersal$pack(n)
  region$calculate_paths(5922)
  paths <- region$get_paths(5922)
  idx <- c(paths$idx$`5922`$cell,
           region$get_aggr()$get_cells(paths$idx$`5922`$aggr))
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_true(sum(new_n) < 2000)
  expect_true(new_n[5922] == 0)
  expect_true(length(which(new_n > 0)) <= sum(new_n))
  expect_true(all(which(new_n > 0) %in% idx))
  n$original[] <- n$remaining[] <- 1000
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_true(sum(new_n[-5922]) <= (1000 - new_n[5922]))
  expect_true(round(new_n[5922]/1000, 1) == 0.5)
  expect_true(all(which(new_n > 0) %in% c(5922, idx)))
  n$original[] <- n$remaining[] <- 2000
  proportion_vect <- n$relocated
  proportion_vect[5922] <- 0.7
  events_vect <- n$relocated
  events_vect[5922] <- 100
  expect_silent(dispersal <- Dispersal(
    region, population_model = population, proportion = proportion_vect,
    events = events_vect, max_distance = 30000))
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_true(sum(new_n[-5922]) <= (2000 - new_n[5922]))
  expect_true(round(new_n[5922]/2000, 1) == 0.3)
  expect_true(all(which(new_n > 0) %in% c(5922, idx)))
  # dynamically linked suitability to impacts
  expect_silent(dispersal <- Dispersal(
    region, population_model = population, proportion = 1,
    density_dependent = TRUE, max_distance = 30000))
  mult <- rep(0, region$get_locations())
  mult[1:5922] <- 1
  n <- dispersal$unpack(n)
  attr(n, "dynamic_mult") <- list(NULL, list(mult))
  attr(attr(n, "dynamic_mult")[[2]], "links") <- "suitability"
  n <- dispersal$pack(n)
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_true(sum(new_n) < 2000)
  expect_true(new_n[5922] == 0)
  expect_true(length(which(new_n > 0)) <= sum(new_n))
  expect_true(all(which(new_n > 0) %in% idx[which(idx <= 5922)]))
})

test_that("disperses staged population in a two-tier raster grid region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  region$set_aggr(aggr_factor = 5, inner_radius = 10000)
  establish_pr = template[region$get_indices()][,1]
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population <- StagedPopulation(region,
                                 growth = stage_matrix,
                                 establish_pr = establish_pr)
  n <- rep(0, region$get_locations())
  n[5922] <- 2000
  n <- population$make(initial = n)
  expect_equal(sum(n[5922,]), 2000)
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       dispersal_stages = 1:2,
                                       proportion = 1, max_distance = 30000))
  expect_silent(n <- dispersal$pack(n))
  region$calculate_paths(5922)
  paths <- region$get_paths(5922)
  idx <- c(paths$idx$`5922`$cell,
           region$get_aggr()$get_cells(paths$idx$`5922`$aggr))
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_true(all(new_n[5922,] == c(0, 0, n$original[3])))
  expect_true(all(colSums(new_n)[1:2] < n$original[,1:2]))
  expect_true(all(which(rowSums(new_n) > 0) %in% c(5922, idx)))
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       dispersal_stages = 1:2,
                                       proportion = 0.7, events = 100,
                                       max_distance = 30000))
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_true(all(new_n[5922, 3] == n$original[3]))
  expect_true(all(colSums(new_n)[1:2] < n$original[,1:2]))
  expect_true(all(colSums(new_n[-5922,]) <= (n$original - new_n[5922,])))
  expect_true(all(round(new_n[5922,1:2]/n$original[,1:2], 1) == 0.3))
  expect_true(all(which(rowSums(new_n) > 0) %in% c(5922, idx)))
})

test_that("disperses population in a patch/network region", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[1] <- TRUE
  expect_silent(dispersal <- Dispersal(region, population_model = population))
  n <- dispersal$pack(n)
  expect_equal(dispersal$unpack(dispersal$disperse(n, tm = 1)),
               dispersal$unpack(n))
  dispersal <- Dispersal(region, population_model = population,
                         proportion = 1)
  expect_true(all(dispersal$unpack(dispersal$disperse(n, tm = 1))))
  dispersal <- Dispersal(region, population_model = population,
                         proportion = 0.5)
  set.seed(4321); new_locs <- sum(stats::rbinom(13, size = 1, prob = 0.5))
  set.seed(4321)
  expect_equal(sum(dispersal$unpack(dispersal$disperse(n, tm = 1))),
               new_locs + 1)
  dispersal <- Dispersal(region, population_model = population,
                         events = 8)
  set.seed(4321); new_locs <- stats::rpois(1, 8); set.seed(4321)
  expect_true(sum(dispersal$unpack(
    dispersal$disperse(n, tm = 1))) <= new_locs + 1)
  dispersal <- Dispersal(region, population_model = population,
                         proportion = 1, max_distance = 200000)
  idx <- region$get_paths(1, max_distance = 200000)$idx[["1"]]
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(sum(new_n), length(idx) + 1)
  expect_true(all(new_n[c(1, idx)]))
})

test_that("disperses in patch/network with distance and direction functions", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[1] <- TRUE
  expect_silent(dispersal <-
                  Dispersal(region, population_model = population,
                            proportion = 1, distance_adjust = FALSE,
                            distance_function = function(d) +(d <= 200000),
                            direction_function = function(d) +(d <= 150)))
  n <- dispersal$pack(n)
  region$calculate_paths(1)
  idx <- region$get_paths(1, max_distance = 200000)$idx[["1"]]
  directions <- region$get_paths(1, directions = TRUE,
                                 max_distance = 200000)$directions[["1"]]
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(sum(new_n), length(which(directions <= 150)) + 1)
  expect_true(all(new_n[idx[which(directions <= 150)]]))
})

test_that("disperses in patch/network with combined function", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[1] <- TRUE
  expect_silent(dispersal <-
                  Dispersal(region, population_model = population,
                            proportion = 1, distance_adjust = FALSE,
                            combined_function = function(x, d) {
                              +(x[[1]] <= 200000 & d <= 150)}))
  n <- dispersal$pack(n)
  region$calculate_paths(1)
  idx <- region$get_paths(1, max_distance = 200000)$idx[["1"]]
  directions <- region$get_paths(1, directions = TRUE,
                                 max_distance = 200000)$directions[["1"]]
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(sum(new_n), length(which(directions <= 150)) + 1)
  expect_true(all(new_n[idx[which(directions <= 150)]]))
})

test_that("disperses in patch/network with attractors", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[1] <- TRUE
  attractor_vect <- rep(0, region$get_locations())
  attractor_vect[1:8] <- 1
  attractor <- Attractor(attractor_vect, region, is_dynamic = TRUE)
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       proportion = 1, max_distance = 200000,
                                       attractors = list(attractor)))
  n <- dispersal$pack(n)
  region$calculate_paths(1)
  idx <- region$get_paths(1, max_distance = 200000)$idx[["1"]]
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(sum(new_n), length(which(idx <= 8)) + 1)
  expect_true(all(new_n[idx[which(idx <= 8)]]))
  # dynamically linked to impacts
  mult <- rep(1, region$get_locations())
  mult[3:4] <- 0
  n <- dispersal$unpack(n)
  attr(n, "dynamic_mult") <- list(NULL, list(mult))
  attr(attr(n, "dynamic_mult")[[2]], "links") <- "attractors"
  expected_n <- rep(0, region$get_locations())
  expected_n[c(1, idx)] <- 1
  expected_n <- expected_n*attractor_vect*mult
  n <- dispersal$pack(n)
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(as.numeric(new_n), expected_n)
})

test_that("disperses in patch/network with permeability", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  region$configure_paths(max_distance = 300000)
  population <- Population(region)
  n <- rep(FALSE, region$get_locations())
  n[1] <- TRUE
  perm_data <- matrix(c( 1,  2, 0.8,
                         1,  3, 0.7,
                         1,  4, 0.7,
                         1,  6, 0.6,
                         1, 10, 0.8,
                         2,  5, 0.6,
                         3, 12, 0.6,
                         3, 14, 0.5,
                         4, 13, 0.5,
                         5,  8, 0.6,
                         6,  9, 0.6,
                         7, 10, 0.5,
                         10, 11, 0.8), ncol = 3, byrow = TRUE)
  colnames(perm_data) <- c("i", "j", "weight")
  permeability <- Permeability(perm_data, region)
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       proportion = 1, max_distance = 200000,
                                       permeability = permeability))
  n <- dispersal$pack(n)
  region$calculate_paths(1)
  paths <- region$get_paths(1, perm_id = 1)
  perm_dist <- paths$perm_dist[["1"]]
  idx <- paths$idx[["1"]][which(perm_dist <= 200000)]
  distances <- paths$distances[["1"]][which(perm_dist <= 200000)]
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_equal(sum(new_n), length(idx) + 1)
  expect_true(all(new_n[idx]))
})

test_that("disperses unstructured population in patch/network", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  establish_pr = c(1, rep(0.8, 4), rep(0.6, 9))
  population <- UnstructPopulation(region, establish_pr = establish_pr)
  n <- rep(0, region$get_locations())
  n[1] <- 2000
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       proportion = 1, max_distance = 200000))
  n <- dispersal$pack(n)
  region$calculate_paths(1)
  paths <- region$get_paths(1)
  idx <- paths$idx$`1`
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_true(sum(new_n) < 2000)
  expect_true(new_n[1] == 0)
  expect_true(all(which(new_n > 0) %in% idx))
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       proportion = 0.7, events = 3,
                                       max_distance = 200000))
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_true(sum(new_n[-1]) <= (2000 - new_n[1]))
  expect_true(round(new_n[1]/2000, 1) == 0.3)
  expect_true(all(which(new_n > 0) %in% c(1, idx)))
  # dynamically linked suitability to impacts
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       proportion = 1, max_distance = 200000))
  mult <- rep(1, region$get_locations())
  mult[3:4] <- 0
  n <- dispersal$unpack(n)
  attr(n, "dynamic_mult") <- list(NULL, list(mult))
  attr(attr(n, "dynamic_mult")[[2]], "links") <- "suitability"
  n <- dispersal$pack(n)
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expected <- rep(0, region$get_locations())
  expected[idx] <- 1
  expected <- expected*mult
  expect_true(sum(new_n) < 2000)
  expect_true(new_n[1] == 0)
  expect_true(all(which(new_n > 0) %in% which(expected > 0)))
})

test_that("disperses staged population in patch/network", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  establish_pr = c(1, rep(0.8, 4), rep(0.6, 9))
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population <- StagedPopulation(region,
                                 growth = stage_matrix,
                                 establish_pr = establish_pr)
  n <- rep(0, region$get_locations())
  n[1] <- 2000
  n <- population$make(initial = n)
  sum(n[1,]) # 2000
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       dispersal_stages = 2:3,
                                       proportion = 1, max_distance = 200000))
  expect_silent(n <- dispersal$pack(n))
  region$calculate_paths(1)
  paths <- region$get_paths(1)
  idx <- paths$idx$`1`
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_true(all(new_n[1,] == c(n$original[1], 0, 0)))
  expect_true(all(colSums(new_n)[2:3] < n$original[,2:3]))
  expect_true(all(which(rowSums(new_n) > 0) %in% c(1, idx)))
  expect_silent(dispersal <- Dispersal(region, population_model = population,
                                       dispersal_stages = 2:3,
                                       proportion = 0.7, events = 100,
                                       max_distance = 200000))
  expect_silent(new_n <- dispersal$unpack(dispersal$disperse(n, tm = 1)))
  expect_true(all(new_n[1,1] == n$original[1]))
  expect_true(all(colSums(new_n)[2:3] < n$original[,2:3]))
  expect_true(all(colSums(new_n[-1,]) <= (n$original - new_n[1,])))
  expect_true(all(new_n[1, 2:3]/n$original[,2:3] > 0.2))
  expect_true(all(which(rowSums(new_n) > 0) %in% c(1, idx)))
})
