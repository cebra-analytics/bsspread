context("Controls")

test_that("initializes with region, population, and other parameters", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "template.tif"))
  region <- Region(template)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  set.seed(1234)
  manage_pr <- runif(region$get_locations())
  expect_error(controls <- Controls(region, population_model,
                                    control_type = "search_destroy"),
               paste("Management control effectiveness values is required for",
                     "'search and destroy' control."))
  expect_error(controls <- Controls(region, population_model,
                                    control_type = "search_destroy",
                                    manage_pr = 1:10),
               paste("The management control effectiveness parameter must be",
                     "a numeric vector with values for each location."))
  expect_error(controls <- Controls(region, population_model,
                                    control_type = "search_destroy",
                                    manage_pr = manage_pr + 0.5),
               paste("Management control effectiveness values should be",
                     "<= 0 and <= 1."))
  expect_error(
    controls <- Controls(region, population_model,
                         control_type = "search_destroy",
                         manage_pr = manage_pr,
                         control_cost = 1:5),
    paste("The control cost parameter must be a numeric vector with values",
          "for each location."))
  control_cost <- 2
  attr(control_cost, "unit") <- "$"
  expect_silent(
    controls <- Controls(region, population_model,
                         control_type = "search_destroy",
                         manage_pr = manage_pr,
                         control_cost = control_cost,
                         stages = 2:3, schedule = 4:6))
  expect_is(controls, "Controls")
  expect_s3_class(controls, "Actions")
  expect_named(controls,
               c(c("get_type", "get_id", "set_id", "get_label", "get_stages",
                   "get_schedule", "include_cost", "get_cost_label",
                   "get_cost_unit", "clear_attributes", "apply",
                   "get_manage_pr_type")))
  expect_equal(controls$get_type(), "control")
  expect_equal(controls$get_label(), "control_search_destroy")
  expect_equal(controls$get_manage_pr_type(), "individual")
  expect_equal(controls$get_stages(), 2:3)
  expect_equal(controls$get_schedule(), 4:6)
  expect_true(controls$include_cost())
  expect_equal(controls$get_cost_label(), "control_search_destroy_cost")
  expect_equal(controls$get_cost_unit(), "$")
  expect_silent(controls$set_id(1))
  expect_equal(controls$get_id(), 1)
  expect_equal(controls$get_label(), "1_control_search_destroy")
  expect_equal(controls$get_label(include_id = FALSE),
               "control_search_destroy")
  expect_equal(controls$get_cost_label(),
               "1_control_search_destroy_cost")
  expect_equal(controls$get_cost_label(include_id = FALSE),
               "control_search_destroy_cost")
  expect_error(
    controls <- Controls(region, population_model,
                         control_type = "growth"),
    paste("Control multiplier is required for growth, spread, and",
          "establishment control."))
  expect_error(
    controls <- Controls(region, population_model,
                         control_type = "growth",
                         control_mult = (0:10)/10),
    paste("Control multiplier should be a vector with a value 0-1 for",
          "each region location."))
  expect_error(
    controls <- Controls(region, population_model,
                         control_type = "growth",
                         control_mult = 2),
    paste("Control multiplier should be a vector with a value 0-1 for",
          "each region location."))
  expect_error(controls <- Controls(region, population_model,
                       control_type = "growth",
                       control_mult = 0.7,
                       exist_control = (0:10)/10),
               paste("Existing control should be a logical or numeric vector",
                     "with a value for each region location."))
  expect_error(
    controls <- Controls(region, population_model,
                         control_type = "growth",
                         control_mult = 0.7,
                         radius = -1),
    "The radius (m) parameter must be numeric and >= 0.", fixed = TRUE)
  expect_error(
    controls <- Controls(region, population_model,
                         control_type = "growth",
                         control_mult = 0.7,
                         apply_to = "dummy"),
    paste("Growth control 'apply to' attribute should be 'reproduction' or",
          "'survival'."))
  expect_silent(
    controls <- Controls(region, population_model,
                         control_type = "growth",
                         control_design = control_design,
                         radius = 1000,
                         control_mult = 0.7,
                         exist_control = (manage_pr > 0.5)*2,
                         control_cost = control_cost,
                         stages = 2:3,
                         apply_to = "survival",
                         schedule = 4:6))
  expect_equal(controls$get_label(), "control_growth")
  expect_true(controls$include_cost())
  expect_equal(controls$get_cost_label(), "control_growth_cost")
  expect_equal(controls$get_cost_unit(), "$")
  expect_silent(controls$set_id(2))
  expect_equal(controls$get_id(), 2)
  expect_equal(controls$get_label(), "2_control_growth")
  expect_equal(controls$get_label(include_id = FALSE), "control_growth")
  expect_equal(controls$get_cost_label(), "2_control_growth_cost")
  expect_equal(controls$get_cost_label(include_id = FALSE),
               "control_growth_cost")
})

test_that("applies stochastic search and destroy controls to population", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "template.tif"))
  region <- Region(template)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  initial_n <- rep(0, region$get_locations())
  initial_n[101:150] <- 11:60
  initializer <- Initializer(initial_n, region = region,
                             population_model = population_model)
  set.seed(1234)
  n <- initializer$initialize()
  set.seed(1234)
  manage_pr <- runif(region$get_locations())
  idx <- 1:region$get_locations()
  set.seed(1234)
  expected_controls <- array(0, dim(n))
  colnames(expected_controls) <- colnames(n)
  expected_controls[101:150, 2:3] <- stats::rbinom(100, size = n[101:150, 2:3],
                                                   manage_pr[101:150])
  expected_controlled <- (rowSums((n - expected_controls)[,2:3]) == 0 &
                            rowSums(n[,2:3]) > 0)
  attr(n, "attachment") <- "extra"
  expected_n <- n - expected_controls
  attr(expected_n, "control_search_destroy") <- expected_controls
  expected_costs <- 2*(manage_pr > 0)
  attr(expected_costs, "unit") <- "$"
  attr(expected_n, "control_search_destroy_cost") <- expected_costs
  control_cost <- 2
  attr(control_cost, "unit") <- "$"
  expect_silent(
    controls <- Controls(region, population_model,
                         control_type = "search_destroy",
                         manage_pr = manage_pr,
                         control_cost = control_cost,
                         stages = 2:3, schedule = 4:6))
  set.seed(1234)
  expect_silent(new_n <- controls$apply(n, 4))
  expect_equal(new_n[idx,], expected_n[idx,])
  expect_equal(attr(new_n, "undetected")[idx,],
               n[idx,] - expected_controls[idx,])
  expect_equal(as.logical(attr(new_n, "control_search_destroy")),
               expected_controlled)
  expect_equal(attr(attr(new_n, "control_search_destroy"), "number")[idx,],
               expected_controls[idx,])
  expect_equal(attr(new_n, "control_search_destroy_cost")[idx],
               expected_costs[idx])
  expect_equal(attr(new_n, "attachment"), "extra")
  attr(n, "attachment") <- NULL
  # duplicate with extra detections and removals
  attr(new_n, "undetected")[101:110, 2:3] <- 0
  rm_idx <- which(new_n[101:110,] > 1)
  new_n[101:110,][rm_idx] <- new_n[101:110,][rm_idx] - 1
  attr(new_n, "undetected")[101:110, 1] <- new_n[101:110, 1]
  n_undetected <- attr(new_n, "undetected")[,]
  n_apply <- list(detected = new_n[,] - n_undetected,
                  undetected = n_undetected)
  idx1 <- which(rowSums(new_n[,2:3]) > 0)
  set.seed(1234)
  expected_controls2 <- lapply(n_apply, function(a) {
    controlled <- expected_controls*0
    controlled[idx1, 2:3] <- stats::rbinom(length(idx1)*2, size = a[idx1, 2:3],
                                           manage_pr[idx1])
    return(controlled)
  })
  expected_controls_plus <-
    expected_controls2$detected + expected_controls2$undetected
  expected_controlled2 <-
    (rowSums((new_n - expected_controls_plus)[,2:3]) == 0 &
       rowSums(new_n[,2:3]) > 0)
  set.seed(1234)
  expect_silent(new_n2 <- controls$apply(new_n, 4))
  expect_equal(new_n2[idx,], new_n[idx,] - expected_controls_plus[idx,])
  expect_equal(attr(new_n2, "undetected")[idx,],
               (attr(new_n, "undetected")[idx,] -
                  expected_controls2$undetected[idx,]))
  expect_equal(as.logical(attr(new_n2, "control_search_destroy")),
               expected_controlled2)
  expect_equal(attr(attr(new_n2, "control_search_destroy"), "number")[idx,],
               expected_controls_plus[idx,])
  expect_equal(attr(new_n2, "control_search_destroy_cost")[idx],
               expected_costs[idx])
  expect_silent(controls$set_id(2))
  set.seed(1234)
  expect_silent(new_n <- controls$apply(n, 4))
  new_n[idx,] ; expected_n[idx,]
  expect_equal(attr(new_n, "undetected")[idx,],
               n[idx,] - expected_controls[idx,])
  expect_equal(as.logical(attr(new_n, "2_control_search_destroy")),
               expected_controlled)
  expect_equal(attr(attr(new_n, "2_control_search_destroy"), "number")[idx,],
               expected_controls[idx,])
  expect_equal(attr(new_n, "2_control_search_destroy_cost")[idx],
               expected_costs[idx])
  # population level effectiveness
  expect_silent(
    controls <- Controls(region, population_model,
                         control_type = "search_destroy",
                         manage_pr = manage_pr,
                         manage_pr_type = "population",
                         control_cost = control_cost,
                         stages = 2:3, schedule = 4:6))
  idx1 <- which(rowSums(n[,2:3]) > 0)
  controlled <- undetected <- destroyed <- rep(0, length(idx1))
  set.seed(1234)
  for (i in 1:1000) {
    new_n <- controls$apply(n, 4)
    controlled <- controlled + (rowSums(new_n[idx1, 2:3]) == 0)
    undetected <- undetected + (
      rowSums(attr(new_n, "undetected")[idx1, 2:3]) > 0)
    destroyed <- destroyed + attr(new_n, "control_search_destroy")[idx1]
  }
  expect_true(all(abs(controlled/1000 - manage_pr[idx1]) < 0.05))
  expect_true(all(abs(undetected/1000 - (1 - manage_pr[idx1])) < 0.05))
  expect_equal(destroyed, controlled)
  # unstructured
  population_model <- UnstructPopulation(region, growth = 1.2)
  n <- rowSums(n)
  set.seed(1234)
  expected_controls <- n*0
  expected_controls[101:150] <- stats::rbinom(50, size = n[101:150],
                                              manage_pr[101:150])
  expected_controlled <- n - expected_controls == 0 & n > 0
  expected_n <- n - expected_controls
  expected_costs <- 2*(manage_pr > 0)
  attr(expected_costs, "unit") <- "$"
  expect_silent(
    controls <- Controls(region, population_model,
                         control_type = "search_destroy",
                         manage_pr = manage_pr,
                         control_cost = control_cost,
                         schedule = 4:6))
  set.seed(1234)
  expect_silent(new_n <- controls$apply(n, 4))
  expect_equal(as.numeric(new_n), expected_n)
  expect_equal(attr(new_n, "undetected"), n - expected_controls)
  expect_equal(as.logical(attr(new_n, "control_search_destroy")),
               expected_controlled)
  expect_equal(attr(attr(new_n, "control_search_destroy"), "number"),
               expected_controls)
  expect_equal(attr(new_n, "control_search_destroy_cost"), expected_costs)
  # duplicate with extra detections and removals
  attr(new_n, "undetected")[101:110] <- 0
  rm_idx <- which(new_n[101:110] > 1)
  new_n[101:110][rm_idx] <- new_n[101:110][rm_idx] - 1
  n_undetected <- as.numeric(attr(new_n, "undetected"))
  n_apply <- list(detected = as.numeric(new_n) - n_undetected,
                  undetected = n_undetected)
  set.seed(1234)
  idx1 <- which(new_n > 0)
  expected_controls2 <- lapply(n_apply, function(a) {
    controlled <- expected_controls*0
    controlled[idx1] <- stats::rbinom(length(idx1), size = a[idx1],
                                      manage_pr[idx1])
    return(controlled)
  })
  expected_controls_plus <-
    expected_controls2$detected + expected_controls2$undetected
  expected_controlled2 <- new_n - expected_controls_plus == 0 & new_n > 0
  set.seed(1234)
  expect_silent(new_n2 <- controls$apply(new_n, 4))
  expect_equal(as.numeric(new_n2), as.numeric(new_n) - expected_controls_plus)
  expect_equal(attr(new_n2, "undetected"),
               attr(new_n, "undetected") - expected_controls2$undetected)
  expect_equal(as.logical(attr(new_n2, "control_search_destroy")),
               expected_controlled2)
  expect_equal(attr(attr(new_n2, "control_search_destroy"), "number"),
               expected_controls_plus)
  expect_equal(attr(new_n2, "control_search_destroy_cost")[idx],
               expected_costs[idx])
  # population level effectiveness
  expect_silent(
    controls <- Controls(region, population_model,
                         control_type = "search_destroy",
                         manage_pr = manage_pr,
                         manage_pr_type = "population",
                         control_cost = control_cost,
                         schedule = 4:6))
  idx1 <- which(n > 0)
  controlled <- undetected <- destroyed <- rep(0, length(idx1))
  set.seed(1234)
  for (i in 1:1000) {
    new_n <- controls$apply(n, 4)
    controlled <- controlled + (new_n[idx1] == 0)
    undetected <- undetected + (attr(new_n, "undetected")[idx1] > 0)
    destroyed <- destroyed + attr(new_n, "control_search_destroy")[idx1]
  }
  expect_true(all(abs(controlled/1000 - manage_pr[idx1]) < 0.05))
  expect_true(all(abs(undetected/1000 - (1 - manage_pr[idx1])) < 0.05))
  expect_equal(destroyed, controlled)
  # presence-only
  population_model <- PresencePopulation(region)
  n <- n > 0
  set.seed(1234)
  expected_controls <- n*0
  expected_controls[101:150] <- stats::rbinom(50, size = n[101:150],
                                              manage_pr[101:150])
  expected_n <- n - expected_controls
  expected_costs <- 2*(manage_pr > 0)
  attr(expected_costs, "unit") <- "$"
  expect_silent(
    controls <- Controls(region, population_model,
                         control_type = "search_destroy",
                         manage_pr = manage_pr,
                         control_cost = control_cost,
                         schedule = 4:6))
  set.seed(1234)
  expect_silent(new_n <- controls$apply(n, 4))
  expect_equal(as.numeric(new_n), expected_n)
  expect_equal(attr(new_n, "undetected"), as.logical(n - expected_controls))
  expect_equal(attr(new_n, "control_search_destroy"),
               as.logical(expected_controls))
  expect_equal(attr(new_n, "control_search_destroy_cost"), expected_costs)
  # duplicate with extra detections and removals
  attr(new_n, "undetected")[101:110] <- FALSE
  new_n[101] <- FALSE
  n_undetected <- as.numeric(attr(new_n, "undetected"))
  n_apply <- list(detected = as.numeric(new_n) - n_undetected,
                  undetected = n_undetected)
  set.seed(1234)
  idx1 <- which(new_n > 0)
  expected_controls2 <- lapply(n_apply, function(a) {
    controlled <- expected_controls*0
    controlled[idx1] <- stats::rbinom(length(idx1), size = a[idx1],
                                      manage_pr[idx1])
    return(controlled)
  })
  expected_controls_plus <-
    expected_controls2$detected + expected_controls2$undetected
  set.seed(1234)
  expect_silent(new_n2 <- controls$apply(new_n, 4))
  expect_equal(as.numeric(new_n2), as.numeric(new_n) - expected_controls_plus)
  expect_equal(attr(new_n2, "undetected"),
               attr(new_n, "undetected") & !expected_controls2$undetected)
  expect_equal(attr(new_n2, "control_search_destroy"),
               as.logical(expected_controls_plus))
  expect_equal(attr(new_n2, "control_search_destroy_cost")[idx],
               expected_costs[idx])
})

test_that("applies stochastic growth/spread/establishment controls", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "template.tif"))
  region <- Region(template)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  initial_n <- rep(0, region$get_locations())
  initial_n[101:150] <- 11:60
  initializer <- Initializer(initial_n, region = region,
                             population_model = population_model)
  set.seed(1234)
  n <- initializer$initialize()
  set.seed(1234)
  manage_pr <- runif(region$get_locations())
  expected_costs <- rep(2, region$get_locations())
  attr(expected_costs, "unit") <- "$"
  idx <- 1:region$get_locations() # idx <- 101:120; idx <- 101:150
  exist_control = +(manage_pr > 0.95)
  exist_mask <- rowSums(n[,2:3])*exist_control > 0
  detected <- n*0
  detected[,2:3] <- round((n[,2:3] > 8)*n[,2:3]*0.7)
  attr(n, "undetected") <- n - detected
  detected_mask <- +(rowSums(detected) > 0 | NA)
  detected_mask <- rowSums(n[,2:3])*terra::buffer(
    region$get_rast(detected_mask), width = 1500)[region$get_indices()][,1] > 0
  control_cost <- 2
  attr(control_cost, "unit") <- "$"
  expect_silent(
    controls <- Controls(region, population_model,
                         control_type = "growth",
                         control_mult = 0.7,
                         exist_control = exist_control,
                         control_cost = control_cost,
                         radius = 1500,
                         stages = 2:3,
                         apply_to = "survival",
                         schedule = 4:6))
  expect_silent(new_n <- controls$apply(n, 4))
  expect_equal(new_n[idx,], n[idx,])
  expect_equal(attr(new_n, "undetected")[idx,], (n - detected)[idx,])
  expect_equal(as.numeric(attr(new_n, "control_growth")),
               ((exist_mask | detected_mask)*0.7 +
                  !(exist_mask | detected_mask)))
  expect_equal(attr(attr(new_n, "control_growth"), "stages"), 2:3)
  expect_equal(attr(attr(new_n, "control_growth"), "apply_to"), "survival")
  expect_equal(attr(new_n, "control_growth_cost"),
               (expected_costs*(exist_control > 0 | detected_mask)))
  establish_mask <- +(rowSums(detected) > 0 | NA)
  establish_mask <- terra::buffer(region$get_rast(establish_mask),
                                  width = 3000)[region$get_indices()][,1] > 0
  expect_silent(
    controls <- Controls(region, population_model,
                         control_type = "establishment",
                         control_mult = 0.7,
                         exist_control = exist_control,
                         control_cost = control_cost,
                         radius = 3000,
                         stages = 2:3,
                         apply_to = "survival",
                         schedule = 4:6))
  expect_silent(new_n <- controls$apply(n, 4))
  expect_equal(attr(new_n, "control_establishment"),
               ((exist_mask | establish_mask)*0.7 +
                  !(exist_mask | establish_mask)))
  expect_equal(attr(new_n, "control_establishment_cost"),
               (expected_costs*(exist_control > 0 | establish_mask)))
  # additional detection and removals
  detected[,2:3] <- round((n[,2:3] > 6)*n[,2:3]*0.7)
  new_n[which(rowSums(detected) >= 8),] <- 0
  attr(new_n, "undetected") <- new_n[,] - detected
  exist_mask <- rowSums(new_n[,2:3])*exist_control > 0
  detected_mask <- +(rowSums(detected) > 0 | NA)
  detected_mask <- terra::buffer(region$get_rast(detected_mask),
                                 width = 3000)[region$get_indices()][,1] > 0
  expect_silent(new_n2 <- controls$apply(new_n, 4))
  expected_mult <- rep(1, region$get_locations())
  expected_mult[which(exist_mask | detected_mask)] <-
    expected_mult[which(exist_mask | detected_mask)]*0.7
  expect_equal(attr(new_n2, "control_establishment"), expected_mult)
  expect_equal(
    attr(new_n2, "control_establishment_cost"),
    expected_costs*(exist_control > 0 | detected_mask))
  expect_silent(controls$set_id(3))
  expect_silent(new_n <- controls$apply(n, 4))
  detected[,2:3] <- round((n[,2:3] > 8)*n[,2:3]*0.7)
  expect_equal(attr(new_n, "undetected")[idx,], (n - detected)[idx,])
  expect_equal(attr(new_n, "3_control_establishment"),
               ((exist_mask | establish_mask)*0.7 +
                  !(exist_mask | establish_mask)))
  expect_equal(attr(new_n, "3_control_establishment_cost"),
               expected_costs*(exist_control > 0 | establish_mask))
  # spatially-implicit
  region <- Region()
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  initializer <- Initializer(50, region = region,
                             population_model = population_model)
  set.seed(1234)
  n <- initializer$initialize()
  detected <- n*0
  detected[,2:3] <- round(n[,2:3]*0.7)
  attr(n, "undetected") <- n - detected
  expect_silent(
    controls <- Controls(region, population_model,
                         control_type = "growth",
                         control_mult = 0.7,
                         control_cost = control_cost,
                         radius = NULL,
                         stages = 2:3,
                         apply_to = "survival",
                         schedule = 4:6))
  expect_error(new_n <- controls$apply(n, 4),
               paste("Cannot calculate spatially implicit control cost",
                     "without area occupied."))
  attr(n, "spread_area") <- 300
  set.seed(1234)
  expect_silent(new_n <- controls$apply(n, 4))
  expect_equal(new_n[,], n[,])
  expect_equal(attr(new_n, "undetected")[,], (n - detected)[,])
  expect_equal(as.numeric(attr(new_n, "control_growth")), 0.7)
  expect_equal(attr(attr(new_n, "control_growth"), "stages"), 2:3)
  expect_equal(attr(attr(new_n, "control_growth"), "apply_to"), "survival")
  expected_cost <- 2*300
  attr(expected_cost,"unit") <- "$"
  expect_equal(attr(new_n, "control_growth_cost"), expected_cost)
})
