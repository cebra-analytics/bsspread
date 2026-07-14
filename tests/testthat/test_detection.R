context("Detection")

test_that("initializes with region, population, and sensitivity", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  sensitivity <- template[region$get_indices()][,1]
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  expect_error(detection <- Detection(region, population_model, 1:10,
                                      stages = 2:3),
               paste("The surveillance sensitivity parameter must be a",
                     "numeric vector with values for each location."))
  expect_error(detection <- Detection(region, population_model, sensitivity + 0.4,
                                      stages = 2:3),
               "Surveillance sensitivity values should be <= 0 and <= 1.")
  expect_error(detection <- Detection(
    region, population_model, sensitivity,
    sensitivity_threshold = 0),
    paste("The sensitivity population size threshold parameter should be a",
          "numeric value > 0."))
  expect_message(detection <- Detection(
    region, population_model, sensitivity,
    sensitivity_type = "individual",
    sensitivity_threshold = 5),
    paste("Ignoring the sensitivity population size threshold value, it is",
          "not used for a individual level sensitivity type."))
  expect_error(detection <- Detection(
    region, population_model, sensitivity,
    sensitivity_type = "population"),
    paste("A sensitivity population size threshold is required for a",
          "population level sensitivity type."))
  expect_message(detection <- Detection(
    region, population_model, sensitivity,
    sensitivity_type = "presence",
    sensitivity_threshold = 5),
    paste("The sensitivity population size threshold value of 5 has been set",
          "to 1 for a presence level sensitivity type."))
  expect_error(detection <- Detection(
    region, population_model, sensitivity, surv_cost = 1:5, stages = 2:3,
    schedule = 4:6),
    paste("The surveillance cost parameter must be a numeric vector with",
          "values for each location."))
  surv_cost <- 1
  attr(surv_cost, "unit") <- "$"
  expect_silent(detection <- Detection(
    region, population_model, sensitivity,
    sensitivity_type = "population",
    sensitivity_threshold = 5,
    surv_cost = surv_cost,
    stages = 2:3,
    schedule = 4:6))
  expect_is(detection, "Detection")
  expect_s3_class(detection, "Actions")
  expect_named(detection,
               c(c("get_type", "get_id", "set_id", "get_label", "get_stages",
                   "get_schedule", "include_cost", "get_cost_label",
                   "get_cost_unit", "get_attributes", "clear_attributes",
                   "apply", "get_sensitivity", "get_sensitivity_type")))
  expect_equal(detection$get_type(), "detection")
  expect_equal(detection$get_label(), "detected")
  expect_equal(detection$get_sensitivity(), sensitivity)
  expect_equal(detection$get_sensitivity_type(), "population")
  expect_equal(detection$get_stages(), 2:3)
  expect_equal(detection$get_schedule(), 4:6)
  expect_true(detection$include_cost())
  expect_equal(detection$get_cost_label(), "surv_cost")
  expect_equal(detection$get_cost_unit(), "$")
  expect_silent(detection$set_id(1))
  expect_equal(detection$get_id(), 1)
  expect_equal(detection$get_label(), "1_detected")
  expect_equal(detection$get_label(include_id = FALSE), "detected")
  expect_equal(detection$get_cost_label(), "1_surv_cost")
  expect_equal(detection$get_cost_label(include_id = FALSE),
               "surv_cost")
})

test_that("applies stochastic detection to invasive population", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  idx <- 5916:5922
  region <- Region(template)
  template[region$get_indices()][idx,] <- c(rep(0.5, 4), 0.5, 0.75, 1)
  idx <- idx[5:7]
  sensitivity <- template[region$get_indices()][,1]
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  initial_n <- rep(0, region$get_locations())
  initial_n[idx] <- (10:12)*5
  initializer <- Initializer(initial_n, region = region,
                             population_model = population_model)
  set.seed(1212)
  n <- initializer$initialize()
  no_detections <- expected_detections <- rep(FALSE, region$get_locations())
  set.seed(1234)
  expected_detected <- array(c(rep(0, 3),
                               stats::rbinom(6, size = n[idx, 2:3],
                                             c(0.5, 0.75, 1))), c(3, 3))
  colnames(expected_detected) <- colnames(n)
  expected_detections[idx] <- rowSums(expected_detected) > 0
  surv_cost <- 2
  attr(surv_cost, "unit") <- "$"
  expect_silent(detection <- Detection(
    region, population_model, sensitivity, surv_cost = surv_cost,
    stages = 2:3, schedule = 4:6))
  expect_silent(new_n <- detection$apply(n, 2))
  expect_equal(as.logical(attr(new_n, "detected")), no_detections)
  expect_equal(attr(attr(new_n, "detected"), "number")[idx,],
               expected_detected*0)
  expect_equal(attr(new_n, "undetected"), n)
  expect_equal(as.numeric(attr(new_n, "surv_cost")),
               rep(0, region$get_locations()))
  expect_equal(attr(attr(new_n, "surv_cost"), "unit"), "$")
  set.seed(1234)
  expect_silent(new_n <- detection$apply(n, 4))
  expect_equal(as.logical(attr(new_n, "detected")), expected_detections)
  expect_equal(attr(attr(new_n, "detected"), "number")[idx,],
               expected_detected)
  expect_equal(attr(new_n, "undetected"),
               n - attr(attr(new_n, "detected"), "number"))
  expect_equal(as.numeric(attr(new_n, "surv_cost")), 2*(sensitivity > 0))
  expect_equal(attr(attr(new_n, "surv_cost"), "unit"), "$")
  set.seed(1234)
  expect_silent(new_n2 <- detection$apply(new_n, 4))
  expect_true(all(attr(attr(new_n2, "detected"), "number") <=
                    attr(attr(new_n, "detected"), "number")))
  expect_equal(attr(new_n2, "undetected"),
               (n - attr(attr(new_n, "detected"), "number") -
                  attr(attr(new_n2, "detected"), "number")))
  expect_equal(as.numeric(attr(new_n2, "surv_cost")), 2*(sensitivity > 0))
  expect_equal(attr(attr(new_n2, "surv_cost"), "unit"), "$")
  expect_equal(detection$get_attributes(new_n2),
               list(detected = attr(new_n2, "detected"),
                    undetected = attr(new_n2, "undetected"),
                    surv_cost = attr(new_n2, "surv_cost")))
  expect_silent(new_n2 <- detection$clear_attributes(new_n2))
  expect_null(attr(new_n2, "detected"))
  expect_null(attr(new_n2, "undetected"))
  expect_null(attr(new_n2, "surv_cost"))
  expect_equal(new_n2, n)
  expect_silent(detection$set_id(1))
  set.seed(1234)
  expect_silent(new_n <- detection$apply(n, 4))
  expect_equal(as.logical(attr(new_n, "1_detected")), expected_detections)
  expect_equal(attr(attr(new_n, "1_detected"), "number")[idx,],
               expected_detected)
  expect_equal(attr(new_n, "undetected"),
               n - attr(attr(new_n, "1_detected"), "number"))
  expect_equal(as.numeric(attr(new_n, "1_surv_cost")), 2*(sensitivity > 0))
  expect_equal(attr(attr(new_n, "1_surv_cost"), "unit"), "$")
  # population level sensitivity
  expect_silent(detection <- Detection(
    region, population_model, sensitivity,
    sensitivity_type = "population",
    sensitivity_threshold = 20,
    stages = 2:3, schedule = 4:6))
  adj_sens <- c(0.5, 0.75, 1)*pmin(rowSums(n[idx, 2:3]), 20)/20
  detected <- rep(0, length(idx))
  set.seed(1234)
  for (i in 1:1000) {
    new_n <- detection$apply(n, 4)
    detected <- detected + attr(new_n, "detected")[idx]
  }
  expect_true(all(abs(detected/1000 - adj_sens) < 0.05))
  # presence/absence level sensitivity
  expect_silent(detection <- Detection(
    region, population_model, sensitivity,
    sensitivity_type = "presence",
    stages = 2:3, schedule = 4:6))
  detected <- rep(0, length(idx))
  set.seed(1234)
  for (i in 1:1000) {
    new_n <- detection$apply(n, 4)
    detected <- detected + attr(new_n, "detected")[idx]
  }
  expect_true(all(abs(detected/1000 - c(0.5, 0.75, 1)) < 0.05))
  expect_equal(detected[3], 1000)
  # unstructured population
  population_model <- UnstructPopulation(region, growth = 1.2)
  n <- rowSums(n)
  set.seed(1234)
  expected_detected <- stats::rbinom(3, size = n[idx], c(0.5, 0.75, 1))
  expect_silent(
    detection <- Detection(region, population_model, sensitivity,
                           surv_cost = surv_cost, schedule = 4:6))
  set.seed(1234)
  expect_silent(new_n <- detection$apply(n, 4))
  expect_equal(as.logical(attr(new_n, "detected")), expected_detections)
  expect_equal(attr(attr(new_n, "detected"), "number")[idx], expected_detected)
  expect_equal(attr(new_n, "undetected"),
               n - attr(attr(new_n, "detected"), "number"))
  expect_equal(as.numeric(attr(new_n, "surv_cost")), 2*(sensitivity > 0))
  expect_equal(attr(attr(new_n, "surv_cost"), "unit"), "$")
  set.seed(1234)
  expect_silent(new_n2 <- detection$apply(new_n, 4))
  expect_true(all(attr(attr(new_n2, "detected"), "number") <=
                    attr(attr(new_n, "detected"), "number")))
  expect_equal(attr(new_n2, "undetected"),
               (n - attr(attr(new_n, "detected"), "number") -
                  attr(attr(new_n2, "detected"), "number")))
  expect_equal(as.numeric(attr(new_n2, "surv_cost")), 2*(sensitivity > 0))
  expect_equal(attr(attr(new_n2, "surv_cost"), "unit"), "$")
  # population level sensitivity
  expect_silent(detection <- Detection(
    region, population_model, sensitivity,
    sensitivity_type = "population",
    sensitivity_threshold = 20,
    schedule = 4:6))
  n[idx] <- 75 - n[idx]
  adj_sens <- c(0.5, 0.75, 1)*pmin(n[idx], 20)/20
  detected <- rep(0, length(idx))
  set.seed(1234)
  for (i in 1:1000) {
    new_n <- detection$apply(n, 4)
    detected <- detected + attr(new_n, "detected")[idx]
  }
  expect_true(all(abs(detected/1000 - adj_sens) < 0.05))
  # presence-only population
  population_model <- PresencePopulation(region)
  n <- n > 0
  set.seed(1234)
  expected_detected <- stats::rbinom(3, size = n[idx], c(0.5, 0.75, 1))
  expect_silent(
    detection <- Detection(region, population_model, sensitivity,
                           surv_cost = surv_cost, schedule = 4:6))
  set.seed(1234)
  expect_silent(new_n <- detection$apply(n, 4))
  expect_equal(attr(new_n, "detected")[idx], as.logical(expected_detected))
  expect_equal(attr(new_n, "undetected"), n & !attr(new_n, "detected"))
  expect_equal(as.numeric(attr(new_n, "surv_cost")), 2*(sensitivity > 0))
  expect_equal(attr(attr(new_n, "surv_cost"), "unit"), "$")
  set.seed(1234)
  expect_silent(new_n2 <- detection$apply(new_n, 4))
  expect_true(sum(attr(new_n2, "detected")) <= sum(attr(new_n, "undetected")))
  expect_equal((attr(new_n, "detected") | attr(new_n2, "detected") |
                  attr(new_n2, "undetected")), n)
  expect_equal(as.numeric(attr(new_n2, "surv_cost")), 2*(sensitivity > 0))
  expect_equal(attr(attr(new_n2, "surv_cost"), "unit"), "$")
})
