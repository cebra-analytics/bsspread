context("Actions")

test_that("initializes with region, population model, stages & schedule", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  expect_error(actions <- Actions(region, "population_model",
                                  stages = 2:3),
               paste("Population model must be a 'Population' or inherited",
                     "class object."))
  template2 <- template*1
  template2[1] <- NA
  expect_error(actions <-
                 Actions(Region(template2), population_model,
                         stages = 2:3),
               "Population model must be compatible with the region object.")
  expect_error(actions <- Actions(region, population_model,
                                  stages = 3:4),
               paste("Stages must be a vector of stage indices consistent",
                     "with the population model."))
  expect_error(actions <- Actions(region, population_model,
                                  schedule = "2"),
               paste("The schedule for applying actions should be a vector of",
                     "numeric simulation time steps."))
  expect_silent(actions <- Actions(region, population_model,
                                   stages = 2:3, schedule = 4:6))
  expect_is(actions, "Actions")
  expect_named(actions,
               c(c("get_type", "get_id", "set_id", "get_label", "get_stages",
                   "get_schedule", "include_cost", "get_cost_label",
                   "get_cost_unit", "get_attributes", "clear_attributes",
                   "apply")))
  expect_equal(actions$get_type(), "detection")
  expect_equal(actions$get_label(), "action")
  expect_equal(actions$get_stages(), 2:3)
  expect_equal(actions$get_schedule(), 4:6)
  expect_false(actions$include_cost())
  expect_equal(actions$get_cost_label(), "action_cost")
  expect_null(actions$get_cost_unit())
  expect_equal(actions$get_attributes(1:10), list())
  expect_equal(actions$clear_attributes(1:10), 1:10) # returns n
  expect_equal(actions$apply(1:10, 4), 1:10) # returns n
  expect_null(actions$get_id())
  expect_error(actions$set_id(1.2), "Actions id should be an integer >= 1.")
  expect_silent(actions$set_id(2))
  expect_equal(actions$get_id(), 2)
  expect_equal(actions$get_label(), "2_action")
  expect_equal(actions$get_label(include_id = FALSE), "action")
  expect_equal(actions$get_cost_label(), "2_action_cost")
  expect_equal(actions$get_cost_label(include_id = FALSE),
               "action_cost")
})
