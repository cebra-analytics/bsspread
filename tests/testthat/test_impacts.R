context("Impacts")

test_that("initializes with region, population model and impact details", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  population_model <-  PresencePopulation(region)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model_staged <- StagedPopulation(region, stage_matrix)
  expect_error(impacts_o <- Impacts("region", population_model),
               "Region model must be a 'Region' or inherited class object.")
  expect_error(impacts_o <- Impacts(region, "population_model"),
               paste("Population model must be a 'Population' or inherited",
                     "class object."))
  expect_error(impacts_o <- Impacts(region, population_model,
                                    impact_type = "area"),
               paste("Area-based impacts are only applicable for spatially",
                     "implicit regions."))
  expect_error(impacts_o <- Impacts(region, population_model,
                                    impact_type = "density"),
               paste("Density-based impacts are only available for",
                     "unstructured or stage-based populations with carrying",
                     "capacity."))
  expect_error(impacts_o <- Impacts(region, population_model,
                                    asset_name = 2),
               "Asset name should be a character string.")
  expect_error(impacts_o <- Impacts(region, population_model,
                                    asset_name = "crop",
                                    asset_value = 1:10),
               paste("Asset value should be either a spatial layer or numeric",
                     "vector compatible with the defined region."))
  expect_error(impacts_o <- Impacts(region, population_model,
                                    asset_name = "asset1",
                                    asset_value = 100*template,
                                    value_unit = 3),
               "Asset value unit should be a character string.")
  expect_error(impacts_o <- Impacts(region, population_model,
                                    asset_name = "asset1",
                                    asset_value = 100*template,
                                    loss_rate = -1),
               "Loss rate should be numeric, >= 0, and <= 1.")
  expect_error(impacts_o <- Impacts(region, population_model,
                                    asset_name = "asset1",
                                    asset_value = 100*template,
                                    loss_rate = 0.2,
                                    discount_rate = "0.3"),
               "Discount rate should be numeric, >= 0, and <= 1.")
  expect_error(impacts_o <- Impacts(region, population_model,
                                    valuation_type = "non-monetary",
                                    asset_name = "asset1",
                                    asset_value = 100*template,
                                    loss_rate = 0.2,
                                    discount_rate = 0.3),
               "Discount rate is only applicable for monetary value impacts.")
  expect_error(impacts_o <- Impacts(region, population_model_staged,
                                    asset_name = "asset1",
                                    asset_value = 100*template,
                                    loss_rate = 0.2,
                                    stages = 3:4),
               paste("Impact stages must be a vector of stage indices",
                     "consistent with the population model."))
  expect_error(impacts_o <- Impacts(region, population_model,
                                    asset_name = "asset1",
                                    asset_value = 100*template,
                                    loss_rate = 0.2,
                                    recovery_delay = -1),
               "Recover delay should a number >= 0.")
  expect_error(impacts_o <- Impacts(region, population_model,
                                    asset_name = "asset1",
                                    asset_value = 100*template,
                                    loss_rate = 0.2,
                                    dynamic_links = "suitability"),
               "Dynamic links are only applicable for dynamic value impacts.")
  expect_error(impacts_o <- Impacts(region, population_model,
                                    valuation_type = "dynamic",
                                    asset_name = "asset1",
                                    asset_value = 100*template,
                                    loss_rate = 0.2,
                                    dynamic_links = c("suitability", "growth")),
               paste("Dynamic links should be a vector of one or more of",
                     "'suitability', 'capacity', and/or 'attractors'."))
  expect_silent(impacts_o <- Impacts(region, population_model,
                                     asset_name = "asset1",
                                     asset_value = 100*template,
                                     loss_rate = 0.2))
  expect_is(impacts_o, "Impacts")

  expect_named(impacts_o,
               c("get_id", "set_id", "get_impact_type", "get_valuation_type",
                 "get_value_unit", "update_recovery_delay", "calculate"))
  expect_equal(impacts_o$get_id(), 1)
  expect_error(impacts_o$set_id(0), "Impacts id should be an integer >= 1.")
  expect_silent(impacts_o$set_id(2))
  expect_equal(impacts_o$get_id(), 2)
  expect_equal(impacts_o$get_impact_type(), "presence")
  expect_equal(impacts_o$get_valuation_type(), "monetary")
  expect_equal(impacts_o$get_value_unit(), "$")
})

test_that("calculates presence and density-based impacts", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  template[region$get_indices()][5922,] <- 0.25
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- StagedPopulation(region, stage_matrix)
  initial_n <- rep(0, region$get_locations())
  initial_n[5915:5922] <- 10:17
  initializer <- Initializer(initial_n, region = region,
                             population_model = population_model)
  n <- initializer$initialize()
  asset_value_1 <- 100*(template > 0.1 & template < 0.3)
  asset_value_2 <- 200*(template > 0.2 & template < 0.4)
  # presence-based impacts
  expect_silent(impacts_1 <- Impacts(region, population_model,
                                     asset_name = "impact1",
                                     asset_value = asset_value_1,
                                     loss_rate = 0.3,
                                     discount_rate = 0.05,
                                     stages = 2:3))
  expect_silent(impacts_1$set_id(1))
  expect_silent(impacts_2 <- Impacts(region, population_model,
                                     asset_name = "impact2",
                                     asset_value = asset_value_2,
                                     loss_rate = 0.4,
                                     discount_rate = 0.05,
                                     stages = 2:3))
  expect_silent(impacts_2$set_id(2))
  expected_impact_1 <- (pmin(rowSums(n[,2:3]), 1)*
                          asset_value_1[region$get_indices()][,1]*0.3/(1.05^4))
  expected_impact_2 <- (pmin(rowSums(n[,2:3]), 1)*
                          asset_value_2[region$get_indices()][,1]*0.4/(1.05^4))
  expect_silent(n <- impacts_1$calculate(n, 4))
  expect_silent(n <- impacts_2$calculate(n, 4))
  expect_true("impacts" %in% names(attributes(n)))
  calc_impacts <- attr(n, "impacts")
  expect_true(is.list(calc_impacts))
  idx <- 1:region$get_locations() # idx <- 5910:5925
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  # density-based impacts
  expect_error(impacts_1 <- Impacts(region, population_model,
                                    impact_type = "density",
                                    asset_name = "impact1",
                                    asset_value = asset_value,
                                    loss_rate = 0.3,
                                    stages = 2:3),
               paste("Density-based impacts are only available for",
                     "unstructured or stage-based populations with carrying",
                     "capacity."))
  population_model <- StagedPopulation(region, stage_matrix,
                                       capacity = (initial_n > 0)*5)
  expect_silent(impacts_1 <- Impacts(region, population_model,
                                     impact_type = "density",
                                     asset_name = "impact1",
                                     asset_value = asset_value_1,
                                     loss_rate = 0.3,
                                     stages = 2:3))
  expect_silent(impacts_1$set_id(1))
  expect_silent(impacts_2 <- Impacts(region, population_model,
                                     impact_type = "density",
                                     asset_name = "impact2",
                                     asset_value = asset_value_2,
                                     loss_rate = 0.4,
                                     stages = 2:3))
  expect_silent(impacts_2$set_id(2))
  n[5920:5922,] <- c(4, 9, 9, 4, 2, 3, 2, 0, 0)
  expected_impact_1 <- (pmin(rowSums(n[,2:3])/5, 1)*
                          asset_value_1[region$get_indices()][,1]*0.3)
  expected_impact_2 <- (pmin(rowSums(n[,2:3])/5, 1)*
                          asset_value_2[region$get_indices()][,1]*0.4)
  expect_silent(n <- impacts_1$calculate(n, 4))
  expect_silent(n <- impacts_2$calculate(n, 4))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
})

test_that("updates recovery delay to prolong presence-based impacts", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template*0)
  population_model <- UnstructPopulation(region, growth = 1.2)
  n <- rep(0, region$get_locations())
  n[5901:6000] <- round(runif(100, 1, 10))
  asset_value_1 <- 100*(template > 0.1 & template < 0.3)
  asset_value_2 <- 200*(template > 0.3)
  expect_silent(impacts_1 <- Impacts(region, population_model,
                                     asset_name = "impact1",
                                     asset_value = asset_value_1,
                                     loss_rate = 0.3,
                                     recovery_delay = 1))
  expect_silent(impacts_1$set_id(1))
  expect_silent(impacts_2 <- Impacts(region, population_model,
                                     asset_name = "impact2",
                                     asset_value = asset_value_2,
                                     loss_rate = 0.4,
                                     recovery_delay = 2))
  expect_silent(impacts_2$set_id(2))
  expected_impact_1 <- (n > 0)*asset_value_1[region$get_indices()][,1]*0.3
  expected_impact_2 <- (n > 0)*asset_value_2[region$get_indices()][,1]*0.4
  expected_recovery_delay_1 <- (n > 0)*1
  expected_recovery_delay_2 <- (n > 0)*2
  idx <- 1:region$get_locations() # idx <- 5901:6010
  expect_silent(n <- impacts_1$calculate(n, 0))
  expect_silent(n <- impacts_2$calculate(n, 0))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  expect_silent(n <- impacts_1$update_recovery_delay(n))
  expect_silent(n <- impacts_2$update_recovery_delay(n))
  expect_equal(attr(n, "recovery_delay")[[1]][idx],
               expected_recovery_delay_1[idx])
  expect_equal(attr(n, "recovery_delay")[[2]][idx],
               expected_recovery_delay_2[idx])
  idx1 <- c(5901:5910, 5951:5955)
  n[idx1] <- 0
  expected_recovery_delay_1[idx1] <- 0
  expected_recovery_delay_2[idx1] <- 1
  expect_silent(n <- impacts_1$calculate(n, 1))
  expect_silent(n <- impacts_2$calculate(n, 1))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  expect_silent(n <- impacts_1$update_recovery_delay(n))
  expect_silent(n <- impacts_2$update_recovery_delay(n))
  expect_equal(attr(n, "recovery_delay")[[1]][idx],
               expected_recovery_delay_1[idx])
  expect_equal(attr(n, "recovery_delay")[[2]][idx],
               expected_recovery_delay_2[idx])
  idx2 <- c(5911:5930, 5956:5960)
  n[idx2] <- 0
  expected_impact_1[idx1] <- 0
  expected_recovery_delay_1[idx2] <- 0
  expected_recovery_delay_2[c(idx1, idx2)] <-
    expected_recovery_delay_2[c(idx1, idx2)] - 1
  expect_silent(n <- impacts_1$calculate(n, 2))
  expect_silent(n <- impacts_2$calculate(n, 2))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  expect_silent(n <- impacts_1$update_recovery_delay(n))
  expect_silent(n <- impacts_2$update_recovery_delay(n))
  expect_equal(attr(n, "recovery_delay")[[1]][idx],
               expected_recovery_delay_1[idx])
  expect_equal(attr(n, "recovery_delay")[[2]][idx],
               expected_recovery_delay_2[idx])
  expected_impact_1[idx2] <- 0
  expected_impact_2[idx1] <- 0
  expected_recovery_delay_2[idx2] <- 0
  expect_silent(n <- impacts_1$calculate(n, 3))
  expect_silent(n <- impacts_2$calculate(n, 3))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  expect_silent(n <- impacts_1$update_recovery_delay(n))
  expect_silent(n <- impacts_2$update_recovery_delay(n))
  expect_equal(attr(n, "recovery_delay")[[1]][idx],
               expected_recovery_delay_1[idx])
  expect_equal(attr(n, "recovery_delay")[[2]][idx],
               expected_recovery_delay_2[idx])
  n[5921:5960] <- round(runif(40, 1, 10))
  expected_impact_1 <- (n > 0)*asset_value_1[region$get_indices()][,1]*0.3
  expected_impact_2 <- (n > 0)*asset_value_2[region$get_indices()][,1]*0.4
  expected_recovery_delay_1 <- (n > 0)*1
  expected_recovery_delay_2 <- (n > 0)*2
  expect_silent(n <- impacts_1$calculate(n, 4))
  expect_silent(n <- impacts_2$calculate(n, 4))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  expect_silent(n <- impacts_1$update_recovery_delay(n))
  expect_silent(n <- impacts_2$update_recovery_delay(n))
  expect_equal(attr(n, "recovery_delay")[[1]][idx],
               expected_recovery_delay_1[idx])
  expect_equal(attr(n, "recovery_delay")[[2]][idx],
               expected_recovery_delay_2[idx])
})

test_that("updates recovery delay to prolong density-based impacts", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template*0)
  template_vect <- template[region$get_indices()][,1]
  idx <- 5901:6000
  population_model <- UnstructPopulation(region, growth = 1.2,
                                         capacity = template_vect*50)
  asset_value_1 <- 100*(template > 0.1 & template < 0.4)
  asset_value_2 <- 200*(template > 0.3)
  expect_silent(impacts_1 <- Impacts(region, population_model,
                                     impact_type = "density",
                                     asset_name = "impact1",
                                     asset_value = asset_value_1,
                                     loss_rate = 0.3,
                                     recovery_delay = 2))
  expect_silent(impacts_1$set_id(1))
  expect_silent(impacts_2 <- Impacts(region, population_model,
                                     impact_type = "density",
                                     asset_name = "impact2",
                                     asset_value = asset_value_2,
                                     loss_rate = 0.4,
                                     recovery_delay = 3))
  expect_silent(impacts_2$set_id(2))
  n_density <- n <- rep(0, region$get_locations())
  n[idx] <- round(runif(100, 1, 10))
  n_orig <- n
  dens_idx <- idx[which(template_vect[idx] > 0)]
  n_density[dens_idx] <- pmin(n[dens_idx]/(template_vect[dens_idx]*50), 1)
  expected_impact_1 <- n_density*asset_value_1[region$get_indices()][,1]*0.3
  expected_impact_2 <- n_density*asset_value_2[region$get_indices()][,1]*0.4
  expect_silent(n <- impacts_1$calculate(n, 0))
  expect_silent(n <- impacts_2$calculate(n, 0))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  expect_silent(n <- impacts_1$update_recovery_delay(n))
  expect_silent(n <- impacts_2$update_recovery_delay(n))
  expect_equal(attr(n, "recovery_delay")[[1]], 2)
  expect_equal(attr(n, "recovery_delay")[[2]], 3)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"), list(n_density))
  expect_equal(attr(attr(n, "recovery_delay"), "max"), 3)
  expect_equal(attr(attr(n, "recovery_delay"), "first"), 1)
  i1 <- c(1:10, 51:55)
  n[idx[i1]] <- 0
  i6 <- c(11:20, 56:57)
  idx_1 <- idx[i6][which(n_density[idx[i6]] < 1)]
  n[idx_1] <- round(n[idx_1]*0.6)
  mask1 <- +(n > 0); mask1[idx[i6]] <- mask1[idx[i6]]*(n/n_orig)[idx[i6]]
  expect_silent(n <- impacts_1$calculate(n, 1))
  expect_silent(n <- impacts_2$calculate(n, 1))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  expect_silent(n <- impacts_1$update_recovery_delay(n))
  expect_silent(n <- impacts_2$update_recovery_delay(n))
  expect_equal(attr(n, "recovery_delay")[[1]], 2)
  expect_equal(attr(n, "recovery_delay")[[2]], 3)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"),
               list(mask1*n_density, n_density))
  expect_equal(attr(attr(n, "recovery_delay"), "max"), 3)
  expect_equal(attr(attr(n, "recovery_delay"), "first"), 1)
  i2 <- c(11:30, 56:60)
  n[idx[i2]] <- 0
  mask2 <- +(n > 0)
  expect_silent(n <- impacts_1$calculate(n, 2))
  expect_silent(n <- impacts_2$calculate(n, 2))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  expect_silent(n <- impacts_1$update_recovery_delay(n))
  expect_silent(n <- impacts_2$update_recovery_delay(n))
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"),
               list(mask2*n_density, mask1*n_density, n_density))
  mask3 <- +(n > 0)
  expected_impact_1[idx[i1]] <- 0
  expected_impact_1[idx[i6]] <- (expected_impact_1*mask1)[idx[i6]]
  expect_silent(n <- impacts_1$calculate(n, 3))
  expect_silent(n <- impacts_2$calculate(n, 3))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  expect_silent(n <- impacts_1$update_recovery_delay(n))
  expect_silent(n <- impacts_2$update_recovery_delay(n))
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"),
               list(mask3*n_density, mask2*n_density, mask1*n_density))
  i4 <- c(21:30, 58:60)
  n[idx[i4]] <- n_orig[idx[i4]]
  mask4 <- +(n > 0)
  expected_impact_1[idx[i6]] <- 0
  expected_impact_2[idx[i1]] <- 0
  expected_impact_2[idx[i6]] <- (expected_impact_2*mask1)[idx[i6]]
  expect_silent(n <- impacts_1$calculate(n, 4))
  expect_silent(n <- impacts_2$calculate(n, 4))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  expect_silent(n <- impacts_1$update_recovery_delay(n))
  expect_silent(n <- impacts_2$update_recovery_delay(n))
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"),
               list(mask4*n_density, mask3*n_density, mask2*n_density))
  mask5 <- +(n > 0)
  expected_impact_2[idx[i6]] <- 0
  expect_silent(n <- impacts_1$calculate(n, 5))
  expect_silent(n <- impacts_2$calculate(n, 5))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts[[1]][idx], expected_impact_1[idx])
  expect_equal(calc_impacts[[2]][idx], expected_impact_2[idx])
  expect_silent(n <- impacts_1$update_recovery_delay(n))
  expect_silent(n <- impacts_2$update_recovery_delay(n))
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"),
               list(mask5*n_density, mask4*n_density, mask3*n_density))
})

test_that("calculates spatially implicit impacts via area occupied", {
})

test_that("applies dynamic presence-based impacts", {
  TEST_DIRECTORY <- test_path("test_inputs")
})

test_that("applies dynamic density-based impacts with recovery", {
  TEST_DIRECTORY <- test_path("test_inputs")
})

test_that("applies dynamic implicit area-based impacts with recovery", {
})
