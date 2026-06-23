context("Impacts")

test_that("initializes with impacts, populations and impact stages", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- bsspread::Region(template)
  context <- bsimpact::Context("My species",
                               impact_scope = c("aspect1", "aspect2"))
  incursion <- bsimpact::Incursion(template*0, region)
  impact_layers <- list(aspect1 = 100*(template > 0.1 & template < 0.3),
                        aspect2 = 200*(template > 0.2 & template < 0.4))
  impacts <- bsimpact::ImpactAnalysis(context, region, incursion,
                                      impact_layers)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- bsspread::StagedPopulation(region, stage_matrix)
  expect_error(impacts_o <- Impacts("impacts", population_model),
               "Impacts must be a 'ImpactAnalysis' or inherited class object.")
  expect_error(impacts_o <- Impacts(impacts, "population_model"),
               paste("Population model must be a 'Population' or inherited",
                     "class object."))
  expect_error(impacts_o <- Impacts(impacts, population_model,
                                    impact_stages = 3:4),
               paste("Impact stages must be a vector of stage indices",
                     "consistent with the population model."))
  expect_error(impacts_o <- Impacts(impacts, population_model,
                                    calc_total = 1),
               "Calculate total indicator should be logical.")
  expect_error(Impacts(impacts, population_model,
                       dynamic_links = "suitability"),
               paste("Dynamic links are only applicable for dynamic impacts",
                     "(i.e. a bsimpact::ValueImpacts object with is_dynamic",
                     "= TRUE)."), fixed = TRUE)
  impacts <- bsimpact::ValueImpacts(context, region, incursion,
                                    impact_layers,
                                    loss_rates = c(aspect1 = 0.3,
                                                   aspect2 = 0.4),
                                    is_dynamic = TRUE)
  expect_error(Impacts(impacts, population_model,
                       dynamic_links = "dummy"),
               paste("Dynamic links should be a vector of one or more of",
                     "'suitability', 'capacity', and/or 'attractors'."))
  expect_error(impacts_o <- Impacts(impacts, population_model,
                                    recovery_delay = "1"),
               "Recover delay should a number >= 0.")
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     impact_stages = 2:3,
                                     dynamic_links = "suitability",
                                     recovery_delay = 2))
  expect_is(impacts_o, "Impacts")
  expect_named(impacts_o,
               c("get_impacts", "get_context", "get_calc_total",
                 "includes_combined", "update_recovery_delay", "calculate"))
  expect_is(impacts_o$get_impacts(), "ValueImpacts")
  expect_is(impacts_o$get_context(), "Context")
  expect_true(impacts_o$includes_combined())
  expect_true(impacts_o$get_calc_total())
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     calc_total = FALSE))
  expect_false(impacts_o$get_calc_total())
  context <- bsimpact::Context("My species",
                               impact_scope = c("aspect1", "aspect2"),
                               valuation_type = "non-monetary")
  impacts <- bsimpact::ImpactAnalysis(context, region, incursion,
                                      impact_layers)
  expect_silent(impacts_o <- Impacts(impacts, population_model))
  expect_false(impacts_o$get_calc_total())
  impacts <- bsimpact::ImpactAnalysis(context, region, incursion,
                                      impact_layers, combine_function = "none")
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     calc_total = TRUE))
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     calc_total = FALSE))
  expect_false(impacts_o$get_calc_total())
  context <- bsimpact::Context("My species", impact_scope = "aspect1",
                               valuation_type = "non-monetary")
  impacts <- bsimpact::ImpactAnalysis(context, region, incursion,
                                      impact_layers[1],
                                      combine_function = "none")
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     calc_total = TRUE))
  expect_true(impacts_o$get_calc_total())
})

test_that("calculates impacts including combined", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- bsspread::Region(template)
  template[region$get_indices()][5922,] <- 0.25
  context <- bsimpact::Context("My species",
                               impact_scope = c("aspect1", "aspect2"))
  incursion <- bsimpact::Incursion(template*0, region, type = "presence",
                                   multiplier = 0.2)
  aspects <- list(aspect1 = "aspect1", aspect2 = "aspect2")
  impact_layers <- list(aspect1 = 100*(template > 0.1 & template < 0.3),
                        aspect2 = 200*(template > 0.2 & template < 0.4))
  loss_rates <- c(aspect1 = 0.3, aspect2 = 0.4)
  impacts <- bsimpact::ValueImpacts(context, region, incursion,
                                    impact_layers, loss_rates = loss_rates,
                                    discount_rate = 0.05)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  population_model <- bsspread::StagedPopulation(region, stage_matrix)
  initial_n <- rep(0, region$get_locations())
  initial_n[5920:5922] <- 10:12
  initializer <- bsspread::Initializer(initial_n, region = region,
                                       population_model = population_model)
  n <- initializer$initialize()
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     impact_stages = 2:3))
  expected_impacts <- lapply(aspects, function(l) {
    impact_incursion <- rowSums(n[,2:3])*0.2
    impact_incursion[which(impact_incursion > 1)] <- 1
    (impact_incursion*impact_layers[[l]][region$get_indices()][,1]*
        loss_rates[l]/(1.05^4))
  })
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expect_silent(n <- impacts_o$calculate(n, 4))
  expect_true("impacts" %in% names(attributes(n)))
  calc_impacts <- attr(n, "impacts")
  expect_named(calc_impacts, c("aspect1", "aspect2", "combined"))
  expect_equal(calc_impacts, expected_impacts)
  # density-based impacts
  incursion <- bsimpact::Incursion(template*0, region, type = "density")
  impacts <- bsimpact::ValueImpacts(context, region, incursion,
                                    impact_layers, loss_rates = loss_rates)
  expect_error(impacts_o <- Impacts(impacts, population_model,
                                    impact_stages = 2:3),
               "Population capacity is required for density-based impacts.")
  population_model <- bsspread::StagedPopulation(region, stage_matrix,
                                                 capacity = (initial_n > 0)*5)
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     impact_stages = 2:3))
  n[5920:5922,] <- c(4, 9, 9, 4, 2, 3, 2, 0, 0)
  expected_impacts <- lapply(aspects, function(l) {
    impact_incursion <- pmin(rowSums(n[,2:3])/5, 1)
    (impact_incursion*impact_layers[[l]][region$get_indices()][,1]*
        loss_rates[l])
  })
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expect_silent(n <- impacts_o$calculate(n, 4))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts, expected_impacts)
})

test_that("updates recovery delay to prolong presence-based impacts", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- bsspread::Region(template*0)
  context <- bsimpact::Context("My species",
                               impact_scope = c("aspect1", "aspect2"))
  incursion <- bsimpact::Incursion(template*0, region)
  aspects <- list(aspect1 = "aspect1", aspect2 = "aspect2")
  impact_layers <- list(aspect1 = 100*(template > 0.1 & template < 0.4),
                        aspect2 = 200*(template > 0.3))
  loss_rates = c(aspect1 = 0.3, aspect2 = 0.4)
  impacts <- bsimpact::ValueImpacts(context, region, incursion,
                                    impact_layers, loss_rates = loss_rates)
  impacts$set_id(3)
  population_model <- bsspread::UnstructPopulation(region, growth = 1.2)
  n <- rep(0, region$get_locations())
  n[5901:6000] <- round(runif(100, 1, 10))
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     recovery_delay = 2))
  expected_impacts <- lapply(aspects, function(l) {
    impact_incursion <- +(n > 0)
    (impact_incursion*impact_layers[[l]][region$get_indices()][,1]*
        loss_rates[l])
  })
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expected_recovery_delay <- (n > 0)*2
  idx <- 1:region$get_locations()
  expect_silent(n <- impacts_o$calculate(n, 0))
  calc_impacts <- attr(n, "impacts")
  expect_equal(lapply(calc_impacts, function(impact) impact[idx]),
               lapply(expected_impacts, function(impact) impact[idx]))
  expect_equal(attr(n, "recovery_delay")[[3]][idx],
               expected_recovery_delay[idx])
  n[5901:5910] <- 0
  expected_recovery_delay[5901:5910] <- 1
  expect_silent(n <- impacts_o$calculate(n, 1))
  calc_impacts <- attr(n, "impacts")
  expect_equal(lapply(calc_impacts, function(impact) impact[idx]),
               lapply(expected_impacts, function(impact) impact[idx]))
  expect_equal(attr(n, "recovery_delay")[[3]][idx],
               expected_recovery_delay[idx])
  n[5911:5930] <- 0
  expected_recovery_delay[5901:5930] <- expected_recovery_delay[5901:5930] - 1
  expect_silent(n <- impacts_o$calculate(n, 2))
  calc_impacts <- attr(n, "impacts")
  expect_equal(lapply(calc_impacts, function(impact) impact[idx]),
               lapply(expected_impacts, function(impact) impact[idx]))
  expect_equal(attr(n, "recovery_delay")[[3]][idx],
               expected_recovery_delay[idx])
  expected_impacts <- lapply(expected_impacts, function(impact) {
    impact[5901:5910] <- 0
    impact
  })
  expected_recovery_delay[5911:5930] <- expected_recovery_delay[5911:5930] - 1
  expect_silent(n <- impacts_o$calculate(n, 3))
  calc_impacts <- attr(n, "impacts")
  expect_equal(lapply(calc_impacts, function(impact) impact[idx]),
               lapply(expected_impacts, function(impact) impact[idx]))
  expect_equal(attr(n, "recovery_delay")[[3]][idx],
               expected_recovery_delay[idx])
  n[5921:5930] <- round(runif(10, 1, 10))
  expected_impacts <- lapply(expected_impacts, function(impact) {
    impact[5911:5920] <- 0
    impact
  })
  expected_recovery_delay[5921:5930] <- 2
  expect_silent(n <- impacts_o$calculate(n, 4))
  calc_impacts <- attr(n, "impacts")
  expect_equal(lapply(calc_impacts, function(impact) impact[idx]),
               lapply(expected_impacts, function(impact) impact[idx]))
  expect_equal(attr(n, "recovery_delay")[[3]][idx],
               expected_recovery_delay[idx])
})

test_that("updates recovery delay to prolong density-based impacts", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- bsspread::Region(template*0)
  template_vect <- template[region$get_indices()][,1]
  idx <- 5901:6000
  context <- bsimpact::Context("My species",
                               impact_scope = c("aspect1", "aspect2"))
  incursion <- bsimpact::Incursion(template*0, region, type = "density")
  aspects <- list(aspect1 = "aspect1", aspect2 = "aspect2")
  impact_layers <- list(aspect1 = 100*(template > 0.1 & template < 0.4),
                        aspect2 = 200*(template > 0.3))
  loss_rates = c(aspect1 = 0.3, aspect2 = 0.4)
  impacts <- bsimpact::ValueImpacts(context, region, incursion,
                                    impact_layers, loss_rates = loss_rates)
  impacts$set_id(2)
  impacts_2 <- bsimpact::ValueImpacts(context, region, incursion,
                                      impact_layers, loss_rates = loss_rates)
  impacts_2$set_id(1)
  population_model <- bsspread::UnstructPopulation(region, growth = 1.2,
                                                   capacity = template_vect*50)
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     recovery_delay = 2))
  expect_silent(impacts_o_2 <- Impacts(impacts_2, population_model,
                                       recovery_delay = 3))
  n_density <- n <- rep(0, region$get_locations())
  n[idx] <- round(runif(100, 1, 10))
  n_orig <- n
  dens_idx <- idx[which(template_vect[idx] > 0)]
  n_density[dens_idx] <- pmin(n[dens_idx]/(template_vect[dens_idx]*50), 1)
  expected_impacts <- lapply(aspects, function(l) {
    n_density*impact_layers[[l]][region$get_indices()][,1]*loss_rates[l]
  })
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expect_silent(n <- impacts_o$calculate(n, 0))
  calc_impacts <- attr(n, "impacts")
  expect_silent(n <- impacts_o_2$update_recovery_delay(n))
  expect_equal(calc_impacts, expected_impacts)
  expect_equal(attr(n, "recovery_delay")[[1]], 3)
  expect_equal(attr(n, "recovery_delay")[[2]], 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"), list(n_density))
  expect_equal(attr(attr(n, "recovery_delay"), "max"), 3)
  expect_equal(attr(attr(n, "recovery_delay"), "first"), 2)
  n[idx[1:10]] <- 0
  idx_1 <- idx[11:20][which(n_density[idx[11:20]] < 1)]
  n[idx_1] <- round(n[idx_1]*0.6)
  mask1 <- +(n > 0)
  mask1[idx[11:20]] <- mask1[idx[11:20]]*(n/n_orig)[idx[11:20]]
  expect_silent(n <- impacts_o$calculate(n, 1))
  calc_impacts <- attr(n, "impacts")
  expect_silent(n <- impacts_o_2$update_recovery_delay(n))
  expect_equal(calc_impacts, expected_impacts)
  expect_equal(attr(n, "recovery_delay")[[1]], 3)
  expect_equal(attr(n, "recovery_delay")[[2]], 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"),
               list(mask1*n_density, n_density))
  n[idx[11:30]] <- 0
  mask2 <- +(n > 0)
  expect_silent(n <- impacts_o$calculate(n, 2))
  calc_impacts <- attr(n, "impacts")
  expect_silent(n <- impacts_o_2$update_recovery_delay(n))
  expect_equal(calc_impacts, expected_impacts)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"),
               list(mask2*n_density, mask1*n_density, n_density))
  mask3 <- +(n > 0)
  expected_impacts <- lapply(expected_impacts, function(impact) {
    impact[idx[1:10]] <- 0
    impact[idx[11:20]] <- (impact*mask1)[idx[11:20]]
    impact
  })
  expect_silent(n <- impacts_o$calculate(n, 3))
  calc_impacts <- attr(n, "impacts")
  expect_silent(n <- impacts_o_2$update_recovery_delay(n))
  expect_equal(calc_impacts, expected_impacts)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"),
               list(mask3*n_density, mask2*n_density, mask1*n_density))
  n[idx[21:30]] <- n_orig[idx[21:30]]
  mask4 <- +(n > 0)
  expected_impacts <- lapply(expected_impacts, function(impact) {
    impact[idx[11:20]] <- 0
    impact
  })
  expect_silent(n <- impacts_o$calculate(n, 0))
  calc_impacts <- attr(n, "impacts")
  expect_silent(n <- impacts_o_2$update_recovery_delay(n))
  expect_equal(calc_impacts, expected_impacts)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"),
               list(mask4*n_density, mask3*n_density, mask2*n_density))
})

test_that("calculates spatially implicit impacts via area occupied", {
  context <- bsimpact::Context("My species",
                               impact_scope = c("aspect1", "aspect2"))
  incursion <- bsimpact::Incursion(0, bsimpact::Region(), type = "area")
  aspects <- list(aspect1 = "aspect1", aspect2 = "aspect2")
  impact_layers <- list(aspect1 = 100, aspect2 = 200)
  loss_rates <- c(aspect1 = 0.3, aspect2 = 0.4)
  impacts <- bsimpact::ValueImpacts(context, bsimpact::Region(), incursion,
                                    impact_layers, loss_rates = loss_rates)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  region <- bsspread::Region()
  population_model <- bsspread::StagedPopulation(region, stage_matrix)
  initializer <- bsspread::Initializer(20, region = region,
                                       population_model = population_model)
  n <- initializer$initialize()
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     impact_stages = 2:3))
  expect_error(n <- impacts_o$calculate(n, 4),
               paste("Cannot calculate spatially implicit impacts without",
                     "area occupied."))
  attr(n, "spread_area") <- 50
  expect_silent(n <- impacts_o$calculate(n, 4))
  calc_impacts <- attr(n, "impacts")
  expect_named(calc_impacts, c("aspect1", "aspect2", "combined"))
  expect_equal(calc_impacts, list(aspect1 = 100*50*0.3, aspect2 = 200*50*0.4,
                                  combined = 100*50*0.3 + 200*50*0.4))
  # with recovery delay
  impacts$set_id(3)
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     impact_stages = 2:3,
                                     recovery_delay = 2))
  n <- initializer$initialize()
  attr(n, "spread_area") <- 50
  expect_silent(n <- impacts_o$calculate(n, 0))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts, list(aspect1 = 100*50*0.3, aspect2 = 200*50*0.4,
                                  combined = 100*50*0.3 + 200*50*0.4))
  expect_equal(attr(n, "recovery_delay")[[3]], 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"), 50)
  attr(n, "spread_area") <- 30
  expect_silent(n <- impacts_o$calculate(n, 1))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts, list(aspect1 = 100*50*0.3, aspect2 = 200*50*0.4,
                                  combined = 100*50*0.3 + 200*50*0.4))
  expect_equal(attr(n, "recovery_delay")[[3]], 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"), c(30, 50))
  n[,2:3] <- 0
  attr(n, "spread_area") <- 20
  expect_silent(n <- impacts_o$calculate(n, 2))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts, list(aspect1 = 100*50*0.3, aspect2 = 200*50*0.4,
                                  combined = 100*50*0.3 + 200*50*0.4))
  expect_equal(attr(n, "recovery_delay")[[3]], 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"), c(0, 30, 50))
  expect_silent(n <- impacts_o$calculate(n, 3))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts, list(aspect1 = 100*30*0.3, aspect2 = 200*30*0.4,
                                  combined = 100*30*0.3 + 200*30*0.4))
  expect_equal(attr(n, "recovery_delay")[[3]], 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"), c(0, 0, 30, 50))
  expect_silent(n <- impacts_o$calculate(n, 4))
  calc_impacts <- attr(n, "impacts")
  expect_equal(calc_impacts, list(aspect1 = 0, aspect2 = 0, combined = 0))
})

test_that("applies dynamic presence-based impacts", {
  TEST_DIRECTORY <- test_path("test_inputs")
  context <- bsimpact::Context("My species",
                               impact_scope = c("aspect1", "aspect2"))
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- bsspread::Region(template*0)
  template_vect <- template[region$get_indices()][,1]
  idx <- 5901:6000
  incursion <- bsimpact::Incursion(template > 0, region)
  aspects <- list(aspect1 = "aspect1", aspect2 = "aspect2")
  impact_layers <- list(aspect1 = 100*(template > 0.1 & template < 0.4),
                        aspect2 = 200*(template > 0.3))
  loss_rates = c(aspect1 = 0.3, aspect2 = 0.4)
  impacts <- bsimpact::ValueImpacts(context, region, incursion, impact_layers,
                                    loss_rates = loss_rates, is_dynamic = TRUE)
  impacts$set_id(2)
  population_model <- bsspread::UnstructPopulation(region, growth = 1.2)
  expect_silent(impacts_o <-
                  Impacts(impacts, population_model,
                          dynamic_links = c("suitability", "capacity")))
  n <- rep(0, region$get_locations())
  n[idx[1:90]] <- round(runif(90, 1, 10))
  n[idx[1:10]] <- 0
  expected_impacts <- lapply(aspects, function(l) {
    impact_incursion <- +(n > 0)
    (impact_incursion*impact_layers[[l]][region$get_indices()][,1]*
        loss_rates[l])
  })
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expected_dynamic_mult <- lapply(aspects,
                                  function(a) (1 - loss_rates[a]*(n > 0)))
  attr(expected_dynamic_mult, "links") <- c("suitability", "capacity")
  expect_silent(n <- impacts_o$calculate(n, 0))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  n[idx[1:10]] <- 1
  n[idx[11:20]] <- 0
  expected_dynamic_mult <- lapply(aspects, function(a) {
    dynamic_mult <- expected_dynamic_mult[[a]]*(1 - loss_rates[a]*(n > 0))
    dynamic_mult[which(n == 0)] <- 1
    dynamic_mult
  })
  attr(expected_dynamic_mult, "links") <- c("suitability", "capacity")
  expected_impacts <-
    lapply(aspects, function(a) {
      (impact_layers[[a]][region$get_indices()][,1]*
         (1 - expected_dynamic_mult[[a]]))})
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expect_silent(n <- impacts_o$calculate(n, 1))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
})

test_that("applies dynamic density-based impacts with recovery", {
  TEST_DIRECTORY <- test_path("test_inputs")
  context <- bsimpact::Context("My species",
                               impact_scope = c("aspect1", "aspect2"))
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- bsspread::Region(template*0)
  template_vect <- template[region$get_indices()][,1]
  idx <- 5901:6000
  incursion <- bsimpact::Incursion(template*0, region, type = "density")
  aspects <- list(aspect1 = "aspect1", aspect2 = "aspect2")
  impact_layers <- list(aspect1 = 100*(template > 0.1 & template < 0.4),
                        aspect2 = 200*(template > 0.3))
  loss_rates = c(aspect1 = 0.3, aspect2 = 0.4)
  impacts <- bsimpact::ValueImpacts(context, region, incursion, impact_layers,
                                    loss_rates = loss_rates, is_dynamic = TRUE)
  impacts$set_id(2)
  population_model <- bsspread::UnstructPopulation(region, growth = 1.2,
                                                   capacity = template_vect*50)
  n_density <- n <- rep(0, region$get_locations())
  n[idx[1:90]] <- round(runif(90, 1, 10))
  dens_idx <- idx[which(population_model$get_capacity()[idx] > 0)]
  n_density[dens_idx] <- pmin(n[dens_idx]/
                                population_model$get_capacity()[dens_idx], 1)
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     recovery_delay = 2))
  expected_dynamic_mult <-
    lapply(aspects, function(a) (1 - loss_rates[a]*n_density))
  expected_impacts <-
    lapply(aspects, function(a) {
      (impact_layers[[a]][region$get_indices()][,1]*loss_rates[a]*n_density)})
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expected_recovery_delay <- (n_density > 0)*2
  expect_silent(n <- impacts_o$calculate(n, 0))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  expect_equal(attr(n, "recovery_delay")[[2]], expected_recovery_delay)
  attr(n, "impacts") <- NULL
  n[idx[1:10]] <- 0
  idx_1 <- idx[11:20][which(n_density[idx[11:20]] < 1)]
  n[idx_1] <- round(n[idx_1]*0.6)
  n_density[dens_idx] <- pmin(n[dens_idx]/
                                population_model$get_capacity()[dens_idx], 1)
  expected_dynamic_mult <- lapply(aspects, function(a) {
    expected_dynamic_mult[[a]]*(1 - loss_rates[a]*n_density)
  })
  expected_impacts <-
    lapply(aspects, function(a) {
      (impact_layers[[a]][region$get_indices()][,1]*
         (1 - expected_dynamic_mult[[a]]))})
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expected_recovery_delay[idx[1:10]] <- 1
  expect_silent(n <- impacts_o$calculate(n, 1))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  expect_equal(attr(n, "recovery_delay")[[2]], expected_recovery_delay)
  attr(n, "impacts") <- NULL
  n[idx[11:30]] <- 0
  n_density[dens_idx] <- pmin(n[dens_idx]/
                                population_model$get_capacity()[dens_idx], 1)
  expected_dynamic_mult <- lapply(aspects, function(a) {
    expected_dynamic_mult[[a]]*(1 - loss_rates[a]*n_density)
  })
  expected_impacts <-
    lapply(aspects, function(a) {
      (impact_layers[[a]][region$get_indices()][,1]*
         (1 - expected_dynamic_mult[[a]]))})
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expected_recovery_delay[idx[1:30]] <-
    pmax(expected_recovery_delay[idx[1:30]] - 1, 0)
  expect_silent(n <- impacts_o$calculate(n, 2))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  expect_equal(attr(n, "recovery_delay")[[2]], expected_recovery_delay)
  attr(n, "impacts") <- NULL
  expected_dynamic_mult <- lapply(aspects, function(a) {
    mult <- expected_dynamic_mult[[a]]*(1 - loss_rates[a]*n_density)
    mult[idx[1:10]] <- 1
    mult
  })
  expected_impacts <- lapply(aspects, function(a) {
    impact <- (impact_layers[[a]][region$get_indices()][,1]*
                 (1 - expected_dynamic_mult[[a]]))
    impact[idx[1:10]] <- 0
    impact
  })
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expected_recovery_delay[idx[11:30]] <-
    pmax(expected_recovery_delay[idx[11:30]] - 1, 0)
  expect_silent(n <- impacts_o$calculate(n, 3))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  expect_equal(attr(n, "recovery_delay")[[2]], expected_recovery_delay)
  attr(n, "impacts") <- NULL
  n[idx[21:30]] <- round(runif(10, 1, 10))
  n_density[dens_idx] <- pmin(n[dens_idx]/
                                population_model$get_capacity()[dens_idx], 1)
  expected_dynamic_mult <- lapply(aspects, function(a) {
    mult <- expected_dynamic_mult[[a]]*(1 - loss_rates[a]*n_density)
    mult[idx[11:20]] <- 1
    mult
  })
  expected_impacts <- lapply(aspects, function(a) {
    impact <- (impact_layers[[a]][region$get_indices()][,1]*
                 (1 - expected_dynamic_mult[[a]]))
    impact[idx[11:20]] <- 0
    impact
  })
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expected_recovery_delay[idx[21:30]] <- (n_density[idx[21:30]] > 0)*2
  expect_silent(n <- impacts_o$calculate(n, 4))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  expect_equal(attr(n, "recovery_delay")[[2]], expected_recovery_delay)
})

test_that("applies dynamic implicit area-based impacts with recovery", {
  context <- bsimpact::Context("My species",
                               impact_scope = c("aspect1", "aspect2"))
  incursion <- bsimpact::Incursion(0, bsimpact::Region(), type = "area")
  aspects <- list(aspect1 = "aspect1", aspect2 = "aspect2")
  impact_layers <- list(aspect1 = 100, aspect2 = 200)
  loss_rates <- c(aspect1 = 0.3, aspect2 = 0.4)
  impacts <- bsimpact::ValueImpacts(context, bsimpact::Region(), incursion,
                                    impact_layers, loss_rates = loss_rates,
                                    is_dynamic = TRUE)
  impacts$set_id(2)
  stage_matrix <- matrix(c(0.0, 2.0, 5.0,
                           0.3, 0.0, 0.0,
                           0.0, 0.6, 0.8),
                         nrow = 3, ncol = 3, byrow = TRUE)
  region <- bsspread::Region()
  population_model <- bsspread::StagedPopulation(region, stage_matrix,
                                                 capacity = 50,
                                                 capacity_area = 1,
                                                 capacity_stages = 2:3)
  expect_silent(impacts_o <-
                  Impacts(impacts, population_model,
                          impact_stages = 2:3,
                          dynamic_links = c("suitability", "capacity"),
                          recovery_delay = 2))
  initializer <- bsspread::Initializer(100, region = region,
                                       population_model = population_model)
  n <- initializer$initialize()
  attr(n, "spread_area") <- 50
  expected_dynamic_mult <- list(aspect1 = 1 - 0.3, aspect2 = 1 - 0.4)
  attr(expected_dynamic_mult, "incursion") <- 50
  attr(expected_dynamic_mult, "links") <- c("suitability", "capacity")
  expected_impacts <- list(aspect1 = 100*50*0.3, aspect2 = 200*50*0.4)
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expect_silent(n <- impacts_o$calculate(n, 0))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  expect_equal(as.numeric(attr(n, "recovery_delay")[[2]]), 2)
  expect_equal(attr(attr(n, "recovery_delay"), "max"), 2)
  expect_equal(attr(attr(n, "recovery_delay"), "first"), 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"), 50)
  expect_equal(attr(attr(n, "recovery_delay")[[2]], "dynamic_mult"),
               list(aspect1 = 0.7, aspect2 = 0.6))
  expect_silent(impacts_o <- Impacts(impacts, population_model,
                                     impact_stages = 2:3,
                                     recovery_delay = 2))
  attr(n, "impacts") <- NULL
  attr(n, "spread_area") <- 30
  expected_dynamic_mult <- lapply(expected_dynamic_mult, function(a) a^2)
  attr(expected_dynamic_mult, "incursion") <- 30
  expected_impacts <- list(aspect1 = 100*max(30*0.51, 50*0.3),
                           aspect2 = 200*max(30*0.64, 50*0.4))
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expect_silent(n <- impacts_o$calculate(n, 1))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  expect_equal(as.numeric(attr(n, "recovery_delay")[[2]]), 2)
  expect_equal(attr(attr(n, "recovery_delay"), "max"), 2)
  expect_equal(attr(attr(n, "recovery_delay"), "first"), 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"), c(30, 50))
  expect_equal(attr(attr(n, "recovery_delay")[[2]], "dynamic_mult"),
               list(aspect1 = c(0.49, 0.7), aspect2 = c(0.36, 0.6)))
  attr(n, "impacts") <- NULL
  attr(n, "spread_area") <- 50
  expected_dynamic_mult <- list(aspect1 = (50 - 30*(1 - 0.49))*0.7/50,
                                aspect2 = (50 - 30*(1 - 0.36))*0.6/50)
  attr(expected_dynamic_mult, "incursion") <- 50
  expected_impacts <- list(aspect1 = 100*max(50*(1 - 0.4858), 30*0.51, 50*0.3),
                           aspect2 = 200*max(50*(1 - 0.3696), 30*0.64, 50*0.4))
  expected_impacts$combined <-
    expected_impacts$aspect1 + expected_impacts$aspect2
  expect_silent(n <- impacts_o$calculate(n, 2))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  expect_equal(as.numeric(attr(n, "recovery_delay")[[2]]), 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"), c(50, 30, 50))
  expect_equal(attr(attr(n, "recovery_delay")[[2]], "dynamic_mult"),
               list(aspect1 = c(0.4858, 0.49, 0.7),
                    aspect2 = c(0.3696, 0.36, 0.6)))
  attr(n, "impacts") <- NULL
  attr(n, "spread_area") <- 0
  attr(expected_dynamic_mult, "incursion") <- 0
  expect_silent(n <- impacts_o$calculate(n, 3))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  expect_equal(as.numeric(attr(n, "recovery_delay")[[2]]), 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"), c(0, 50, 30, 50))
  expect_equal(attr(attr(n, "recovery_delay")[[2]], "dynamic_mult"),
               list(aspect1 = c(0.4858, 0.4858, 0.49, 0.7),
                    aspect2 = c(0.3696, 0.3696, 0.36, 0.6)))
  attr(n, "impacts") <- NULL
  expect_silent(n <- impacts_o$calculate(n, 4))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  expect_equal(as.numeric(attr(n, "recovery_delay")[[2]]), 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"),
               c(0, 0, 50, 30, 50))
  expect_equal(attr(attr(n, "recovery_delay")[[2]], "dynamic_mult"),
               list(aspect1 = c(0.4858, 0.4858, 0.4858, 0.49, 0.7),
                    aspect2 = c(0.3696, 0.3696, 0.3696, 0.36, 0.6)))
  attr(n, "impacts") <- NULL
  expected_dynamic_mult <- list(aspect1 = 1, aspect2 = 1)
  attr(expected_dynamic_mult, "incursion") <- 0
  expected_impacts <- list(aspect1 = 0, aspect2 = 0, combined = 0)
  expect_silent(n <- impacts_o$calculate(n, 5))
  expect_equal(attr(n, "impacts"), expected_impacts)
  expect_equal(attr(n, "dynamic_mult")[[2]], expected_dynamic_mult)
  expect_equal(as.numeric(attr(n, "recovery_delay")[[2]]), 2)
  expect_equal(attr(attr(n, "recovery_delay"), "incursions"),
               c(0, 0, 0, 50, 30, 50))
  expect_equal(attr(attr(n, "recovery_delay")[[2]], "dynamic_mult"),
               list(aspect1 = c(1, 0.4858, 0.4858, 0.4858, 0.49, 0.7),
                    aspect2 = c(1, 0.3696, 0.3696, 0.3696, 0.36, 0.6)))
})
