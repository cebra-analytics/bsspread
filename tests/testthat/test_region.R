context("Region")

test_that("initializes with planar CRS raster", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  expect_silent(region <- Region(template))
  expect_equal(region$get_type(), "grid")
  expect_equal(region$get_locations(), length(which(is.finite(template[]))))
  expect_false(region$spatially_implicit())
  expect_true(region$is_compatible(region$get_template()))
  expect_equal(region$get_indices(), which(is.finite(template[])))
  expect_equal(region$get_res(), 1000)
  expect_equal(region$is_included(5922:5925), c(TRUE, FALSE, FALSE, TRUE))
})

test_that("initializes with lonlat raster", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb_wgs84.tif"))
  expect_silent(region <- Region(template))
  coord <- terra::xyFromCell(template, 1113)
  mid_res <- terra::distance(coord, coord + 0.025, lonlat = TRUE)/sqrt(2)
  expect_true(abs(region$get_res() - mid_res)/mid_res < 0.05)
  expect_equal(region$is_included(1113:1116), c(TRUE, FALSE, FALSE, TRUE))
})

test_that("gets a raster with specified values", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  expect_silent(region <- Region(template))
  expect_error(region$get_rast(1:10),
               paste("Values should be a single value or vector of length",
                     "matching the number of region locations."))
  expect_silent(value_rast <- region$get_rast(10))
  expect_equal(value_rast[region$get_indices()][,1],
               rep(10, region$get_locations()))
  expect_silent(value_rast <- region$get_rast(1:region$get_locations()))
  expect_equal(value_rast[region$get_indices()][,1], 1:region$get_locations())
  expect_true(terra::is.int(value_rast))
  region <- Region(terra::as.int(template*10))
  expect_true(terra::is.int(region$get_template()))
  expect_silent(value_rast <- region$get_rast(
    round(runif(region$get_locations()), 1)))
  expect_false(terra::is.int(value_rast))
  expect_true(all(
    unique(value_rast[region$get_indices()][,1]) %in% ((0:10)/10)))
  cat_values <- factor(c("a", rep(c("b", "c"), 5915)), c("a", "b", "c", "d"))
  expect_silent(value_rast <- region$get_rast(cat_values))
  expect_true(terra::is.factor(value_rast))
  expect_equal(terra::cats(value_rast)[[1]],
               data.frame(ID = 1:4, category = c("a", "b", "c", "d")))
  expect_equal(value_rast[region$get_indices()][,1], cat_values)
  locs <- rep(NA, region$get_locations())
  locs[5900:5922] <- 1
  expected_rast <- terra::buffer(region$get_rast(locs), width = 5000)
  expect_silent(idx <- region$get_nearby(5900:5922, 5000))
  expect_equal(idx, which(expected_rast[region$get_indices()][,1] == 1))
})

test_that("initializes with CSV data", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  expect_silent(region <- Region(locations))
  expect_equal(region$get_type(), "patch")
  expect_equal(region$get_locations(), nrow(locations))
  expect_false(region$spatially_implicit())
  expect_true(region$is_compatible(1:nrow(locations)))
  expect_false(region$two_tier())
  expect_equal(region$get_coords(), locations[, c("lon", "lat")])
  expect_equal(region$get_coords(extra_cols = TRUE),
               locations[, c("lon", "lat", "name")])
  expect_silent(idx <- region$get_nearby(1:3, 180000))
  expect_equal(idx, which(rowSums(
    terra::distance(locations[,2:1], locations[1:3,2:1],
                    lonlat = TRUE) <= 180000) > 0))
})

test_that("empty initialize for single aspatial patch", {
  expect_silent(region <- Region())
  expect_equal(region$get_type(), "patch")
  expect_equal(region$get_locations(), 1)
  expect_true(region$spatially_implicit())
  expect_null(region$get_max_implicit_area())
  expect_error(region$set_max_implicit_area(0),
               "The maximum spatially implicit area must be numeric and > 0.")
  expect_silent(region$set_max_implicit_area(1e8))
  expect_equal(region$get_max_implicit_area(), 1e8)
  expect_silent(region$configure_paths())
  expect_silent(region$calculate_paths(1))
  expect_silent(paths <- region$get_paths(1))
  expect_length(paths$idx[["1"]], 0)
  expect_length(paths$distances[["1"]], 0)
})

test_that("creates two tier aggregation", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  expect_is(region, "Region")
  expect_false(region$two_tier())
  expect_silent(region$set_aggr(aggr_factor = 5, inner_radius = 10000))
  expect_true(region$two_tier())
  aggr <- region$get_aggr()
  expect_true(all(unlist(lapply(aggr$cells, length)) <= 25))
  expect_is(aggr$rast, "SpatRaster")
  expect_is(aggr$pts, "SpatVector")
  expect_is(aggr$get_cells, "function")
})

test_that("creates paths to single tier cells", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  expect_silent(region$configure_paths(directions = TRUE))
  expect_silent(region$calculate_paths(5922))
  expect_silent(paths <- region$get_paths(5922, directions = TRUE))
  expect_true(all(paths$idx$`5922`$cell %in% 1:region$get_locations()))
  expect_equal(length(paths$idx$`5922`$cell), region$get_locations() - 1)
  expect_true(all(paths$directions$`5922`$cell <= 360))
  region <- Region(template)
  expect_silent(region$configure_paths(max_distance = 20000))
  expect_silent(region$calculate_paths(5922))
  expect_silent(paths <- region$get_paths(5922, max_distance = 20000))
  expect_equal(length(paths$idx$`5922`$cell), 937)
  expect_true(all(paths$distances$`5922`$cell <= 20000))
})

test_that("creates paths to two tier cells and aggregate cells", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  region$set_aggr(aggr_factor = 5, inner_radius = 10000)
  aggr <- region$get_aggr()
  expect_silent(region$configure_paths(directions = TRUE))
  expect_silent(region$calculate_paths(5922))
  paths <- region$get_paths(5922, directions = TRUE, max_distance = 30000)
  expect_true(all(paths$distances$`5922`$cell <= 10000 + 5000 + 1000))
  expect_true(all(paths$idx$`5922`$aggr %in% 1:length(aggr$indices)))
  expect_true(all(paths$distances$`5922`$aggr <= 30000))
  expect_true(all(paths$directions$`5922`$aggr <= 360))
})

test_that("modifies paths via permeability layer", {
  TEST_DIRECTORY <- test_path("test_inputs")
  template <- terra::rast(file.path(TEST_DIRECTORY, "greater_melb.tif"))
  region <- Region(template)
  region$set_aggr(aggr_factor = 5, inner_radius = 10000)
  aggr <- region$get_aggr()
  perm_rast <- region$get_template()
  perm_rast[region$get_indices()] <- round(region$get_indices()/1560)/10
  perm <- Permeability(perm_rast, region)
  expect_null(perm$get_id())
  expect_silent(region$configure_paths(permeability = perm))
  expect_equal(perm$get_id(), 1)
  expect_silent(region$calculate_paths(5922))
  paths1 <- region$get_paths(5922, max_distance = 40000)
  paths2 <- region$get_paths(5922, max_distance = 40000, perm_id = 1)
  expect_true(length(paths2$idx$`5922`$cell) < length(paths1$idx$`5922`$cell))
  expect_true(length(paths2$idx$`5922`$aggr) < length(paths1$idx$`5922`$aggr))
  expect_true(all(paths2$perm_dist$`5922`$cell > paths2$distances$`5922`$cell))
  expect_true(all(paths2$perm_dist$`5922`$aggr > paths2$distances$`5922`$aggr))
  expect_true(all(paths2$perm_dist$`5922`$cell <= 40000))
  expect_true(all(paths2$perm_dist$`5922`$aggr <= 40000))
})

test_that("creates paths to city patches", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
  expect_silent(region$configure_paths(directions = TRUE,
                                       max_distance = 150000))
  expect_silent(region$calculate_paths(1))
  expect_silent(paths <- region$get_paths(1, directions = TRUE))
  expect_equal(locations$name[paths$idx$`1`],
               c("Geelong", "Ballarat", "Bendigo"))
  expect_true(all(paths$distances$`1` <= 150000))
  expect_true(all(paths$directions$`1` <= 360))
})

test_that("creates permeability paths to city patches", {
  TEST_DIRECTORY <- test_path("test_inputs")
  locations <- utils::read.csv(file.path(TEST_DIRECTORY, "vic_cities.csv"))
  region <- Region(locations)
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
  perm <- Permeability(perm_data, region)
  expect_silent(region$configure_paths(permeability = perm))
  expect_silent(region$calculate_paths(1:2))
  paths1 <- region$get_paths(1:2, max_distance = 400000)
  paths2 <- region$get_paths(1, max_distance = 400000, perm_id = 1)
  expect_true(length(paths2$idx$`1`) < length(paths1$idx$`1`))
  expect_true(all(paths2$perm_dist$`1` > paths2$distances$`1`))
  geelong <- which(paths2$idx$`1` == 2)
  expect_true(paths2$perm_dist$`1`[geelong] ==
                round(paths2$distances$`1`[geelong]/0.8))
  warrn1 <- which(paths1$idx$`2` == 5)
  warrn2 <- which(paths2$idx$`1` == 5)
  expect_true(paths2$perm_dist$`1`[warrn1] ==
                (round(paths2$distances$`1`[geelong]/0.8) +
                   round(paths1$distances$`2`[warrn2]/0.6)))
})
