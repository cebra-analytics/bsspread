#' Region class builder
#'
#' Builds a class to represent a spatial region for a spread simulation,
#' defined by a raster layer with active (non-NA) and/or suitable (non-zero)
#' cells, or by a network of locations or patches.
#'
#' @param x A \code{raster::RasterLayer} or \code{terra::SpatRaster}
#'   object representing the spatial region (template) for spread simulations.
#' @param suitability Logical to indicate whether only suitable locations
#'   (cells) with values greater than zero should be used in the simulations.
#' @param ... Additional parameters.
#' @return A \code{Region} class object (list) containing functions for
#'   accessing attributes and checking compatibility of objects with the
#'   region:
#'   \describe{
#'     \item{\code{get_template(empty = FALSE)}}{Get the spatial template with
#'      either zeros in non-NA locations (default), or with no values when
#'      \code{empty = TRUE}.}
#'     \item{\code{get_indices()}}{Get cell indices of grid locations that are
#'       included in the simulation.}
#'     \item{\code{get_locations()}}{Get the number of locations (cells or
#'       patches) that are included in the simulation.}
#'     \item{\code{get_type()}}{Get the type of spatial representation: grid
#'       (raster cells) or patch (network).}
#'     \item{\code{is_compatible(y)}}{Check the compatibility of object
#'       \code{y} with the region defined by \code{x}.}
#'   }
#' @export
Region <- function(x, suitability = FALSE, ...) {
  UseMethod("Region")
}

#' @name Region
#' @export
Region.Raster <- function(x, ...) {
  # Call the terra version of the function
  Region(terra::rast(x), ...)
}

#' @name Region
#' @export
Region.SpatRaster <- function(x, suitability = FALSE, ...) {

  # Extract non-NA cell data
  region_df <- terra::as.data.frame(x, cells = TRUE, na.rm = TRUE)

  # Set environment variables
  if (suitability) { # non-zero only
    region_df <- dplyr::filter(region_df, region_df[,2] > 0)
  }
  indices <- region_df$cell

  # Remove cell data from the function environment
  rm(region_df)

  # Create a class structure
  self <- structure(list(), class = "Region")

  # Get spatial template with zero/NA or empty values
  self$get_template <- function(empty = FALSE) {
    if (empty) {
      template <- terra::rast(x)
    } else {
      template <- x*0
      names(template) <- "value"
    }
    return(template)
  }

  # Get cell indices
  self$get_indices <- function() {
    return(indices)
  }

  # Get the number of active cell locations
  self$get_locations <- function() {
    return(length(indices))
  }

  # Get the spatial region type
  self$get_type <- function() {
    return("grid")
  }

  # Check compatibility of a spatial raster y with the region defined by x
  self$is_compatible <- function(y) {
    y <- terra::rast(y)
    return(terra::crs(y) == terra::crs(x) &&
             terra::ext(y) == terra::ext(x) &&
             all(terra::res(y) == terra::res(x)))
  }

  # Check, get, and set long-distance dispersal aggregation (list)
  aggr <- NULL
  self$use_aggr <- function() {
    return(!is.null(aggr))
  }
  self$get_aggr <- function() { # ???
    return(aggr)
  }
  self$set_aggr <- function(aggr_factor, inner_radius) {
    aggr <<- list(factor = aggr_factor, inner_radius = inner_radius)
    aggr$rast <<- terra::aggregate(self$get_template(), fact = aggr_factor,
                                  fun = "mean", na.rm = TRUE)
  }

  # Lists of distances, directions, and graphs for dispersal calculations
  distances <- list()
  directions <- list()
  graphs <- list()

  # Other variables dynamically created for dispersal calculations
  region_pts <- NULL
  any_dispersal_variables <- function() {
    return(!is.null(region_pts))
  }
  set_dispersal_variables <- function() {
    region_pts <<- terra::as.points(self$get_template())
  }

  # Get distances for a region cell
  self$get_distances <- function(cell_i, max_distance = NULL) {

    # Set dispersal variables when required
    if (!any_dispersal_variables()) {
      set_dispersal_variables()
    }

    # Resolve all reachable (local) cells
    if (self$use_aggr()) {
      # TODO

    } else {
      if (!is.null(max_distance)) {
        inner_vect <- terra::buffer(region_pts[cell_i,], width = max_distance,
                                    quadsegs = 180)
        cell_idx <- terra::cells(x, inner_vect, touches = TRUE)[, "cell"]
        cell_idx <- which(indices %in% cell_idx)
      } else {
        # TODO
      }
    }

    # IN PROGRESS

    distance_list[[as.character(cell_i)]] <- list(cell = as.list(
      as.integer(round(as.numeric(
        terra::distance(region_pts[cell_i], region_pts[cell_idx]))))))

    distance_list[[as.character(cell_i)]][["aggr"]] <- as.list(
      as.integer(round(as.numeric(
        terra::distance(region_pts[cell_i], aggr_region_pts[aggr_idx])))))
  }

  # Get directions
  self$get_directions()

  # Get graph (shortest path) weights
  self$get_graph_weights()

  # IN PROGRESS

  return(self)
}
