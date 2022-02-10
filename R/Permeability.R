#' Permeability class builder
#'
#' Builds a class to represent the permeability of cells/locations within the
#' spatial region of a spread simulation. The inverse of the permeability (0-1)
#' of cells is used to (linearly) scale the actual distance between adjacent
#' cells. For example, a permeability value of \code{0.5} results in an
#' effective distance twice that of the actual distance, a value of \code{0}
#' prevents spread to or through a cell, and a value of \code{1} does not
#' modify the effective distance. The class may also be used to calculate an
#' optional aggregate (courser resolution) version of the permeability layer.
#'
#' @param x A \code{raster::RasterLayer} or \code{terra::SpatRaster}
#'   object representing the spatial permeability of the spread simulation
#'   region.
#' @param region A \code{Region} or inherited class object defining the spatial
#'   locations included in the spread simulations.
#' @param ... Additional parameters.
#' @return A \code{Permeability} class object (list) containing functions for
#'   accessing attributes:
#'   \describe{
#'     \item{\code{get_id()}}{Get the object numeric identifier.}
#'     \item{\code{set_id(id)}}{Set the object numeric identifier.}
#'     \item{\code{get_rast()}}{Get the permeability \code{terra::SpatRaster}
#'       object.}
#'     \item{\code{get_aggr_rast()}}{Get the aggregate permeability
#'       \code{terra::SpatRaster} object.}
#'     \item{\code{calculate_aggregate(aggr_factor)}}{Calculate the aggregate
#'       permeability \code{terra::SpatRaster} object via a specified
#'       aggregation factor.}
#'     \item{\code{has_aggregate()}}{Check the (logical) presence of an
#'       aggregate permeability.}
#'   }
#' @include Region.R
#' @export
Permeability <- function(x, region, ...) {
  UseMethod("Permeability")
}

#' @name Permeability
#' @export
Permeability.Raster <- function(x, ...) {
  # Call the terra version of the function
  Permeability(terra::rast(x), ...)
}

#' @name Permeability
#' @export
Permeability.SpatRaster <- function(x, region, ...) {

  # Check region
  if (!inherits(region, "Region")) {
    stop("Region model must be a 'Region' or inherited class object.",
         call. = FALSE)
  }
  if (!region$is_compatible(x)) {
    stop("The spatial object x should be compatible with that defining the ",
         "region.", call. = FALSE)
  }

  # Create a class structure
  self <- structure(list(), class = "Permeability")

  # Permeability object id
  id <- NULL

  # Aggregate permeability raster
  aggr_rast <- NULL

  # Get the permeability object id
  self$get_id <- function() {
    return(x)
  }

  # Set the permeability object id
  self$set_id <- function(id) {
    id <<- id
  }

  # Get the permeability raster
  self$get_rast <- function() {
    return(x)
  }

  # Check the presence of an aggregate permeability
  self$has_aggregate <- function() {
    return(!is.null(aggr_rast))
  }

  # Get the aggregate permeability raster
  self$get_aggr_rast <- function() {
    return(aggr_rast)
  }

  # Calculate the aggregate permeability
  self$calculate_aggregate <- function(aggr_factor) {
    aggr_rast <<- terra::aggregate(x, fact = aggr_factor, fun = "mean",
                                   na.rm = TRUE)
  }

  return(self)
}
