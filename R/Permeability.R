#' Permeability class builder
#'
#' Builds a class to represent the permeability of cells within a grid-based
#' spatial region, or paths between patches in a patch-based region, of a
#' spread simulation. The inverse of the permeability (0-1) of cells/paths is
#' used to (linearly) scale the actual distance between adjacent cells/patches.
#' For example, a permeability value of \code{0.5} results in an effective
#' distance twice that of the actual distance, a value of \code{0} prevents
#' spread to or through a cell, or spread along a network path to a patch, and
#' a value of \code{1} does not modify the effective distance.
#'
#' @param x A \code{raster::RasterLayer} or \code{terra::SpatRaster}
#'   object representing the grid-based spatial permeability, or a 3-column
#'   data frame (or matrix) of patch index pairs and permeability values
#'   (weights) for connected (adjacent) paths between patches, consistent with
#'   the spread simulation region.
#' @param region A \code{Region} or inherited class object defining the spatial
#'   locations included in the spread simulations.
#' @param ... Additional parameters.
#' @return A \code{Permeability} class object (list) containing functions for
#'   accessing attributes:
#'   \describe{
#'     \item{\code{get_id()}}{Get the object numeric identifier.}
#'     \item{\code{set_id(id)}}{Set the object numeric identifier.}
#'     \item{\code{get_rast()}}{Get the permeability \code{terra::SpatRaster}
#'       object when the simulation region is grid-based.}
#'     \item{\code{get_data()}}{Get the path data specifying the weights
#'       between patches when the simulation region is patch-based.}
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

  # Get class structure from common internal function
  self <- Permeability_internal(x, region, ...)

  # Get the permeability raster
  self$get_rast <- function() {
    return(x)
  }

  return(self)
}

#' @name Permeability
#' @export
Permeability.matrix <- function(x, ...) {

  # Call the data frame version of the function
  Permeability(as.data.frame(x), ...)
}

#' @name Permeability
#' @export
Permeability.data.frame <- function(x, region, ...) {

  # Get class structure from common internal function
  self <- Permeability_internal(x, region, ...)

  # Get the permeability data
  self$get_data <- function() {
    return(x)
  }

  return(self)
}

Permeability_internal <- function(x, region, ...) {

  # Check region
  if (!inherits(region, "Region")) {
    stop("Region model must be a 'Region' or inherited class object.",
         call. = FALSE)
  }
  if (!region$is_compatible(x)) {
    stop("The spatial object x should be compatible with that defining the ",
         "region.", call. = FALSE)
  }
  rm(x)

  # Create a class structure
  self <- structure(list(), class = "Permeability")

  # Permeability object id
  id <- NULL

  # Get the permeability object id
  self$get_id <- function() {
    return(id)
  }

  # Set the permeability object id
  self$set_id <- function(id) {
    id <<- id
  }

  return(self)
}

