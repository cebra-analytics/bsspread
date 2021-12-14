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
#'     \item{\code{get_template()}}{Get the spatial template (with zeros in
#'       non-NA locations).}
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

  # Get spatial template
  self$get_template <- function() {
    template <- x*0
    names(template) <- "value"
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

  return(self)
}
