#' Attractor class builder
#'
#' Builds a class to represent spatially distributed weights, or actual or
#' relative probabilities assigned to cells/locations to represent an
#' "attractor" that influences the likelihood of dispersal from or to other
#' locations.
#'
#' @param x A \code{raster::RasterLayer} or \code{terra::SpatRaster}
#'   object representing spatial attractor values across the simulation region.
#' @param region A \code{Region} or inherited class object defining the spatial
#'   locations included in the spread simulations.
#' @param ... Additional parameters.
#' @return A \code{Attractor} class object (list) containing functions for
#'   accessing attributes:
#'   \describe{
#'     \item{\code{get_rast()}}{Get the attractor \code{terra::SpatRaster}
#'       object.}
#'     \item{\code{get_values(cells = NULL)}}{Get the attractor values for
#'       optionally specified region (non-NA) cell indices. Default is all.}
#'   }
#' @include Region.R
#' @export
Attractor <- function(x, region, ...) {
  UseMethod("Attractor")
}

#' @name Attractor
#' @export
Attractor.Raster <- function(x, ...) {
  # Call the terra version of the function
  Attractor(terra::rast(x), ...)
}

#' @name Attractor
#' @export
Attractor.SpatRaster <- function(x, region, ...) {

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
  self <- structure(list(), class = "Attractor")

  # Get the attractor raster
  self$get_rast <- function() {
    return(x)
  }

  # Get the values for specified region (non-NA) cell indices
  self$get_values <- function(cells = NULL) {
    if (is.numeric(cells)) {
      return(x[region$get_indices()[cells]])
    } else {
      return(x[region$get_indices()])
    }
  }

  return(self)
}
