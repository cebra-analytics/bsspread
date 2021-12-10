#' Dispersal model abstract class builder
#'
#' Builds a class for dispersal representation and functionality for the
#' spread simulations.
#'
#' @param region A \code{raster::RasterLayer}, \code{terra::SpatRaster}, or
#'   \code{Region} or inherited class object representing the spatial region
#'   (template) for spread simulations.
#' @param ... Additional parameters.
#' @return An \code{Dispersal} class object (list) containing functions for
#'   accessing attributes (of the function environment) TODO:
#'   \describe{
#'     \item{\code{TODO()}}{TODO.}
#'   }
#' @include Region.R
#' @export
Dispersal <- function(region, ...) {
  UseMethod("Dispersal")
}

#' @name Dispersal
#' @export
Dispersal.Raster <- function(region, ...) {
  # Call Region class version
  Dispersal(Region(region), ...)
}

#' @name Dispersal
#' @export
Dispersal.SpatRaster <- function(region, ...) {
  # Call Region class version
  Dispersal(Region(region), ...)
}

#' @name Dispersal
#' @export
Dispersal.Region <- function(region, ...) {

  # Create a class structure
  self <- structure(list(), class = "Dispersal")

  return(self)
}
