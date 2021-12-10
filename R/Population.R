#' Population model class builder
#'
#' Builds a class for population representation and growth functionality for
#' the spread simulations.
#'
#' @param region A \code{raster::RasterLayer}, \code{terra::SpatRaster}, or
#'   \code{Region} or inherited class object representing the spatial region
#'   (template) for spread simulations.
#' @param type One of \code{"presence_only"} (default), \code{"unstructured"},
#'   or \code{"stage_structured"} to indicate how populations are represented.
#' @param growth Numeric intrinsic growth rate or lambda (e.g. 1.2 for 20%
#'   growth per time step) when \code{type = "unstructured"}, or
#'   age/stage-based transition rates (matrix) for each time step when
#'   \code{type = "stage_structured"}. Default is \code{NULL} for when
#'   \code{type = "presence_only"}.
#' @param capacity A data frame, \code{raster::RasterLayer}, or
#'   \code{terra::SpatRaster} with values for the carrying capacity of the
#'   invasive species at each location. Default is \code{NULL} for when
#'   \code{type = "presence_only"}.
#' @param ... Additional parameters.
#' @return An \code{Population} class object (list) containing functions for
#'   accessing attributes (of the function environment) TODO:
#'   \describe{
#'     \item{\code{TODO()}}{TODO.}
#'   }
#' @include Region.R
#' @export
Population <- function(region,
                       type = c("presence_only", "unstructured",
                                "stage_structured"),
                       growth = NULL,
                       capacity = NULL, ...) {
  UseMethod("Population")
}

#' @name Population
#' @export
Population.Raster <- function(region, ...) {
  # Call Region class version
  Population(Region(region), ...)
}

#' @name Population
#' @export
Population.SpatRaster <- function(region, ...) {
  # Call Region class version
  Population(Region(region), ...)
}

#' @name Population
#' @export
Population.Region <- function(region,
                              type = c("presence_only", "unstructured",
                                       "stage_structured"),
                              growth = NULL,
                              capacity = NULL, ...) {

  # Create a class structure
  self <- structure(list(), class = "Population")

  return(self)
}
