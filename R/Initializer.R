#' Initializer class builder
#'
#' Builds a class to initialize a spread simulation, given a specified spatial
#' region and population representation, via either a predefined initial
#' population distribution, or via stochastic generation of invasive species
#' incursions.
#'
#' @param x A \code{raster::RasterLayer} or \code{terra::SpatRaster} defining
#'   the initial population distribution, or an \code{Incursions} or inherited
#'   class object for generating initial, and optionally continued, invasive
#'   species incursions.
#' @param region A \code{Region} or inherited class object defining the spatial
#'   locations included in the spread simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation and growth functionality for the
#'   spread simulations.
#' @param ... Additional parameters.
#' @return An \code{Initializer} class object (list) containing functions for
#'   accessing attributes (of the function environment) TODO:
#'   \describe{
#'     \item{\code{TODO()}}{TODO.}
#'   }
#' @export
Initializer <- function(x,
                        region = NULL,
                        population_model = NULL, ...) {
  UseMethod("Initializer")
}

#' @name Initializer
#' @export
Initializer.Raster <- function(x, ...) {
  # Call the terra version of the function
  Initializer(terra::rast(x), ...)
}

#' @name Initializer
#' @export
Initializer.SpatRaster <- function(x,
                                   region = NULL,
                                   population_model = NULL, ...) {

  # Create a class structure
  self <- structure(list(), class = "Initializer")

  return(self)
}

#' @name Initializer
#' @export
Initializer.Incursions <- function(x,
                                   region = NULL,
                                   population_model = NULL, ...) {

  # Create a class structure
  self <- structure(list(), class = "Initializer")

  return(self)
}
