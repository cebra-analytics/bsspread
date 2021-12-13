#' Incursions class builder
#'
#' Builds a class to generate invasive species incursions, either by selecting
#' a single incursion location via stochastic sampling based on weightings for
#' each available spatial location, or by sampling potentially multiple
#' locations via incursion or arrival probabilities defined for each location,
#' which optionally may continue to generate incursions at specified time
#' intervals during a spread simulation.
#'
#' @param x A \code{raster::RasterLayer} or \code{terra::SpatRaster}
#'   object representing the incursion weightings or arrival probabilities at
#'   each spatial location.
#' @param type One of \code{"weight"} or \code{"prob"} to indicate if the
#'   values in \code{x} represent single incursion weightings or arrival
#'   probabilities respectively.
#' @param continued Logical to indicate whether incursions based on arrival
#'   probabilities (when \code{type = "prob"}) continue to be generated during
#'   spread simulations.
#' @param time_steps The time interval in simulation steps for generating
#'   subsequent incursions when \code{type = "prob"} and
#'   \code{continued = TRUE}.
#' @param ... Additional parameters.
#' @return An \code{Incursions} class object (list) containing functions for
#'   accessing attributes (of the function environment) TODO:
#'   \describe{
#'     \item{\code{TODO()}}{TODO.}
#'   }
#' @export
Incursions <- function(x,
                       type = c("weight", "prob"),
                       continued = FALSE,
                       time_steps = 1, ...) {
  UseMethod("Incursions")
}

#' @name Incursions
#' @export
Incursions.Raster <- function(x, ...) {
  # Call the terra version of the function
  Incursions(terra::rast(x), ...)
}

#' @name Incursions
#' @export
Incursions.SpatRaster <- function(x,
                                  type = c("weight", "prob"),
                                  continued = FALSE,
                                  time_steps = 1, ...) {
  type <- match.arg(type)

  # Create a class structure
  self <- structure(list(), class = "Incursions")

  return(self)
}
