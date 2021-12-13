#' Dispersal model abstract class builder
#'
#' Builds a class for dispersal representation and functionality for the
#' spread simulations.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for spread simulations.
#' @param class Character class name for inherited classes. Default is empty.
#' @param ... Additional parameters.
#' @return A \code{Dispersal} class object (list) containing functions for
#'   accessing attributes (of the function environment) and performing
#'   dispersal:
#'   \describe{
#'     \item{\code{get_region()}}{Get the spatial region object.}
#'     \item{\code{disperse(x)}}{Perform location dispersal on population
#'       \code{x} array (stages by locations), and return transformed array.}
#'     \item{\code{TODO()}}{TODO.}
#'   }
#' @include Region.R
#' @export
Dispersal <- function(region,
                      class = character(), ...) {
  UseMethod("Dispersal")
}

#' @name Dispersal
#' @export
Dispersal.Region <- function(region,
                             class = character(), ...) {

  # Create a class structure
  self <- structure(list(), class = c(class, "Dispersal"))

  # Get region object
  self$get_region <- function() {
    return(region)
  }

  # Generic disperse method (overridden in subclass)
  self$disperse <- function(x) return(x) # no change

  return(self)
}
