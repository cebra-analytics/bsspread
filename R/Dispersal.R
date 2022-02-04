#' Dispersal model abstract class builder
#'
#' Builds a class for dispersal representation and functionality for spread
#' simulations.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for spread simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the spread simulations.
#' @param class Character class name for inherited classes. Default is empty.
#' @param ... Additional parameters.
#' @return A \code{Dispersal} class object (list) containing functions for
#'   accessing attributes (of the function environment) and performing
#'   dispersal:
#'   \describe{
#'     \item{\code{get_region()}}{Get the spatial region object.}
#'     \item{\code{get_population_model()}}{Get the population model object.}
#'     \item{\code{disperse(x)}}{Perform location dispersal on a list \code{x}
#'       of population vectors (for each location) or arrays (locations by
#'       stages), representing the \code{original}, \code{remaining} and
#'       \code{relocated} populations, and return the transformed list of
#'       vectors or arrays. The separation of original, remaining and relocated
#'       populations enables multiple models for different dispersal vectors to
#'       run in sequence.}
#'     \item{\code{TODO()}}{TODO.}
#'   }
#' @include Region.R
#' @export
Dispersal <- function(region, population_model,
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

  # Get population model
  self$get_population_model <- function() {
    return(population_model)
  }

  # Generic disperse method (overridden in subclass)
  self$disperse <- function(x) return(x) # no change

  return(self)
}
