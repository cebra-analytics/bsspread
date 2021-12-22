#' Dispersal model abstract class builder
#'
#' Builds a class for dispersal representation and functionality for spread
#' simulations.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for spread simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the spread simulations.
#' @param establish_pr An optional vector of probability values (0-1) to
#'   represent the likelihood of establishment at each location specified by
#'   the \code{region}. This may be used to avoid transient/unsuccessful
#'   migrations from being presented in the simulation results, and/or from
#'   subsequently contributing to spread in presence-only models. Default is
#'   \code{NULL}.
#' @param class Character class name for inherited classes. Default is empty.
#' @param ... Additional parameters.
#' @return A \code{Dispersal} class object (list) containing functions for
#'   accessing attributes (of the function environment) and performing
#'   dispersal:
#'   \describe{
#'     \item{\code{get_region()}}{Get the spatial region object.}
#'     \item{\code{get_population_model()}}{Get the population model object.}
#'     \item{\code{get_establish_pr()}}{Get the establishment probability
#'       vector.}
#'     \item{\code{disperse(x)}}{Perform location dispersal on population
#'       \code{x} vector (for each location) or array (locations by stages),
#'       and return transformed vector or array.}
#'     \item{\code{TODO()}}{TODO.}
#'   }
#' @include Region.R
#' @export
Dispersal <- function(region, population_model,
                      establish_pr = NULL,
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

  # Get establishment probability
  self$get_establish_pr <- function() {
    return(establish_pr)
  }

  # Generic disperse method (overridden in subclass)
  self$disperse <- function(x) return(x) # no change

  return(self)
}
