#' Unstructured population model class builder
#'
#' Builds a class for representing unstructured populations in spread
#' simulations. Simulates logistic (capacity-limited) growth.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the spread simulations.
#' @param growth Numeric intrinsic growth rate or lambda (e.g. 1.2 for 20%
#'   growth per time step). Default is \code{1.0} for no increase.
#' @param capacity A vector of carrying capacity values of the invasive species
#'   at each location specified by the \code{region}. Default is \code{NULL}
#'   for when growth is not capacity-limited.
#' @param establish_pr An optional vector of probability values (0-1) to
#'   represent the likelihood of establishment at each location specified by
#'   the \code{region}. This may be used to avoid transient/unsuccessful
#'   incursions or migrations from being presented in the simulation results.
#'   Default is \code{NULL}.
#' @param incursion_values Defines how incursion locations are populated.
#'   Either via a fixed single value or vector of values for each location,
#'   or via a list specifying the range of values via \code{min} and \code{max}
#'   (single values or vectors) for uniform random sampling of values.
#' @param ... Additional parameters.
#' @return A \code{UnstructPopulation} class object (list) containing functions
#'   for accessing attributes and simulating growth:
#'   \describe{
#'     \item{\code{get_type()}}{Get the population representation type.}
#'     \item{\code{get_growth()}}{Get the growth rate.}
#'     \item{\code{get_capacity()}}{Get the carrying capacity as a vector of
#'       values for each location.}
#'     \item{\code{get_establish_pr()}}{Get the establishment probability as a
#'       vector of values for each location.}
#'     \item{\code{make(initial, current, incursion)}}{Make a population vector
#'       via using vectors of the \code{initial} or \code{current} and
#'       \code{incursion} population at each region location.}
#'     \item{\code{grow(x)}}{Performs logistic (capacity-limited) growth on the
#'       population \code{x} vector, and returns the transformed vector.}
#'   }
#' @include Population.R
#' @export
UnstructPopulation <- function(region,
                               growth = 1.0,
                               capacity = NULL,
                               establish_pr = NULL,
                               incursion_values = NULL, ...) {

  # Build via base class
  self <- Population(region,
                     type = "unstructured",
                     growth = growth,
                     capacity = capacity,
                     establish_pr = establish_pr,
                     incursion_values = incursion_values,
                     class = "UnstructPopulation")

  # Grow method - override for logistic (capacity-limited) growth
  self$grow <- function(x) {

    # Indices of occupied locations
    indices <- which(x > 0)

    # Calculate logistic growth rates
    if (is.numeric(capacity)) { # capacity-limited

      # Remove populations at locations having zero capacity
      if (any(capacity[indices] <= 0)) {
        zero_idx <- indices[which(capacity[indices] <= 0)]
        x[zero_idx] <- 0
        indices <- indices[!indices %in% zero_idx]
      }

      # Calculate capacity-limited growth rates for each occupied location
      r <- exp(log(growth)*(1 - x[indices]/capacity[indices]))

    } else {
      r <- growth
    }

    # Sample the new population values via the Poisson distribution
    if (length(indices)) {
      x[indices] <- stats::rpois(length(indices), r*x[indices])
    }

    return(x)
  }

  return(self)
}