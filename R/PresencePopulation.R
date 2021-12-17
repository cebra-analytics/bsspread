#' Presence-only population model class builder
#'
#' Builds a class for representing presence-only populations in spread
#' simulations. Growth is implicitly simulated via a delay in the ability of
#' newly established populations to spread.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for spread simulations.
#' @param spread_delay Number of simulation time steps before an incursion at
#'   a newly established location can spread to other locations. This provides
#'   an implicit mechanism for growth. Default \code{NULL} assumes no delay.
#' @param ... Additional parameters.
#' @return A \code{PresencePopulation} class object (list) containing functions
#'   for accessing attributes and simulating (implicit) growth:
#'   \describe{
#'     \item{\code{get_type()}}{Get the population representation type.}
#'     \item{\code{get_spread_delay()}}{Get the delay in spread ability.}
#'     \item{\code{grow(x)}}{Performs implicit growth on the population
#'       \code{x} vector via spread delay when present, and returns the
#'       transformed vector.}
#'   }
#' @include Population.R
#' @export
PresencePopulation <- function(region,
                               spread_delay = NULL, ...) {

  # Build via base class
  self <- Population(region,
                     type = "presence_only",
                     class = "PresencePopulation")

  # Setup spread delay counter for each location
  if (is.numeric(spread_delay)) {
    time_since_established <- vector("integer", region$get_locations())*NA
  }

  # Get spread delay
  self$get_spread_delay <- function() {
    return(spread_delay)
  }

  # Grow method - override for implicit delay-based growth when present
  if (is.numeric(spread_delay)) {
    self$grow <- function(x) {
      # TODO
      print("presence only")
      print(head(time_since_established))
      return(x) # no change
    }
  }

  return(self)
}
