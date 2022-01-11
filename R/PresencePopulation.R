#' Presence-only population model class builder
#'
#' Builds a class for representing presence-only populations in spread
#' simulations. Growth is implicitly simulated via a delay in the ability of
#' newly established populations to spread.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the spread simulations.
#' @param establish_pr An optional vector of probability values (0-1) to
#'   represent the likelihood of establishment at each location specified by
#'   the \code{region}. This may be used to avoid transient/unsuccessful
#'   incursions or migrations from being presented in the simulation results,
#'   and/or from subsequently contributing to spread in presence-only models.
#'   Default is \code{NULL}.
#' @param spread_delay Number of simulation time steps before an incursion at
#'   a newly established location can spread to other locations. This provides
#'   an implicit mechanism for growth. Default \code{NULL} assumes no delay.
#' @param ... Additional parameters.
#' @return A \code{PresencePopulation} class object (list) containing functions
#'   for accessing attributes and simulating (implicit) growth:
#'   \describe{
#'     \item{\code{get_type()}}{Get the population representation type.}
#'     \item{\code{get_establish_pr()}}{Get the establishment probability as a
#'       vector of values for each location.}
#'     \item{\code{get_spread_delay()}}{Get the delay in spread ability.}
#'     \item{\code{make(initial, current, incursion)}}{Make a population vector
#'       using vectors of the \code{initial} or \code{current} and
#'       \code{incursion} population at each region location.}
#'     \item{\code{grow(x)}}{Performs implicit growth on the population
#'       \code{x} vector via spread delay when present, and returns the
#'       transformed vector.}
#'   }
#' @include Population.R
#' @export
PresencePopulation <- function(region,
                               establish_pr = NULL,
                               spread_delay = NULL, ...) {

  # Build via base class
  self <- Population(region,
                     type = "presence_only",
                     establish_pr = establish_pr,
                     class = "PresencePopulation")

  # Get spread delay
  self$get_spread_delay <- function() {
    return(spread_delay)
  }

  # Make method (extends inherited function from Population class)
  inherited_make <- self$make
  self$make <- function(initial = NULL, current = NULL, incursion = NULL) {

    # Run inherited function
    values <- inherited_make(initial = initial, current = current,
                            incursion = incursion)

    # Set attributes for type and spread delay
    if (is.null(current)) {
      attributes(values) <- list(type = self$get_type(),
                                spread_delay = NULL)
      if (is.numeric(spread_delay)) {
        attr(values, "spread_delay") <- rep(NA, region$get_locations())
      }
    } else {
      attributes(values) <- attributes(current)
    }

    return(values)
  }

  # Grow method - override for implicit (spread) delay-based growth
  if (is.numeric(spread_delay)) {
    self$grow <- function(x) {

      # Decrement spread delay on previously occupied cells
      idx <- which(attr(x, "spread_delay") > 0)
      attr(x, "spread_delay")[idx] <- attr(x, "spread_delay")[idx] - 1

      # Add spread delay to newly occupied positions
      if (is.numeric(spread_delay)) {
        attr(x, "spread_delay")[which(
          x & is.na(attr(x, "spread_delay")))] <- spread_delay
      }

      return(x)
    }
  }

  return(self)
}
