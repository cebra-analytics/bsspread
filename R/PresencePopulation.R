#' Presence-only population model class builder
#'
#' Builds a class for representing presence-only populations in spread
#' simulations. Growth is implicitly simulated via a delay in the ability of
#' newly established populations to spread.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the spread simulations.
#' @param establish_pr An optional (static) vector or matrix (containing
#'   temporal columns) of probability values (0-1) to represent the likelihood
#'   of establishment at each location (row) specified by the \code{region}.
#'   This may be used to avoid transient/unsuccessful incursion or migration
#'   arrivals from being presented in the simulation results, and/or from
#'   subsequently contributing to spread. Default is \code{NULL}. The number
#'   of columns for temporal capacity should coincide with the number of
#'   simulation time steps, or a cyclic pattern (e.g. 12 columns for seasonal
#'   variation with monthly time steps).
#' @param spread_delay Number of simulation time steps before an incursion at
#'   a newly established location can spread to other locations. This provides
#'   an implicit mechanism for growth. Default \code{NULL} assumes no delay.
#' @param ... Additional parameters.
#' @return A \code{PresencePopulation} class object (list) containing functions
#'   for accessing attributes and simulating (implicit) growth:
#'   \describe{
#'     \item{\code{get_region()}}{Get the region object.}
#'     \item{\code{get_type()}}{Get the population representation type.}
#'     \item{\code{get_establish_pr(cells = NULL, tm = NULL)}}{Get the
#'       establishment probability as a vector of values for each region
#'       location or optionally specified region locations \code{cells}
#'       (indices) at (optional) simulation time step \code{tm} (for
#'       temporally defined establishment probability).}
#'     \item{\code{get_spread_delay()}}{Get the delay in spread ability.}
#'     \item{\code{make(initial, current, incursion, tm)}}{Make a population
#'       vector using vectors of the \code{initial} or \code{current} and
#'       \code{incursion} population at each region location at simulation time
#'       step \code{tm}.}
#'     \item{\code{grow(x, tm)}}{Performs implicit growth on the population
#'       \code{x} vector at simulation time step \code{tm} via spread delay
#'       when present, and returns the transformed vector.}
#'   }
#' @references
#'   Bradhurst, R., Spring, D., Stanaway, M., Milner, J., & Kompas, T. (2021).
#'   A generalised and scalable framework for modelling incursions,
#'   surveillance and control of plant and environmental pests.
#'   \emph{Environmental Modelling & Software}, 139, N.PAG.
#'   \doi{10.1016/j.envsoft.2021.105004}
#'
#'   Jongejans, E., Skarpaas, O., & Shea, K. (2008). Dispersal, demography and
#'   spatial population models for conservation and control management.
#'   \emph{Perspectives In Plant Ecology Evolution And Systematics}, 9(3–4),
#'   153–170. \doi{10.1016/j.ppees.2007.09.005}
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
  self$make <- function(initial = NULL, current = NULL, incursion = NULL,
                        tm = NULL) {

    # Convert to initial values to logical
    if (!is.null(initial)) {
      initial_age <- attr(initial, "age")
      initial <- as.logical(initial)
      if (!is.null(initial_age)) { # re-attach
        attr(initial, "age") <- initial_age + initial*0 # as vector
      }
    }

    # Run inherited function
    values <- inherited_make(initial = initial, current = current,
                             incursion = incursion, tm = tm)

    # Set attributes for type and spread delay
    if (is.null(current)) {
      attr(values, "type") <- self$get_type()
      if (is.numeric(spread_delay)) {
        attr(values, "spread_delay") <- rep(NA, region$get_locations())
        if (!is.null(attr(values, "age"))) {
          attr(values, "spread_delay")[which(values)] <-
            pmax(0, spread_delay - attr(values, "age")[which(values)] + 1)
        }
      }
    } else {
      attributes(values) <- attributes(current)
    }

    return(values)
  }

  # Grow method - override for implicit (spread) delay-based growth
  if (is.numeric(spread_delay)) {
    self$grow <- function(x, tm) {

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
