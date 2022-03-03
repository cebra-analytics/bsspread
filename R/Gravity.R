#' Gravity dispersal model class builder
#'
#' Builds a class for representing gravity dispersal for spread simulations
#' with functionality inherited from the generic dispersal model.
#'
#' Gravity dispersal simulates the spread of populations between locations
#' based on attributes or "attractors" at each location, and the distance
#' between source and destination locations. This approach is analogous to
#' Newton's law of gravitation and has the general form for calculating the
#' (relative) probability of dispersal:
#' \code{dispersal = f(attractors)/distance^beta}, where \code{f} is a function
#' of attractor variables (typically products) and \code{beta} is constant
#' (often 1). Gravity dispersal may be simulated for presence-only,
#' unstructured or stage-based \code{populations}. For unstructured or
#' stage-based populations, a specified \code{proportion} of the population at
#' each occupied (cell) location is selected (sampled) for dispersal at each
#' simulation time step. Presence-only populations may disperse via specifying
#' a (mean) number of dispersal \code{events}. Dispersal events are generated
#' for each occupied location. Unstructured and stage-based population
#' dispersal may also utilize specified events, otherwise an event is assigned
#' to each dispersing individual. Presence-only population dispersal may also
#' be configured without \code{events} when \code{proportion} is specified to
#' represent a scaling multiplier for (presumed) actual (rather than relative)
#' dispersal probabilities. Under these circumstances destination locations are
#' sampled from all reachable destinations using the scaled probabilities. The
#' probability of gravity dispersal may optionally be adjusted via
#' \code{direction} functions. An optional establishment likelihood (layer),
#' which is configured via the population model, may be applied to each
#' dispersal, resulting in potential "deaths" of individuals or unsuccessful
#' presence-only dispersal events. The dispersal functionality utilizes a
#' wrapped population list of separate values for the \code{original},
#' \code{remaining}, and \code{relocated} populations. This separation enables
#' multiple dispersal models, representing different dispersal vectors, to be
#' run in sequence.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for spread simulations. The region object
#'   contains functionality for calculating path distances and directions,
#'   permeability graphs, and structures to facilitate two-tier dispersal.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the spread simulations.
#' @param attractors List containing \code{Attractor} (or inherited) class
#'   objects for spatially weighted dispersal.
#' @param attractor_function An optional function of form
#'   \code{function(attractors)} for combining attractors specified as a list
#'   of \code{Attractor} class objects and returning a list of transformed
#'   \code{Attractor} class objects. Default is \code{NULL}, which results in
#'   the combination (product) of attractor values for source and destination
#'   locations implemented in the generic \code{Dispersal} base class.
#' @param beta Numeric constant for shaping the effect of distance within
#'   gravity dispersal, i.e. \code{dispersal = f(attractors)/distance^beta}.
#'   Default is 1.
#' @param dispersal_stages Numeric vector of population stages (indices) that
#'   disperse. Default is all stages (when set to \code{NULL}).
#' @param proportion The proportion of the (unstructured or staged) population
#'   that disperses at each time step. This parameter is also used to scale the
#'   the number of dispersal destinations selected when the population is
#'   presence-only and the number of dispersal \code{events} is not defined.
#'   Default is \code{NULL} (producing no dispersal unless the population is
#'   presence-only and \code{events} is defined).
#' @param events The mean number of dispersal events generated via a Poisson
#'   distribution for each location at each time step. A dispersal destination
#'   (location) is selected for each dispersal event. Default is \code{NULL}
#'   (resulting in destinations being selected for each individual within
#'   unstructured or staged populations, or stochastic sampling of destinations
#'   for presence-only populations).
#' @param direction_function A function (or kernel) in the form
#'   \code{function(directions)}, that calculates the (relative) probability of
#'   dispersal for each direction (0-360 degrees) specified as an integer
#'   vector. Default is none.
#' @param ... Additional parameters.
#' @return A \code{Gravity} class object (list), containing inherited and
#'   extended functions from the generic \code{Dispersal} class for accessing
#'   attributes (of the function environment) and performing gravity
#'   dispersal:
#'   \describe{
#'     \item{\code{pack(n)}}{Packs a population vector or matrix \code{n} into
#'       a list containing a vector of occupied \code{cells} (indices), the
#'       \code{original} population values at the occupied locations only, the
#'       \code{remaining} occupied values (initially a duplicate
#'       \code{original}), and a vector or matrix for the \code{relocated}
#'       population values at all locations (initially all zero).}
#'     \item{\code{unpack(n)}}{Unpacks a population list by combining the
#'       \code{remaining} and \code{relocated} population values to form a
#'       new post-dispersal population vector or matrix.}
#'     \item{\code{disperse(n)}}{Perform location dispersal on a list \code{n}
#'       of vectors or matrices, representing the occupied \code{cells}
#'       (indices), the \code{original} occupied populations, the
#'       \code{remaining} occupied populations, and the \code{relocated}
#'       populations (at all region cells), and return the transformed list of
#'       vectors or matrices. The separation of original, remaining and
#'       relocated populations enables multiple models for different dispersal
#'       vectors to run in sequence.}
#'   }
#' @include Region.R
#' @include Attractor.R
#' @include Dispersal.R
#' @export
Gravity <- function(region, population_model, attractors,
                    attractor_function = NULL,
                    beta = 1,
                    dispersal_stages = NULL,
                    proportion = NULL,
                    events = NULL,
                    direction_function = NULL, ...) {

  # Check the attractors
  if (!is.list(attractors) ||
      length(attractors) && !all(sapply(1:length(attractors), function(i) {
        inherits(attractors[[i]], "Attractor")
      }))) {
    stop("Attractors must be a list containing zero or more 'Attractor' or ",
         "inherited class objects.", call. = FALSE)
  }

  # Check attractor function and beta parameters
  if (!is.null(attractor_function) && !is.function(attractor_function)) {
    stop("The attractor function should be defined as a function.",
         call. = FALSE)
  }
  if (!is.null(beta) && (!is.numeric(beta) || beta <= 0)) {
    stop("The beta parameter must be numeric and > 0.", call. = FALSE)
  }

  # Run the attractor function to calculate a transformed list of Attractor
  # class objects
  if (!is.null(attractor_function)) {

    # Apply function
    attractors <- attractor_function(attractors)

    # Check the transformed list of attractors again
    if (!is.list(attractors) ||
        length(attractors) && !all(sapply(1:length(attractors), function(i) {
          inherits(attractors[[i]], "Attractor")
        }))) {
      stop("The attractor function should transform a list of 'Attractor' or ",
           "inherited class objects and return another list of class objects.",
           call. = FALSE)
    }
  }

  # Set the distance function as per gravity: f(attractors)/distance^beta
  distance_function <- function(distances) {
    return(1/distances^beta)
  }

  # Build via base class
  self <- Dispersal(region, population_model,
                    dispersal_stages = dispersal_stages,
                    proportion = proportion,
                    events = events,
                    distance_function = distance_function,
                    direction_function = direction_function,
                    attractors = attractors,
                    class = "Gravity", ...)

  return(self)
}
