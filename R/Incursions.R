#' Incursions class builder
#'
#' Builds a class to generate invasive species incursions, either by selecting
#' a single incursion location via stochastic sampling based on weightings for
#' each available spatial location, or by sampling potentially multiple
#' locations via incursion or arrival probabilities defined for each location,
#' which optionally may continue to generate incursions at specified time
#' intervals during a spread simulation.
#'
#' @param x A \code{raster::RasterLayer}, \code{terra::SpatRaster}, or numeric
#'   vector defining the incursion weightings or arrival probabilities at each
#'   spatial location.
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for spread simulations.
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
#'   accessing attributes and generating sampled incursions:
#'   \describe{
#'     \item{\code{get_type()}}{Get the incursion type.}
#'     \item{\code{get_continued()}}{Get the continued indicator when type is
#'       "prob".}
#'     \item{\code{get_time_steps()}}{Get the continued incursion time steps
#'       when applicable.}
#'     \item{\code{generate()}}{Generates incursion locations via sampling.}
#'   }
#' @references
#'   Bradhurst, R., Spring, D., Stanaway, M., Milner, J., & Kompas, T. (2021).
#'   A generalised and scalable framework for modelling incursions,
#'   surveillance and control of plant and environmental pests.
#'   \emph{Environmental Modelling & Software}, 139, N.PAG.
#'   \doi{10.1016/j.envsoft.2021.105004}
#'
#'   García Adeva, J. J., Botha, J. H., & Reynolds, M. (2012). A simulation
#'   modelling approach to forecast establishment and spread of Bactrocera
#'   fruit flies. \emph{Ecological Modelling}, 227, 93–108.
#'   \doi{10.1016/j.ecolmodel.2011.11.026}
#' @export
Incursions <- function(x,
                       region = NULL,
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
                                  region = NULL, ...) {

  # Resolve values and call default (vector) version
  if (!is.null(region)) {

    # Check region
    if (!inherits(region, "Region")) {
      stop("Region model must be a 'Region' or inherited class object.",
           call. = FALSE)
    }
    if (!region$is_compatible(x)) {
      stop("The spatial object x should be compatible with that defining the ",
           "region.", call. = FALSE)
    }

    # Extract values from locations defined by region
    Incursions(as.matrix(x[region$get_indices()]), region = region, ...)

  } else { # Use all values
    Incursions(as.matrix(x[]), region = region, ...)
  }
}

#' @name Incursions
#' @export
Incursions.default <- function(x,
                               region = NULL,
                               type = c("weight", "prob"),
                               continued = FALSE,
                               time_steps = 1, ...) {

  # Check region and x
  if (!is.null(region)) {
    if (!inherits(region, "Region")) {
      stop("Region model must be a 'Region' or inherited class object.",
           call. = FALSE)
    }
    if (length(x) != region$get_locations()) {
      stop("Vector x length must be equal to the number of region locations.",
           call. = FALSE)
    }
  }

  # Create a class structure
  self <- structure(list(), class = "Incursions")

  # Get type
  type <- match.arg(type)
  self$get_type <- function() {
    return(type)
  }

  # Get continued
  self$get_continued <- function() {
    return(continued)
  }

  # Get time steps when continued is TRUE
  if (continued) {
    self$get_time_steps <- function() {
      return(time_steps)
    }
  }

  # Generate incursion locations via sampling
  self$generate <- function() {
    incursions <- vector("logical", length(unlist(x)))
    if (type == "weight") { # sample single location
      incursions[sample(1:length(x), 1, prob = unlist(x))] <- TRUE
    } else if (type == "prob") { # binomial sampling at each location
      incursions <- as.logical(stats::rbinom(length(x), 1, prob = unlist(x)))
    }
    attr(incursions, "type") <- type
    return(incursions)
  }

  return(self)
}
