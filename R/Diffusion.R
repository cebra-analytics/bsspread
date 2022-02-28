#' Diffusion dispersal model class builder
#'
#' Builds a class for representing diffusive dispersal for spread simulations
#' with functionality inherited from the generic dispersal model.
#'
#' Diffusion simulates the local spread of populations into neighboring
#' locations (cells) a speed specified by the \code{diffusion_rate} parameter.
#' Diffusion may be simulated for presence-only, unstructured or stage-based
#' \code{populations}. For unstructured or stage-based populations, a specified
#' \code{proportion} of the population at each occupied (cell) location is
#' selected (sampled) for diffusive dispersal at each simulation time step.
#' Presence-only populations may diffuse into a specified \code{proportion} of
#' neighboring locations that are reachable at the specified diffusion rate
#' within each simulation time step. The probability of diffusive dispersal may
#' optionally be adjusted via \code{direction} functions and/or
#' \code{attractor} (layer) values, as well as \code{permeability} or
#' constraint layers, which are used to adjust (effective) distances between
#' cells and/or omit unreachable destinations. An optional establishment
#' likelihood (layer), which is configured via the population model, may be
#' applied to each diffusive dispersal, resulting in potential "deaths" of
#' individuals or unsuccessful presence-only dispersal events. The dispersal
#' functionality utilizes a wrapped population list of separate values for the
#' \code{original}, \code{remaining}, and \code{relocated} populations. This
#' separation enables multiple dispersal models, representing different
#' dispersal vectors, to be run in sequence. For example, local diffusion may
#' be combined with human-mediated and/or wind-based long-distant "jump"
#' dispersal.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for spread simulations. The region object
#'   contains functionality for calculating path distances and directions,
#'   permeability graphs, and structures to facilitate two-tier dispersal.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the spread simulations.
#' @param dispersal_stages Numeric vector of population stages (indices) that
#'   disperse. Default is all stages (when set to \code{NULL}).
#' @param diffusion_rate The speed of diffusion (in m per time step). Default
#'   is \code{NULL} (resulting in spread to adjacent locations or cells at each
#'   time step when other sampling probabilities result in selection).
#' @param proportion The proportion of the (unstructured or staged) population
#'   that disperses at each time step, or the proportion of local presence-only
#'   destinations selected for diffusive dispersal. Default is \code{NULL}
#'   (producing no diffusive dispersal).
#' @param direction_function A function (or kernel) in the form
#'   \code{function(directions)}, that calculates the (relative) probability of
#'   dispersal for each direction (0-360 degrees) specified as an integer
#'   vector. Default is none.
#' @param attractors List containing \code{Attractor} (or inherited) class
#'   objects for spatially weighted dispersal, and/or containing a named
#'   numeric multiplier \code{source_density = <numeric>}, which (when present)
#'   is used to dynamically create an "attractor" based on the population
#'   density (number/capacity) at each dispersal location source (origin) at
#'   each simulation time step. Default is empty.
#' @param permeability A \code{Permeability} class (or inherited) class object
#'   for representing spatial permeability or constraints. Default is none.
#' @param ... Additional parameters.
#' @return A \code{Diffusion} class object (list), containing inherited and
#'   extended functions from the generic \code{Dispersal} class for accessing
#'   attributes (of the function environment) and performing diffusive
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
#' @include Dispersal.R
#' @export
Diffusion <- function(region, population_model,
                      dispersal_stages = NULL,
                      diffusion_rate = NULL,
                      proportion = NULL,
                      direction_function = NULL,
                      attractors = list(),
                      permeability = NULL, ...) {

  # Check diffusion rate and set its default (when empty)
  if (!is.null(diffusion_rate) &&
      (!is.numeric(diffusion_rate) || diffusion_rate <= 0)) {
    stop("The diffusion rate must be numeric and > 0.", call. = FALSE)
  }
  if (is.null(diffusion_rate)) { # include adjacent cells
    diffusion_rate <- region$get_res()*sqrt(2)
  }

  # Configure maximum distance and define distance function for diffusion
  diag_dist <- region$get_res()*sqrt(2)/2
  max_distance <- max(diffusion_rate, diag_dist) + diag_dist
  distance_function <- function(distances) {
    return(pmin(diffusion_rate/distances, 1))
  }

  # Build via base class
  self <- Dispersal(region, population_model,
                    dispersal_stages = dispersal_stages,
                    proportion = proportion,
                    events = NULL,
                    distance_function = distance_function,
                    distance_adjust = FALSE,
                    direction_function = direction_function,
                    attractors = attractors,
                    permeability = permeability,
                    max_distance = max_distance,
                    class = "Diffusion", ...)

  return(self)
}
