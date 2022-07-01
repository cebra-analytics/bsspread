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
#' dispersal. Spatially implicit diffusion is also provided via simple radial
#' diffusion for a presence-only single-patch population, or via
#' reaction-diffusion for unstructured or stage-based single patch populations.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for spread simulations. The region object
#'   contains functionality for calculating path distances and directions,
#'   permeability graphs, and structures to facilitate two-tier dispersal.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the spread simulations.
#' @param dispersal_stages Numeric vector of population stages (indices) that
#'   disperse. Default is all stages (when set to \code{NULL}).
#' @param diffusion_rate The (asymptotic) speed of diffusion (in m per time
#'   step). Default is \code{NULL} (resulting in no diffusion).
#' @param diffusion_coeff Optional coefficient of diffusion (in m^2 per time
#'   step) for spatially implicit reaction-diffusion. Only used when
#'   unstructured or stage-structured populations are defined within a single
#'   patch region, and the asymptotic \code{diffusion_rate} is unknown (thus
#'   not defined). Default is \code{NULL}.
#' @param diffusion_threshold Optional diffusion boundary threshold for
#'   spatially implicit reaction-diffusion. Defined as a fraction (> 0 and
#'   <= 1) of the initial population outside the diffusion boundary. Only used
#'   when unstructured or stage-structured populations are defined within a
#'   single patch region. Default is 1, resulting in diffusion limited to the
#'   asymptotic diffusion rate, defined via the \code{diffusion_rate}
#'   parameter, or calculated via the \code{diffusion_coeff} parameter.
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
#'       vectors to run in sequence. Spatially implicit diffusion attaches a
#'       \code{diffusion_radius} attribute to \code{n} at each time step.}
#'     \item{\code{get_diffusion_rate()}}{Get the defined or calculated
#'       (asymptotic) speed of diffusion.}
#'   }
#' @include Region.R
#' @include Dispersal.R
#' @export
Diffusion <- function(region, population_model,
                      dispersal_stages = NULL,
                      diffusion_rate = NULL,
                      diffusion_coeff = NULL,
                      diffusion_threshold = 1,
                      proportion = NULL,
                      direction_function = NULL,
                      attractors = list(),
                      permeability = NULL, ...) {

  # Check diffusion rate and set its default (when empty)
  if (!is.null(diffusion_rate) &&
      (!is.numeric(diffusion_rate) || diffusion_rate < 0)) {
    stop("The diffusion rate must be numeric and >= 0.", call. = FALSE)
  }
  if (is.null(diffusion_rate) && is.null(diffusion_coeff)) {
    diffusion_rate <- 0
  }

  # Check diffusion coefficient and threshold
  if (!is.null(diffusion_coeff) &&
      (!is.numeric(diffusion_coeff) || diffusion_coeff <= 0)) {
    stop("The diffusion coefficient must be numeric and > 0.", call. = FALSE)
  }
  if (!is.null(diffusion_threshold) &&
      (!is.numeric(diffusion_threshold) || diffusion_threshold <= 0 ||
       diffusion_threshold > 1)) {
    stop("The diffusion threshold must be numeric and > 0 and <= 1.",
         call. = FALSE)
  }

  # Configure maximum distance and define combined function for diffusion
  if (region$get_type() == "grid") {

    # Select in-range inner cells plus a buffer for stochastic outer cells
    max_distance <- max(region$get_res(), diffusion_rate) + region$get_res()

    # Combined function assigns probabilities to cells
    combined_function <- function(distances, directions) {

      # Set diffusion probability to one when within (minimum) range
      region_res <- region$get_res()
      diff_pr <- +(distances[[1]] <= max(region_res, diffusion_rate))

      # Calculate probability of outer cells within one cell of the range
      outer_idx <- which(diff_pr == 0)
      width <- pmax(abs(cos(directions[outer_idx]*pi/180)),
                    abs(sin(directions[outer_idx]*pi/180)))
      diff_pr[outer_idx] <- (pmax(0, max(region_res, diffusion_rate) +
                                    region_res*width -
                                    distances[[1]][outer_idx])/
                               distances[[1]][outer_idx])

      # Adjust probability when permeable distance is out of range
      if (length(distances) >= 2) {
        perm_outer_idx <-
          which(distances[[2]] > max(region_res, diffusion_rate))
        diff_pr[perm_outer_idx] <-
          (diff_pr*distances[[1]]/distances[[2]])[perm_outer_idx]
      }

      # Scale probability when rate is smaller than cell size
      if (diffusion_rate < region_res) { # TODO: investigate further ####
        diff_pr <- (diff_pr*(diffusion_rate/region_res)^2)
      }
      return(diff_pr)
    }

  } else { # patch

    # Select in-range patches only
    max_distance <- diffusion_rate
    combined_function <- NULL
  }

  # Build via base class
  self <- Dispersal(region, population_model,
                    dispersal_stages = dispersal_stages,
                    proportion = proportion,
                    events = NULL,
                    distance_function = NULL,
                    direction_function = direction_function,
                    combined_function = combined_function,
                    distance_adjust = FALSE,
                    attractors = attractors,
                    permeability = permeability,
                    max_distance = max_distance,
                    class = "Diffusion", ...)

  # Spatially implicit (single patch)
  if (region$spatially_implicit()) {

    # Radial diffusion or reaction-diffusion model
    if (population_model$get_type() == "presence_only") { # radial diffusion

      # Override disperse function
      self$disperse <- function(n) {

        # Attach attribute for diffusion radius
        if (is.numeric(attr(n$relocated, "diffusion_radius"))) {
          attr(n$relocated, "diffusion_radius") <-
            attr(n$relocated, "diffusion_radius") + diffusion_rate
        } else {
          attr(n$relocated, "diffusion_radius") <- diffusion_rate
        }

        return(n)
      }

    } else { # reaction-diffusion

      # Get intrinsic growth rate
      if (population_model$get_type() == "stage_structured") {
        intrinsic_r <- log(population_model$get_growth_r())
      } else { # unstructured
        intrinsic_r <- log(population_model$get_growth())
      }

      # Calculate diffusion rate and coefficient
      if (is.null(diffusion_rate) && is.numeric(diffusion_coeff)) {
        diffusion_rate <- sqrt(4*intrinsic_r*diffusion_coeff)
      }
      if (is.null(diffusion_coeff)) {
        if (intrinsic_r > 0) {
          diffusion_coeff <- diffusion_rate^2/(4*intrinsic_r)
        } else {
          diffusion_coeff <- 0
        }
      }

      # Override disperse function
      self$disperse <- function(n) {

        # Extract initial population size
        if (is.numeric(attr(n$relocated, "initial_n"))) {
          initial_n <- sum(attr(n$relocated, "initial_n"))
        } else {
          stop(paste("The initial population needs to be set as an attribute",
                     "for reaction-diffusion calculations."), call. = FALSE)
        }

        # Extract time step
        if (is.numeric(attr(n$relocated, "tm"))) {
          tm <- attr(n$relocated, "tm")
        } else {
          stop(paste("The current time step needs to be set as an attribute",
                     "for reaction-diffusion calculations."), call. = FALSE)
        }

        # Calculate radius via reaction-diffusion (Okubo & Kareiva, 2001)
        # m' = initial_n*exp(intrinsic_r*tm - Radius^2/(4*diffusion_coeff*tm))
        m_dash <- initial_n*diffusion_threshold
        diffusion_radius <-
          sqrt(4*diffusion_coeff*tm*
                 log(max(initial_n*exp(intrinsic_r*tm)/m_dash, 1)))

        # Attach attribute for diffusion radius
        attr(n$relocated, "diffusion_radius") <- diffusion_radius

        return(n)
      }
    }
  }

  # Get the (asymptotic) speed of diffusion
  self$get_diffusion_rate <- function() {
    return(diffusion_rate)
  }

  return(self)
}
