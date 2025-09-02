#' Diffusion dispersal model class builder
#'
#' Builds a class for representing diffusive dispersal for spread simulations
#' with functionality inherited from the generic dispersal model.
#' Diffusion simulates the local spread of populations into neighbouring
#' locations (cells) a speed specified by the \code{diffusion_rate} parameter.
#' Diffusion may be simulated for presence-only, unstructured or stage-based
#' \code{populations}. For unstructured or stage-based populations, a specified
#' \code{proportion} of the population at each occupied (cell) location is
#' selected (sampled) for diffusive dispersal at each simulation time step.
#' Presence-only populations may diffuse into a specified \code{proportion} of
#' neighbouring locations that are reachable at the specified diffusion rate
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
#' @param allow_contraction Optional logical indicator to allow the contraction
#'   of the occupied area when spatially implicit reaction-diffusion models
#'   undergo population decline. Default is \code{TRUE} to allow area
#'   contraction. Set to \code{FALSE} when declining threats are likely to
#'   retain an occupied area, albeit sparsely.
#' @param proportion The proportion of the (unstructured or staged) population
#'   that disperses from each occupied location at each time step, or the
#'   proportion of local presence-only destinations selected for diffusive
#'   dispersal. It may be specified as a single numeric value or, when
#'   applicable, with spatial and/or temporal variation via a matrix of spatial
#'   (rows) and/or temporal (columns). Spatial values should be specified via a
#'   row for each location, else a single row may specify temporal variation
#'   only. Likewise, a single column may specify spatial variation only. The
#'   number of columns for temporal variation should either coincide with the
#'   number of simulation time steps, or be a cyclic pattern (e.g. 12 columns
#'   for seasonal variation with monthly time steps).. Default is \code{NULL}
#'   (producing no diffusive dispersal).
#' @param density_dependent Logical to indicate that dispersal is density
#'   dependent, whereby the proportion dispersing is scaled by the
#'   (unstructured or staged) population density (number/capacity) at each
#'   occupied location at each simulation time step. Default is \code{FALSE}
#'   for no density dependence.
#' @param direction_function A function (or kernel) in the form
#'   \code{function(directions)}, that calculates the (relative) probability of
#'   dispersal for each direction (0-360 degrees) specified as an integer
#'   vector. Default is none.
#' @param attractors List containing \code{Attractor} (or inherited) class
#'   objects for spatially weighted dispersal to destination locations. Default
#'   is empty.
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
#'     \item{\code{disperse(n, tm)}}{Perform location dispersal at simulation
#'       time step \code{tm} on a list \code{n} of vectors or matrices,
#'       representing the occupied \code{cells} (indices), the \code{original}
#'       occupied populations, the \code{remaining} occupied populations, and
#'       the \code{relocated} populations (at all region cells), and return
#'       the transformed list of vectors or matrices. The separation of
#'       original, remaining and relocated populations enables multiple models
#'       for different dispersal vectors to run in sequence. Spatially implicit
#'       diffusion attaches a \code{diffusion_radius} attribute to \code{n} at
#'       each time step.}
#'     \item{\code{get_diffusion_rate()}}{Get the defined or calculated
#'       (asymptotic) speed of diffusion.}
#'   }
#' @references
#'   Andow, D. A., Kareiva, P. M., Levin, S. A., & Okubo, A. (1990). Spread of
#'   invading organisms. \emph{Landscape Ecology}, 4(2–3), 177.
#'   \doi{10.1007/bf00132860}
#'
#'   Bradhurst, R., Spring, D., Stanaway, M., Milner, J., & Kompas, T. (2021).
#'   A generalised and scalable framework for modelling incursions,
#'   surveillance and control of plant and environmental pests.
#'   \emph{Environmental Modelling & Software}, 139, N.PAG.
#'   \doi{10.1016/j.envsoft.2021.105004}
#'
#'   Fisher, R. A. (1937). The wave of advance of advantageous genes.
#'   \emph{Ann. Eugenics} 7, 355–369. \doi{10.1111/j.1469-1809.1937.tb02153.x}
#'
#'   Jongejans, E., Skarpaas, O., & Shea, K. (2008). Dispersal, demography and
#'   spatial population models for conservation and control management.
#'   \emph{Perspectives In Plant Ecology Evolution And Systematics}, 9(3–4),
#'   153–170. \doi{10.1016/j.ppees.2007.09.005}
#'
#'   Okubo, A., & Kareiva, P. (2001). Some Examples of Animal Diffusion. In
#'   A. Okubo & S. A. Levin (Eds.) \emph{Diffusion and Ecological Problems:}
#'   \emph{Modern Perspectives} (pp. 238-267). Springer New York.
#'   \doi{10.1007/978-1-4757-4978-6}
#'
#'   Robinet, C., Kehlenbeck, H., Kriticos, D. J., Baker, R. H. A.,
#'   Battisti, A., Brunel, S., Dupin, M., Eyre, D., Faccoli, M., Ilieva, Z.,
#'   Kenis, M., Knight, J., Reynaud, P., Yart, A., & van der Werf, W. (2012).
#'   A Suite of Models to Support the Quantitative Assessment of Spread in Pest
#'   Risk Analysis. \emph{PLoS ONE}, 7(10), 1–18.
#'   \doi{10.1371/journal.pone.0043366}
#'
#'   Skellam, J. G. (1951). Random Dispersal in Theoretical Populations.
#'   \emph{Biometrika}, 38(1/2), 196–218. \doi{10.2307/2332328}
#' @include Region.R
#' @include Dispersal.R
#' @export
Diffusion <- function(region, population_model,
                      dispersal_stages = NULL,
                      diffusion_rate = NULL,
                      diffusion_coeff = NULL,
                      allow_contraction = TRUE,
                      proportion = NULL,
                      density_dependent = FALSE,
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

  # Check diffusion coefficient
  if (!is.null(diffusion_coeff) &&
      (!is.numeric(diffusion_coeff) || diffusion_coeff <= 0)) {
    stop("The diffusion coefficient must be numeric and > 0.", call. = FALSE)
  }

  # Check allow contraction indicator
  if (!is.logical(allow_contraction)) {
    stop("The allow contraction indicator must be logical TRUE or FALSE.",
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
    if (!is.null(diffusion_rate) && diffusion_rate > 0) {
      max_distance <- diffusion_rate
    } else {
      max_distance <- NULL
    }
    combined_function <- NULL
  }

  # Build via base class
  self <- Dispersal(region, population_model,
                    dispersal_stages = dispersal_stages,
                    proportion = proportion,
                    events = NULL,
                    density_dependent = density_dependent,
                    distance_function = NULL,
                    direction_function = direction_function,
                    combined_function = combined_function,
                    distance_adjust = FALSE,
                    attractors = attractors,
                    permeability = permeability,
                    max_distance = max_distance,
                    class = "Diffusion", ...)

  # Resolve dispersal stages
  if (population_model$get_type() %in% c("presence_only", "unstructured")) {
    dispersal_stages <- 1
  } else if (population_model$get_type() == "stage_structured" &&
             is.null(dispersal_stages)) {
    dispersal_stages <- 1:population_model$get_stages()
  }

  # Spatially implicit (single patch)
  if (region$spatially_implicit()) {

    # Radial diffusion or reaction-diffusion model
    if (population_model$get_type() == "presence_only") { # radial diffusion

      # Override disperse function
      self$disperse <- function(n, tm) {

        # Calculate diffusion radius
        if (length(n$original) > 0 && n$original) {
          radius_increase <- diffusion_rate
        } else {
          radius_increase <- 0
        }
        if (is.numeric(attr(n$relocated, "diffusion_radius"))) {
          diffusion_radius <-
            attr(n$relocated, "diffusion_radius") + radius_increase
        } else {
          diffusion_radius <- radius_increase
        }

        # Limit via maximum area when applicable
        if (is.numeric(region$get_max_implicit_area())) {
          max_radius <- sqrt(region$get_max_implicit_area()/pi)
          diffusion_radius <- min(diffusion_radius, max_radius)
        }

        # Attach attribute for diffusion radius
        attr(n$relocated, "diffusion_radius") <- diffusion_radius

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
      self$disperse <- function(n, tm) {

        # Extract initial population size
        if (is.numeric(attr(n$relocated, "initial_n"))) {
          initial_n <- attr(n$relocated, "initial_n")
        } else {
          stop(paste("The initial population needs to be set as an attribute",
                     "for reaction-diffusion calculations."), call. = FALSE)
        }

        # Calculate radius via reaction-diffusion (Okubo & Kareiva, 2001)
        # m' = initial_n*exp(intrinsic_r*tm - Radius^2/(4*diffusion_coeff*tm))
        # where n_t = initial_n*exp(intrinsic_r*tm) for exponential growth
        m_dash <- sum(initial_n) # since initial radius is zero
        diffusion_radius <-
          sqrt(4*diffusion_coeff*tm*log(max(sum(n$original)/m_dash, 1)))

        # Limit via maximum area when applicable
        if (is.numeric(region$get_max_implicit_area())) {
          max_radius <- sqrt(region$get_max_implicit_area()/pi)
          diffusion_radius <- min(diffusion_radius, max_radius)
        }

        # Attach attribute for diffusion radius when grown
        if (is.null(attr(n$relocated, "diffusion_radius")) ||
            allow_contraction ||
            (is.numeric(attr(n$relocated, "diffusion_radius")) &&
             diffusion_radius > attr(n$relocated, "diffusion_radius"))) {
          attr(n$relocated, "diffusion_radius") <- diffusion_radius
        }

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
