#' Population model class builder
#'
#' Builds a class for population representation and growth functionality for
#' spread simulations.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for spread simulations.
#' @param type One of \code{"presence_only"} (default), \code{"unstructured"},
#'   or \code{"stage_structured"} to indicate how populations are represented.
#' @param growth Numeric growth rate or lambda (e.g. 1.2 for 20% growth per
#'   time step) when \code{type = "unstructured"}, or age/stage-based
#'   transition rates (matrix) for each time step when
#'   \code{type = "stage_structured"}. Default is \code{NULL} for when
#'   \code{type = "presence_only"}.
#' @param capacity A (static) vector or matrix (containing temporal columns) of
#'   carrying capacity values of the invasive species at each location (row)
#'   specified by the \code{region}, or per unit area defined by
#'   \code{capacity_area} when the \code{region} is spatially implicit. Default
#'   is \code{NULL} for when \code{type = "presence_only"} or growth is not
#'   capacity-limited. The number of columns for temporal capacity should
#'   coincide with the number of simulation time steps, or a cyclic pattern
#'   (e.g. 12 columns for seasonal variation with monthly time steps).
#' @param capacity_area Carrying capacity area (in m^2) specified for
#'   capacity-limited growth in a spatially implicit (single patch)
#'   \code{region}. For example, use a value of \code{1e+06} when
#'   \code{capacity} is specified in km^2. Default is \code{NULL}, although a
#'   value is required when \code{capacity} is specified for a spatially
#'   implicit region.
#' @param establish_pr An optional (static) vector or matrix (containing
#'   temporal columns) of probability values (0-1) to represent the likelihood
#'   of establishment at each location (row) specified by the \code{region}.
#'   This may be used to avoid transient/unsuccessful incursion or migration
#'   arrivals from being presented in the simulation results, and/or from
#'   subsequently contributing to spread in presence-only models. Default is
#'   \code{NULL}. The number of columns for temporal capacity should
#'   coincide with the number of simulation time steps, or a cyclic pattern
#'   (e.g. 12 columns for seasonal variation with monthly time steps).
#' @param incursion_mean Numeric mean population size for unstructured and
#'   stage structured populations applied at incursion locations. The
#'   population size is sampled from the Poisson distribution for each
#'   incursion location.
#' @param class Character class name for inherited classes. Default is
#'   \code{NULL}.
#' @param ... Additional parameters.
#' @return A \code{Population} class object (list) containing functions for
#'   accessing attributes, conforming data to the appropriate population
#'   representation, and simulating growth:
#'   \describe{
#'     \item{\code{get_region()}}{Get the region object.}
#'     \item{\code{get_type()}}{Get the population representation type.}
#'     \item{\code{get_growth()}}{Get the growth rate or age/stage transition
#'       matrix.}
#'     \item{\code{get_stages()}}{Get the number of stages in a age/stage-based
#'       representation.}
#'     \item{\code{get_capacity(cells = NULL, tm = NULL)}}{Get the carrying
#'       capacity as a vector of values for each region location or optionally
#'       specified region locations \code{cells} (indices) at (optional)
#'       simulation time step \code{tm} (for temporally defined capacity).}
#'     \item{\code{get_establish_pr(cells = NULL, tm = NULL)}}{Get the
#'       establishment probability as a vector of values for each region
#'       location or optionally specified region locations \code{cells}
#'       (indices) at (optional) simulation time step \code{tm} (for
#'       temporally defined establishment probability).}
#'     \item{\code{set_incursion_mean(m)}}{Set the incursion mean.}
#'     \item{\code{make(initial, current, incursion, tm)}}{Make a population
#'       vector or array (rows:stages x columns:locations) via the defined
#'       population representation using vectors or arrays of the
#'       \code{initial} or \code{current} and \code{incursion} population at
#'       each region location at simulation time step \code{tm}.}
#'     \item{\code{grow(x, tm)}}{Performs growth or age/stage-based transitions
#'       on population \code{x} vector/array (stages by locations) at
#'       simulation time step \code{tm}, and returns the transformed
#'       vector/array.}
#'   }
#' @references
#'   Beverton, R. J. H., & Holt, S. J. (1957). On the dynamics of exploited
#'   fish populations. \emph{Fisheries Investigations}, 19, 1-533.
#'
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
#'
#'   Jongejans, E., Skarpaas, O., & Shea, K. (2008). Dispersal, demography and
#'   spatial population models for conservation and control management.
#'   \emph{Perspectives In Plant Ecology Evolution And Systematics}, 9(3–4),
#'   153–170. \doi{10.1016/j.ppees.2007.09.005}
#'
#'   Ricker, W. E. (1958). Handbook of computations for biological statistics
#'   of fish populations. \emph{Bulletin (Fisheries Research Board of Canada)},
#'   119. \url{https://waves-vagues.dfo-mpo.gc.ca/Library/1485.pdf}
#' @include Region.R
#' @export
Population <- function(region,
                       type = c("presence_only", "unstructured",
                                "stage_structured"),
                       growth = NULL,
                       capacity = NULL,
                       capacity_area = NULL,
                       establish_pr = NULL,
                       incursion_mean = NULL,
                       class = character(), ...) {
  UseMethod("Population")
}

#' @name Population
#' @export
Population.Region <- function(region,
                              type = c("presence_only", "unstructured",
                                       "stage_structured"),
                              growth = NULL,
                              capacity = NULL,
                              capacity_area = NULL,
                              establish_pr = NULL,
                              incursion_mean = NULL,
                              class = character(), ...) {

  # Validate based on type and set (number of) stages
  type <- match.arg(type)
  if (type == "presence_only") {
    stages = NULL
  } else if (type == "unstructured") {
    if (length(growth) > 1) {
      stop("Population growth should be a single value (e.g. 1.2).",
           call. = FALSE)
    }
    stages = NULL
  } else if (type == "stage_structured") {
    if (is.numeric(growth)) growth <- as.matrix(growth)
    if (!length(growth) >= 4 || nrow(growth) != ncol(growth)) {
      stop("Population growth should be a square matrix (at least 2 x 2).",
           call. = FALSE)
    }
    stages = nrow(growth)
  } else {
    stop("Population type should be 'presence_only', 'unstructured', or ",
         "'stage_structured'.", call. = FALSE)
  }

  # Validate capacity and area via region
  if (!is.null(capacity) && (!is.numeric(capacity) ||
                             length(capacity) != region$get_locations())) {
    stop("Population capacity should be a vector with a value for each ",
         "region location.", call. = FALSE)
  }
  if (!is.null(capacity_area) && (!is.numeric(capacity_area) ||
                                  capacity_area <= 0)) {
    stop("Population capacity area should be a numeric value > 0.",
         call. = FALSE)
  }

  # Capacity area required when capacity specified and spatially implicit
  if (is.null(capacity_area) && is.numeric(capacity) &&
      region$spatially_implicit()) {
    stop("Population capacity area is required when capacity is specified ",
         "and the region is spatially implicit.", call. = FALSE)
  }
  if (is.numeric(capacity_area)) { # attach area
    attr(capacity, "area") <- capacity_area
  }

  # Validate establishment probability via region
  if (!is.null(establish_pr) &&
      (!is.numeric(establish_pr) ||
       length(establish_pr) != region$get_locations())) {
    stop("Establishment probability should be a vector with a value for each ",
         "region location.", call. = FALSE)
  }

  # Validate incursion mean
  if (!is.null(incursion_mean) && (!is.numeric(incursion_mean) ||
                                   incursion_mean <= 0)) {
    stop("Incursion mean population size should be a numeric value greater ",
         "than zero.", call. = FALSE)
  }

  # Create a class structure
  self <- structure(list(), class = c(class, "Population"))

  # Get the region object get_region
  self$get_region <- function() {
    return(region)
  }

  # Get type
  self$get_type <- function() {
    return(type)
  }

  # Get growth
  self$get_growth <- function() {
    return(growth)
  }

  # Get stages
  self$get_stages <- function() {
    return(stages)
  }

  # Get carrying capacity for specified region (non-NA) cell indices
  # at specified time step when capacity is temporal (matrix)
  self$get_capacity <- function(cells = NULL, tm = NULL) {
    if (is.numeric(capacity) && is.numeric(cells)) { # TODO temporal
      return(capacity[cells])
    } else {
      return(capacity)
    }
  }

  # Get establishment probability for specified region (non-NA) cell indices
  # at specified time step when capacity is temporal (matrix)
  self$get_establish_pr <- function(cells = NULL, tm = NULL) { # TODO temporal
    if (is.numeric(establish_pr) && is.numeric(cells)) {
      return(establish_pr[cells])
    } else {
      return(establish_pr)
    }
  }

  # Set the incursion mean
  self$set_incursion_mean <- function(m) {
    incursion_mean <<- m
  }

  # Generic make method (extended/overridden in subclasses)
  self$make <- function(initial = NULL, current = NULL, incursion = NULL,
                        tm = NULL) { # TODO temporal

    # Initial values only (ignore current and incursion)
    if (!is.null(initial)) {

      # Check consistency with locations and stages
      if (!nrow(as.matrix(initial)) %in% c(1, region$get_locations())) {
        stop("Initial population values must be consistent with the number ",
             "of region locations.", call. = FALSE)
      }
      if (!ncol(as.matrix(initial)) %in% c(1, stages)) {
        stop("Initial population values must be consistent with the number ",
             "of stages.", call. = FALSE)
      }

      # Repeat single (stage) values for each location
      if (nrow(as.matrix(initial)) == 1) {
        if (ncol(as.matrix(initial)) > 1) { # stages
          initial <- apply(as.matrix(initial), 2, rep, region$get_locations())
        } else {
          initial <- rep(initial, region$get_locations())
        }
      }

      return(initial)
    }

    # Initial or ongoing incursions (combined with current in subclasses)
    if (is.logical(incursion)) {

      # Apply stochastic establishment to arrival probability-based incursions
      if (is.numeric(establish_pr) && attr(incursion, "type") == "prob") {

        # Indices of incursion locations
        indices <- which(incursion)

        # Select incursions via binomial sampling
        incursion[indices] <- as.logical(
          stats::rbinom(length(indices), size = 1,
                        prob = establish_pr[indices]))
      }

      # Generate population size values
      if (!is.null(incursion_mean)) {

        # Indices of incursion locations
        indices <- which(incursion)

        if (type != "presence_only") { # to integer
          incursion <- as.integer(incursion)
        }

        if (length(indices)) { # incursions present

          # Sample incursion population sizes via the Poisson distribution
          values <- stats::rpois(length(indices), incursion_mean)

          # Fill values
          incursion[indices] <- values
        }
      }

      # Combine with current if columns match, else error
      if (!is.null(current)) {
        if (ncol(as.matrix(current)) == ncol(as.matrix(incursion))) {
          if (is.logical(incursion) && is.logical(current)) {
            incursion <- incursion | current
          } else { # numeric
            incursion <- incursion + current
          }
        } else {
          stop("Cannot combine incursion with current population array as ",
               "the columns are inconsistent.", call. = FALSE)
        }
      }

      return(incursion)
    }
  }

  # Generic grow method (overridden in inherited classes)
  self$grow <- function(x, tm) return(x) # no change

  return(self)
}
