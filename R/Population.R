#' Population model class builder
#'
#' Builds a class for population representation and growth functionality for
#' spread simulations.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for spread simulations.
#' @param type One of \code{"presence_only"} (default), \code{"unstructured"},
#'   or \code{"stage_structured"} to indicate how populations are represented.
#' @param growth Numeric intrinsic growth rate or lambda (e.g. 1.2 for 20%
#'   growth per time step) when \code{type = "unstructured"}, or
#'   age/stage-based transition rates (matrix) for each time step when
#'   \code{type = "stage_structured"}. Default is \code{NULL} for when
#'   \code{type = "presence_only"}.
#' @param capacity A vector of carrying capacity values of the invasive species
#'   at each location specified by the \code{region}. Default is \code{NULL}
#'   for when \code{type = "presence_only"} or growth is not capacity-limited.
#' @param establish_pr An optional vector of probability values (0-1) to
#'   represent the likelihood of establishment at each location specified by
#'   the \code{region}. This may be used to avoid transient/unsuccessful
#'   incursion or migration arrivals from being presented in the simulation
#'   results, and/or from subsequently contributing to spread in presence-only
#'   models. Default is \code{NULL}.
#' @param incursion_values Defines how incursion locations are populated.
#'   Either via a fixed single value, vector, or matrix of values for each
#'   location (vector or matrix rows) and/or stage (matrix columns), or via a
#'   list specifying the range of values via \code{min} and \code{max} (single,
#'   vector or matrix) for uniform random sampling of values.
#' @param class Character class name for inherited classes. Default is
#'   \code{NULL}.
#' @param ... Additional parameters.
#' @return A \code{Population} class object (list) containing functions for
#'   accessing attributes, conforming data to the appropriate population
#'   representation, and simulating growth:
#'   \describe{
#'     \item{\code{get_type()}}{Get the population representation type.}
#'     \item{\code{get_growth()}}{Get the growth rate or age/stage transition
#'       matrix.}
#'     \item{\code{get_stages()}}{Get the number of stages in a age/stage-based
#'       representation.}
#'     \item{\code{get_capacity()}}{Get the carrying capacity as a vector of
#'       values for each location.}
#'     \item{\code{get_establish_pr()}}{Get the establishment probability as a
#'       vector of values for each location.}
#'     \item{\code{make(initial, current, incursion)}}{Make a population vector
#'       or array (rows:stages x columns:locations) via the defined population
#'       representation using vectors or arrays of the \code{initial} or
#'       \code{current} and \code{incursion} population at each region
#'       location.}
#'     \item{\code{grow(x)}}{Performs growth or age/stage-based transitions on
#'       population \code{x} vector/array (stages by locations), and returns
#'       the transformed vector/array.}
#'   }
#' @include Region.R
#' @export
Population <- function(region,
                       type = c("presence_only", "unstructured",
                                "stage_structured"),
                       growth = NULL,
                       capacity = NULL,
                       establish_pr = NULL,
                       incursion_values = NULL,
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
                              establish_pr = NULL,
                              incursion_values = NULL,
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

  # Validate capacity via region
  if (!is.null(capacity) && (!is.numeric(capacity) ||
                             length(capacity) != region$get_locations())) {
    stop("Population capacity should be a vector with a value for each ",
         "region location.", call. = FALSE)
  }

  # Validate establishment probability via region
  if (!is.null(establish_pr) &&
      (!is.numeric(establish_pr) ||
       length(establish_pr) != region$get_locations())) {
    stop("Establishment probability should be a vector with a value for each ",
         "region location.", call. = FALSE)
  }

  # Validate incursion values
  if (!is.null(incursion_values)) {

    # Check list
    if (is.list(incursion_values)) {
      if (length(incursion_values) != 2) {
        stop("A range of incursion values defined by a list should have two ",
             "entries named: 'min' and 'max'.", call. = FALSE)
      }
      if (!all(c("min", "max") %in% names(incursion_values))) {
        names(incursion_values) <- c("min", "max")
      }
    }

    # Check consistency with locations and stages
    for (values in
         if (is.list(incursion_values)) incursion_values
         else list(incursion_values)) {
      if (!nrow(as.matrix(values)) %in% c(1, region$get_locations())) {
        stop("Incursion values must be consistent with the number of region ",
             "locations.", call. = FALSE)
      }
      if (!ncol(as.matrix(values)) %in% c(1, stages)) {
        stop("Incursion values must be consistent with the number of stages.",
             call. = FALSE)
      }
    }
  }

  # Create a class structure
  self <- structure(list(), class = c(class, "Population"))

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

  # Get carrying capacity
  self$get_capacity <- function() {
    return(capacity)
  }

  # Get establishment probability
  self$get_establish_pr <- function() {
    return(establish_pr)
  }

  # Generic make method (extended/overridden in subclasses)
  self$make <- function(initial = NULL, current = NULL, incursion = NULL) {

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

      # Apply stochastic establishment to incursions
      if (is.numeric(establish_pr)) {

        # Indices of incursion locations
        indices <- which(incursion)

        # Select incursions via binomial sampling
        incursion[indices] <- as.logical(
          stats::rbinom(length(indices), size = 1,
                        prob = establish_pr[indices]))
      }

      if (!is.null(incursion_values)) {

        # Indices of incursion locations
        indices <- which(incursion)

        # Sample across range
        if (is.list(incursion_values)) {

          # Columns to output
          cols <- max(unlist(lapply(incursion_values,
                                    function(v) ncol(as.matrix(v)))))

          if (length(indices)) { # incursions present

            # Shape values for incursions
            values <- list()
            for (key in c("min", "max")) {
              values[[key]] <- as.matrix(incursion_values[[key]])
              if (nrow(values[[key]]) == 1) {
                values[[key]] <- apply(values[[key]], 2, rep, length(indices))
              } else { # a row per location
                values[[key]] <- values[[key]][indices, , drop = FALSE]
              }
            }

            # Sample values
            values <- stats::runif(length(indices)*cols, min = values$min,
                                   max = values$max)
            if (all(unlist(lapply(incursion_values, is.integer))) ||
                type == "stage_structured") {
              values <- as.integer(round(values))
            }
          }

        } else { # fixed incursion values

          # Columns to output
          cols <- ncol(as.matrix(incursion_values))

          if (length(indices)) { # incursions present

            # Shape values for incursions
            values <- as.matrix(incursion_values)
            if (nrow(values) == 1) {
              values <- apply(values, 2, rep, length(indices))
            } else { # a row per location
              values <- values[indices, , drop = FALSE]
            }
          }
        }

        # Shape incursions
        if (cols > 1) {
          incursion <- array(FALSE, c(length(incursion), cols))
        } else {
          incursion <- vector("logical", length(incursion))
        }

        if (length(indices)) { # incursions present

          # Fill values
          if (cols > 1) {
            incursion[indices,] <- values
          } else {
            incursion[indices] <- values
          }
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
  self$grow <- function(x) return(x) # no change

  return(self)
}
