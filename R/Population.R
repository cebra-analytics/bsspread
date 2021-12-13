#' Population model class builder
#'
#' Builds a class for population representation and growth functionality for
#' the spread simulations.
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
#'   for when \code{type = "presence_only"}.
#' @param class Character class name for inherited classes. Default is
#'   \code{NULL}.
#' @param ... Additional parameters.
#' @return A \code{Population} class object (list) containing functions for
#'   accessing attributes (when present) and simulating growth:
#'   \describe{
#'     \item{\code{get_type()}}{Get the population representation type.}
#'     \item{\code{get_growth()}}{Get the growth rate or age/stage transition
#'       matrix.}
#'     \item{\code{get_stages()}}{Get the number of stages in a age/stage-based
#'       representation.}
#'     \item{\code{get_capacity()}}{Get the carrying capacity as a vector of
#'       values for each location.}
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
                              class = character(), ...) {
  # Validate based on type
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

  # Create a class structure
  self <- structure(list(), class = c(class, "Population"))

  # Get type
  self$get_type <- function() {
    return(type)
  }

  # Get growth when present
  if (!is.null(growth)) {
    self$get_growth <- function() {
      return(growth)
    }
  }

  # Get stages when present
  if (!is.null(stages)) {
    self$get_stages <- function() {
      return(stages)
    }
  }

  # Get carrying capacity when present
  if (!is.null(capacity)) {
    self$get_capacity <- function() {
      return(capacity)
    }
  }

  # Generic grow method (overridden in subclasses)
  self$grow <- function(x) return(x) # no change

  return(self)
}
