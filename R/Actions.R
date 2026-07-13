#' Actions class builder
#'
#' Builds a generic class for simulating the application of invasive species
#' management actions, such as detection, control, and removal.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the population spread model simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the population spread model
#'   simulations.
#' @param type One of \code{"detection"} (default), \code{"control"},
#'   or \code{"removal"} to indicate the type of actions applied.
#' @param stages Numeric vector of population stages (indices) to which
#'   actions are applied. Default is all stages (when set to \code{NULL}).
#' @param schedule Vector of discrete simulation time steps (t = 0, 1, 2, ...)
#'   in which to apply actions. Default is all time steps (when set to
#'   \code{NULL}).
#' @param ... Additional parameters.
#' @return A \code{Actions} class object (list) containing a function for
#'   applying simulated actions:
#'   \describe{
#'     \item{\code{get_type()}}{Get the actions type.}
#'     \item{\code{get_id()}}{Get the actions numeric identifier.}
#'     \item{\code{set_id(id)}}{Set the actions numeric identifier. Should be
#'       an integer >= 1).}
#'     \item{\code{get_label(include_id = TRUE)}}{Get the actions label used
#'       in simulation results. Set \code{include_id} to include the action
#'       \code{id} as a label prefix (default is \code{TRUE}).}
#'     \item{\code{get_stages()}}{Get the population stages to which actions
#'       are applied.}
#'     \item{\code{get_schedule()}}{Get the scheduled simulation time steps in
#'       which actions are applied.}
#'     \item{\code{include_cost()}}{Logical indication of a cost parameter
#'       having a value.}
#'     \item{\code{get_cost_label(include_id = TRUE)}}{Get the action cost
#'       label used in simulation results. Set \code{include_id} to include
#'       the action \code{id} as a label prefix (default is \code{TRUE}).}
#'     \item{\code{get_cost_unit()}}{Get the unit of action cost.}
#'     \item{\code{clear_attributes(n)}}{Clear attached attributes associated
#'       with this action from a simulated population vector or matrix
#'       \code{n}, and return \code{n} without the attached attributes.}
#'     \item{\code{apply(n, tm)}}{Apply actions to a simulated population
#'       vector or matrix \code{n}, potentially with attached attributes
#'       relating to previously applied actions, providing the time step
#'       \code{tm} is in the \code{schedule}, and return the resulting
#'       population \code{n} along with attached attributes relating to the
#'       newly applied actions.}
#'   }
#' @export
Actions <- function(region, population_model,
                    type = c("detection", "control", "removal"),
                    stages = NULL,
                    schedule = NULL,
                    class = character(), ...) {
  UseMethod("Actions")
}

#' @name Actions
#' @export
Actions.Region <- function(region, population_model,
                           type = c("detection", "control", "removal"),
                           stages = NULL,
                           schedule = NULL,
                           class = character(), ...) {

  # Check the population model
  if (!is.null(population_model) &&
      !inherits(population_model, "Population")) {
    stop("Population model must be a 'Population' or inherited class object.",
         call. = FALSE)
  } else if (population_model$get_region()$get_locations() !=
             region$get_locations()) {
    stop("Population model must be compatible with the region object.",
         call. = FALSE)
  }

  # Check the applied stages
  population_type <- population_model$get_type()
  if (population_type == "presence_only" ||
      population_type == "unstructured") {
    stages <- 1
  } else if (population_type == "stage_structured") {
    if (is.null(stages)) {
      stages <- 1:population_model$get_stages()
    } else if (!is.numeric(stages) ||
               !all(stages %in% 1:population_model$get_stages())) {
      stop(paste("Stages must be a vector of stage indices consistent with",
                 "the population model."), call. = FALSE)
    }
  }

  # Check the time step schedule for action application
  if (!is.null(schedule) && !is.numeric(schedule)) {
    stop(paste("The schedule for applying actions should be a vector of",
               "numeric simulation time steps."), call. = FALSE)
  }

  type <- match.arg(type)

  # Create a class structure
  self <- structure(list(), class = c(class, "Actions"))

  # Get type
  self$get_type <- function() {
    return(type)
  }

  # Id for tracking multiple actions
  id <- NULL

  # Get the actions id
  self$get_id <- function() {
    return(id)
  }

  # Set the actions id
  self$set_id <- function(id) {
    if (!is.numeric(id) || trunc(id) != id || id < 1) {
      stop("Actions id should be an integer >= 1.", call. = FALSE)
    }
    id <<- id
  }

  # Get results label (overridden in inherited classes)
  self$get_label <- function(include_id = TRUE) {
    if (!is.null(id) && include_id) {
      return(paste0(id, "_", "action"))
    } else {
      return("action")
    }
  }

  # Get stages
  self$get_stages <- function() {
    return(stages)
  }

  # Get schedule
  self$get_schedule <- function() {
    return(schedule)
  }

  # Does cost parameter (named) having a value?
  self$include_cost <- function() {
    # overridden in inherited classes
    return(FALSE)
  }

  # Get cost results label (overridden in inherited classes)
  self$get_cost_label <- function(include_id = TRUE) {
    if (!is.null(id) && include_id) {
      return(paste0(id, "_", "action_cost"))
    } else {
      return("action_cost")
    }
  }

  # Get the unit of action cost
  self$get_cost_unit <- function() { # overridden in inherited classes
    return(NULL)
  }

  # Clear attached attributes
  self$clear_attributes <- function(n) { # overridden in inherited classes
    return(n)
  }

  # Generic apply method (overridden in inherited classes)
  self$apply <- function(x, tm) return(x) # no change

  return(self)
}
