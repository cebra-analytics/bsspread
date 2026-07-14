#' Controls class builder
#'
#' Builds a class for simulating the application of invasive species management
#' controls, such as combined surveillance and removal ("search & destroy"),
#' and reducing or suppressing growth, reproduction, spread, or establishment.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the population spread model simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the population spread model
#'   simulations.
#' @param control_type One of \code{"search_destroy"} (default),
#'   \code{"growth"} (including reproduction), \code{"spread"}, or
#'   \code{"establishment"} to indicate the type of control applied.
#' @param manage_pr A vector of probability of management control success (or
#'   effectiveness) values at each spatial location for combined surveillance
#'   and removal approaches ("search & destroy", traps, etc.). Required when
#'   \code{"search_destroy"} control type is selected. Default is \code{NULL}
#'   for other control types.
#' @param manage_pr_type Configures how the management control effectiveness
#'   (probability of control success) values are specified. May be specified
#'   at either: the code{"individual"} level (default), whereby effectiveness
#'   values denote the probability of control for each (controllable) invasive
#'   individual; the \code{"population"} level, whereby effectiveness values
#'   denote the probability of controlling all (controllable) individuals
#'   within the (local) invasive population. Only utilised when
#'   \code{"search_destroy"} control type is selected.
#' @param control_mult Control multipliers or rates required for control types
#'   \code{"growth"}, \code{"spread"}, or \code{"establishment"}. May be a
#'   single value (0-1), or vector of values for each location specified by
#'   the \code{region}. Control multipliers are applied to reduce (with values
#'   < 1) or suppress (with values = 0) growth, spread, or establishment at
#'   locations (and optionally surroundings) where the invasive species has
#'   been detected (when specified via a population attribute), as well as at
#'   existing, known, or scheduled treatment locations when
#'   \code{"exist_control"} is specified.
#' @param exist_control A vector to specify the distribution (logical mask) of
#'   existing, known, or scheduled treatment locations for the reduction or
#'   suppression of growth, spread, or establishment (spatially-explicit only).
#'   Default is \code{NULL}. Only utilised when \code{"growth"},
#'   \code{"spread"}, or \code{"establishment"} control type is selected.
#' @param control_cost Numeric vector of distributed control costs (combined
#'   resource and fixed costs) or a single cost value for each location where
#'   control is applied. For spatially-implicit area-based regions, cost
#'   should be specified as either cost for the entire spatially-implicit area
#'   for control type \code{"search_destroy"}, or cost per metres squared for
#'   control type \code{"growth"}, \code{"spread"}, or \code{"establishment"}.
#'   Costs are accumulated for each application of the control at each
#'   (scheduled) simulation time step. The cost unit may be added as an
#'   attribute (\code{attr(control_cost, "unit")}). Default is \code{NULL}
#'   when costs are unavailable.
#' @param radius Optional radius (m) of the applied control for types
#'   \code{"growth"}, \code{"spread"}, or \code{"establishment"} only.
#'   Control is applied to all surrounding locations within the specified
#'   radius of each location where the invasive species has been detected.
#'   Default is \code{NULL}, indicating that control is only applied at
#'   existing, known, or scheduled treatment locations, and/or at detected
#'   locations (when specified via a population attribute).
#' @param stages Numeric vector of population stages (indices) to which
#'   search & destroy or growth controls are applied. Default is all stages
#'   (when set to \code{NULL}).
#' @param apply_to Optional label for growth control to indicate that it should
#'   be only be applied to reproduction or survival rates. If applicable, set
#'   to either \code{"reproduction"} or \code{"survival"}. If not specified,
#'   control is applied to both reproduction and survival. The \code{stages}
#'   and \code{apply_to} parameters may be utilised together, for example, to
#'   control the seasonal survival rates of specified life stages.
#' @param schedule Vector of discrete simulation time steps (t = 0, 1, 2, ...)
#'   in which to apply controls. Default is all time steps (when set to
#'   \code{NULL}).
#' @param ... Additional parameters.
#' @return A \code{Controls} class object (list) containing a function
#'   for accessing attributes and applying simulated controls:
#'   \describe{
#'     \item{\code{get_type()}}{Get the type of action ("control").}
#'     \item{\code{get_id()}}{Get the actions numeric identifier.}
#'     \item{\code{set_id(id)}}{Set the actions numeric identifier. Should be
#'       an integer >= 1).}
#'     \item{\code{get_label(include_id = TRUE)}}{Get the actions label used
#'       in simulation results (e.g."<id>_control_growth" or "control_growth").
#'       Set \code{include_id} to include the action \code{id} as a label
#'       prefix (default is \code{TRUE}).}
#'     \item{\code{get_manage_pr_type()}}{Get the control effectiveness
#'       (probability of management success) type.}
#'     \item{\code{get_stages()}}{Get the population stages to which controls
#'       are applied.}
#'     \item{\code{get_schedule()}}{Get the scheduled simulation time steps in
#'       which controls are applied.}
#'     \item{\code{include_cost()}}{Logical indication of a cost parameter
#'       having a value.}
#'     \item{\code{get_cost_label(include_id = TRUE)}}{Get the control
#'       cost label used in simulation results ("<id>_control_growth_cost" or
#'       "control_growth_cost"). Set \code{include_id} to include the action
#'       \code{id} as a label prefix (default is \code{TRUE}).}
#'     \item{\code{get_cost_unit()}}{Get the unit of control cost.}
#'     \item{\code{get_attributes(n)}}{Get attached attributes associated
#'       with this action from a simulated population vector or matrix
#'       \code{n}, and return a list of the attached attributes.}
#'     \item{\code{clear_attributes(n)}}{Clear attached attributes associated
#'       with this action from a simulated population vector or matrix
#'       \code{n}, and return \code{n} without the attached attributes.}
#'     \item{\code{apply(n, tm)}}{Apply controls to a simulated population
#'       vector or matrix \code{n}, potentially with attached attributes
#'       relating to previously applied actions, providing the time step
#'       \code{tm} is in the \code{schedule}, and return the resulting
#'       population \code{n} along with attached attributes relating to the
#'       newly applied controls.}
#'   }
#' @export
Controls <- function(region, population_model,
                     control_type = c("search_destroy",
                                      "growth",
                                      "spread",
                                      "establishment"),
                     manage_pr = NULL,
                     manage_pr_type = c("individual",
                                        "population"),
                     control_mult = NULL,
                     exist_control = NULL,
                     control_cost = NULL,
                     radius = NULL,
                     stages = NULL,
                     apply_to = NULL,
                     schedule = NULL, ...) {
  UseMethod("Controls")
}

#' @name Controls
#' @export
Controls.Region <- function(region, population_model,
                            control_type = c("search_destroy",
                                             "growth",
                                             "spread",
                                             "establishment"),
                            manage_pr = NULL,
                            manage_pr_type = c("individual",
                                               "population"),
                            control_mult = NULL,
                            exist_control = NULL,
                            control_cost = NULL,
                            radius = NULL,
                            stages = NULL,
                            apply_to = NULL,
                            schedule = NULL, ...) {

  # Build via base class
  self <- Actions(region = region,
                  population_model = population_model,
                  type = "control",
                  stages = stages,
                  schedule = schedule,
                  class = "Controls")

  control_type <- match.arg(control_type)

  # Check and process the management control effectiveness
  if (control_type == "search_destroy") {
    if (is.null(manage_pr)) {
      stop(paste("Management control effectiveness values is required for",
                 "'search and destroy' control."), call. = FALSE)
    }
    if (!is.numeric(manage_pr) ||
        !length(manage_pr) %in% c(1, region$get_locations())) {
      stop(paste("The management control effectiveness parameter must be a",
                 "numeric vector with values for each location."),
           call. = FALSE)
    }
    if (is.numeric(manage_pr) && (any(manage_pr < 0) || any(manage_pr > 1))) {
      stop("Management control effectiveness values should be <= 0 and <= 1.",
           call. = FALSE)
    }
    if (length(manage_pr) == 1) {
      manage_pr <- rep(manage_pr, region$get_locations())
    }
  }

  # Check management control effectiveness type
  manage_pr_type <- match.arg(manage_pr_type)

  # Validate control multiplier and existing control
  if (control_type %in% c("growth", "spread", "establishment")) {
    if (is.null(control_mult)) {
      stop(paste("Control multiplier is required for growth, spread, and",
                 "establishment control."), call. = FALSE)
    } else if (!is.null(control_mult) &&
               (!is.numeric(control_mult) ||
                any(control_mult < 0) || any(control_mult > 1) ||
                !length(control_mult) %in% c(1, region$get_locations()))) {
      stop(paste("Control multiplier should be a vector with a value 0-1",
                 "for each region location."), call. = FALSE)
    }
    if (!is.null(exist_control) &&
        (!(is.numeric(exist_control) || is.logical(exist_control)) ||
         !length(exist_control) == region$get_locations())) {
      stop(paste("Existing control should be a logical or numeric vector with",
                 "a value for each region location."), call. = FALSE)
    } else if (!is.null(exist_control)) {
      exist_control <- as.logical(exist_control)
    }
  }

  # Validate radius, and growth apply to
  if (!is.null(radius) && (!is.numeric(radius) || radius < 0)) {
    stop("The radius (m) parameter must be numeric and >= 0.", call. = FALSE)
  }
  if (!is.null(apply_to)) {
    if (!apply_to %in% c("reproduction", "survival")) {
      stop(paste("Growth control 'apply to' attribute should be",
                 "'reproduction' or 'survival'."), call. = FALSE)
    }
  } else {
    apply_to <- "both"
  }

  # Check and process control cost
  if (!is.null(control_cost)) {
    if (!is.numeric(control_cost) ||
        !length(control_cost) %in% c(1, region$get_locations())) {
      stop(paste("The control cost parameter must be a numeric vector with",
                 "values for each location."), call. = FALSE)
    }
    cost_unit <- attr(control_cost, "unit")
    if (length(control_cost) == 1) {
      if (control_type == "search_destroy") {
        control_cost <- control_cost*(manage_pr > 0)
      } else if (control_type %in% c("growth", "spread", "establishment")) {
        control_cost <-
          control_cost*(control_mult*rep(1, region$get_locations()) < 1)
      }
    }
    attr(control_cost, "unit") <- cost_unit
  }

  # Get results label
  self$get_label <- function(include_id = TRUE) {
    if (!is.null(self$get_id()) && include_id) {
      return(paste0(self$get_id(), "_", "control_", control_type))
    } else {
      return(paste0("control_", control_type))
    }
  }

  # Get the control effectiveness type
  self$get_manage_pr_type <- function() {
    return(manage_pr_type)
  }

  # Does the control cost parameter have a value?
  self$include_cost <- function() {
    return(is.numeric(control_cost))
  }

  # Get control cost results label
  self$get_cost_label <- function(include_id = TRUE) {
    return(paste0(self$get_label(include_id = include_id), "_cost"))
  }

  # Get the unit of control cost
  self$get_cost_unit <- function() {
    return(attr(control_cost, "unit"))
  }

  # Get attached attributes
  self$get_attributes <- function(n) {
    attr_list <- list()
    if (self$get_label() %in% names(attributes(n))) {
      attr_list[[self$get_label()]] <- attr(n, self$get_label())
    }
    if ("undetected" %in% names(attributes(n))) {
      attr_list[["undetected"]] <- attr(n, "undetected")
    }
    if (self$get_cost_label() %in% names(attributes(n))) {
      attr_list[[self$get_cost_label()]] <- attr(n, self$get_cost_label())
    }
    return(attr_list)
  }

  # Clear attached attributes
  self$clear_attributes <- function(n) {
    attr(n, self$get_label()) <- NULL
    attr(n, "undetected") <- NULL
    attr(n, self$get_cost_label()) <- NULL
    return(n)
  }

  # Control apply method
  self$apply <- function(n, tm) {

    # Search and destroy (combined detect and remove)
    if (control_type == "search_destroy") {

      # Initial zero removals
      controlled <- list(detected = as.numeric(n*0),
                         undetected = as.numeric(n*0))
      if (population_model$get_type() == "stage_structured") {
        controlled <- lapply(controlled, function(controlled_i) {
          controlled_i <- array(controlled_i, dim(n))
          colnames(controlled_i) <-
            attr(population_model$get_growth(), "labels")
          return(controlled_i)
        })
      }

      # Scheduled time step?
      if (is.null(schedule) || tm %in% schedule) {

        # Partition n into detected and undetected
        if ("undetected" %in% names(attributes(n))) {
          n_apply <- list(detected = as.numeric(n - attr(n, "undetected")),
                          undetected = as.numeric(attr(n, "undetected")))
        } else {
          n_apply <- list(detected = as.numeric(n*0),
                          undetected = as.numeric(n))
        }

        # Occupied locations
        idx <- which(rowSums(
          as.matrix(n)[,self$get_stages(), drop = FALSE]) > 0)
        if (length(idx) > 0) {

          # Get control effectiveness (probability)
          control_pr <- manage_pr[idx]

          # Transform population level effectiveness values
          if (manage_pr_type == "population" &&
              population_model$get_type() %in% c("unstructured",
                                                 "stage_structured")) {

            # Controllable occupied population sizes
            if (population_model$get_type() == "stage_structured") {
              n_idx <- rowSums(n[idx, self$get_stages(), drop = FALSE])
            } else {
              n_idx <- n[idx]
            }

            # Transform overall to individual probabilities
            control_pr <- control_pr^(1/n_idx)
          }

          # Sample number controlled and remove when search & destroy
          if (population_model$get_type() == "stage_structured") {
            n_apply <- lapply(n_apply, function(a) array(a, dim(n)))
            for (i in c("detected", "undetected")) {
              for (j in self$get_stages()) {
                controlled[[i]][idx, j] <-
                  stats::rbinom(length(idx), size = n_apply[[i]][idx, j],
                                prob = control_pr)
              }
            }
            if ("undetected" %in% names(attributes(n))) {
              attr(n, "undetected")[idx,] <-
                attr(n, "undetected")[idx,] - controlled$undetected[idx,]
            } else {
              attr(n, "undetected") <- n[,] - controlled$undetected
            }
            controlled <- controlled$detected + controlled$undetected
            n[idx,] <- n[idx,] - controlled[idx,]
          } else {
            for (i in c("detected", "undetected")) {
              controlled[[i]][idx] <-
                stats::rbinom(length(idx), size = n_apply[[i]][idx],
                              prob = control_pr)
            }
            if (population_model$get_type() == "presence_only") {
              controlled <- lapply(controlled, as.logical)
              if ("undetected" %in% names(attributes(n))) {
                attr(n, "undetected")[idx] <-
                  attr(n, "undetected")[idx] & !controlled$undetected[idx]
              } else {
                attr(n, "undetected") <- as.logical(n) & !controlled$undetected
              }
              controlled <- controlled$detected | controlled$undetected
              n[idx] <- n[idx] & !controlled[idx]
            } else {
              if ("undetected" %in% names(attributes(n))) {
                attr(n, "undetected")[idx] <-
                  attr(n, "undetected")[idx] - controlled$undetected[idx]
              } else {
                attr(n, "undetected") <- as.numeric(n) - controlled$undetected
              }
              controlled <- controlled$detected + controlled$undetected
              n[idx] <- n[idx] - controlled[idx]
            }
          }
        } else {
          controlled <- controlled[[1]] # zeros
        }
      } else {
        controlled <- controlled[[1]] # zeros
      }

      # Attach all undetected when zero controlled
      if (!"undetected" %in% names(attributes(n))) {
        if (population_model$get_type() == "stage_structured") {
          attr(n, "undetected") <- n[,]
        } else if (population_model$get_type() == "presence_only") {
          attr(n, "undetected") <- as.logical(n)
        } else {
          attr(n, "undetected") <- as.numeric(n)
        }
      }

      # Attach control
      if (population_model$get_type() == "presence_only") {
        attr(n, self$get_label()) <- as.logical(controlled)
      } else {
        if (population_model$get_type() == "stage_structured") {
          attr(n, self$get_label()) <-
            (rowSums(controlled) > 0 &
               rowSums(n[,self$get_stages(), drop = FALSE]) == 0)
        } else {
          attr(n, self$get_label()) <- (controlled > 0 & n == 0)
        }
        if (manage_pr_type == "individual") {
          attr(attr(n, self$get_label()), "number") <- controlled
        }
      }

      # Attach control costs as an attribute via label
      if (!is.null(control_cost)) {
        if (is.null(schedule) || tm %in% schedule) {
          attr(n, self$get_cost_label()) <- control_cost
        } else {
          attr(n, self$get_cost_label()) <- control_cost*0
        }
      }
    }

    # Growth, spread, or establishment control
    if (control_type %in% c("growth", "spread", "establishment")) {

      # Initial no application locations
      n_apply <- rep(FALSE, region$get_locations())

      # Initial no control cost locations
      if (!is.null(control_cost)) {
        cost_apply <- rep(FALSE, region$get_locations())
      }

      # Scheduled time step?
      if (is.null(schedule) || tm %in% schedule) {

        # Apply control to existing/known/scheduled treatment locations
        if (!is.null(exist_control)) {

          # Occupied locations
          idx <- which(rowSums(
            as.matrix(n)[,self$get_stages(), drop = FALSE]) > 0)
          if (length(idx) > 0) {

            # Get control application locations
            n_apply[idx] <- exist_control[idx] > 0
          }

          # Control cost locations
          if (!is.null(control_cost)) {
            cost_apply <- exist_control > 0
          }
        }

        # Detection-based control
        if ("undetected" %in% names(attributes(n))) {

          # Detected locations
          idx <- which(rowSums(as.matrix(n - attr(n, "undetected"))) > 0)

          # Expand control locations via radius
          if (is.numeric(radius) && length(idx) > 0 &&
              region$get_type() %in% c("grid", "patch")) {
            idx <- region$get_nearby(idx, radius)
            if (control_type %in% c("growth", "spread")) { # occupied only
              idx <- idx[which(rowSums(
                as.matrix(n)[,self$get_stages(), drop = FALSE])[idx] > 0)]
            }
          }

          # Set/update control locations
          if (length(idx) > 0) {
            n_apply[idx] <- TRUE
          }

          # Set/update control cost locations
          if (!is.null(control_cost) && length(idx) > 0) {
            cost_apply[idx] <- TRUE
          }
        }
      }

      # Attach control as attributes to process in growth or spread models
      attr(n, self$get_label()) <- n_apply*control_mult + !n_apply
      if (population_model$get_type() == "stage_structured") {
        if (control_type == "growth") {
          attr(attr(n, self$get_label()), "stages") <- self$get_stages()
          attr(attr(n, self$get_label()), "apply_to") <- apply_to
        }
      }

      # Attach control costs as an attribute via label
      if (!is.null(control_cost)) {
        if (is.null(schedule) || tm %in% schedule) {
          if (region$spatially_implicit()) { # cost/m2
            if (is.numeric(attr(n, "diffusion_radius"))) {
              total_area <- pi*(attr(n, "diffusion_radius"))^2
            } else if (is.numeric(attr(n, "spread_area"))) {
              total_area <- attr(n, "spread_area")
            } else {
              stop(paste("Cannot calculate spatially implicit control cost",
                         "without area occupied."), call. = FALSE)
            }
            attr(n, self$get_cost_label()) <-
              control_cost*total_area*cost_apply
          } else {
            attr(n, self$get_cost_label()) <- control_cost*cost_apply
          }
        } else {
          attr(n, self$get_cost_label()) <- control_cost*0
        }
      }
    }

    return(n)
  }

  return(self)
}
