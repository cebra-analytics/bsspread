#' Removals class builder
#'
#' Builds a class for simulating the application of management removals or
#' eradication of an invasive species.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the population spread model simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the population spread model
#'   simulations.
#' @param removal_pr A single probability value (0-1), or vector of values for
#'   each location specified by the \code{region}, to represent the likelihood
#'   of removal at locations where the invasive species has been detected, as
#'   indicated via a population attribute when present, else removal is applied
#'   at all locations. Default is \code{1}, indicating all (detected)
#'   occurrences are removed.
#' @param removal_pr_type Configures how the removal probability values are
#'   specified. May be specified at either: the code{"individual"} level
#'   (default), whereby removal probability values denote the probability of
#'   removal for each (removable) invasive individual; the \code{"population"}
#'   level, whereby removal probability values denote the probability of
#'   removing all (removable) individuals within the (local) invasive
#'   population.
#' @param remove_always A logical indication of whether or not removal is
#'   always to applied to locations where invasive species are present (even
#'   when they are not detected by explicit surveillance). Default is
#'   \code{FALSE}, indicating that removal is dependent on detection when
#'   surveillance actions are present. Set to \code{TRUE} for locations where
#'   some removal is likely regardless of explicit management efforts (e.g.
#'   via the general public or property owners).
#' @param detected_only A logical indication of whether or not removal is only
#'   applied to detected individuals (e.g. via traps). Default is \code{FALSE},
#'   indicating that removal is applied to all individuals at locations where
#'   the invasive species has been detected (e.g. via treatment).
#' @param removal_cost Numeric vector of distributed removal costs (combined
#'   resource and fixed costs) or a single cost value for each location where
#'   removal is applied. For spatially-implicit area-based regions, cost
#'   should be specified as cost per metres squared. Costs are accumulated for
#'   each application of the removal at each (scheduled) simulation time step.
#'   The cost unit may be added as an attribute
#'   (\code{attr(removal_cost, "unit")}). Default is \code{NULL} when costs
#'   are unavailable.
#' @param radius Optional radius (m) of the removal. Stochastic removal is
#'   applied to all locations within the specified radius of each location
#'   where the invasive species has been detected. Default is \code{NULL},
#'   indicating that removal is only applied at detected locations (when
#'   specified via a population attribute).
#' @param stages Numeric vector of population stages (indices) to which
#'   removals are applied. Default is all stages (when set to \code{NULL}).
#' @param schedule Vector of discrete simulation time steps (t = 0, 1, 2, ...)
#'   in which to apply removals. Default is all time steps (when set to
#'   \code{NULL}).
#' @param ... Additional parameters.
#' @return A \code{Removals} class object (list) containing a function
#'   for accessing attributes and applying simulated removals:
#'   \describe{
#'     \item{\code{get_type()}}{Get the type of action ("removal").}
#'     \item{\code{get_id()}}{Get the actions numeric identifier.}
#'     \item{\code{set_id(id)}}{Set the actions numeric identifier. Should be
#'       an integer >= 1).}
#'     \item{\code{get_label(include_id = TRUE)}}{Get the actions label used
#'       in simulation results ("<id>_removed" or "removed"). Set
#'       \code{include_id} to include the action \code{id} as a label prefix
#'       (default is \code{TRUE}).}
#'     \item{\code{get_removal_pr_type()}}{Get the removal probability type.}
#'     \item{\code{get_stages()}}{Get the population stages to which removals
#'       are applied.}
#'     \item{\code{get_schedule()}}{Get the scheduled simulation time steps in
#'       which removals are applied.}
#'     \item{\code{include_cost()}}{Logical indication of a cost parameter
#'       having a value.}
#'     \item{\code{get_cost_label(include_id = TRUE)}}{Get the removal
#'       cost label used in simulation results ("<id>_removal_cost" or
#'       "removal_cost"). Set \code{include_id} to include the action \code{id}
#'       as a label prefix (default is \code{TRUE}).}
#'     \item{\code{get_cost_unit()}}{Get the unit of removal cost.}
#'     \item{\code{get_attributes(n)}}{Get attached attributes associated
#'       with this action from a simulated population vector or matrix
#'       \code{n}, and return a list of the attached attributes.}
#'     \item{\code{clear_attributes(n)}}{Clear attached attributes associated
#'       with this action from a simulated population vector or matrix
#'       \code{n}, and return \code{n} without the attached attributes.}
#'     \item{\code{apply(n, tm)}}{Apply removals to a simulated population
#'       vector or matrix \code{n}, potentially with attached attributes
#'       relating to previously applied actions, providing the time step
#'       \code{tm} is in the \code{schedule}, and return the resulting
#'       population \code{n} along with attached attributes relating to the
#'       newly applied removals.}
#'   }
#' @export
Removals <- function(region, population_model,
                     removal_pr = 1,
                     removal_pr_type = c("individual",
                                         "population"),
                     remove_always = FALSE,
                     detected_only = FALSE,
                     removal_cost = NULL,
                     radius = NULL,
                     stages = NULL, schedule = NULL, ...) {
  UseMethod("Removals")
}

#' @name Removals
#' @export
Removals.Region <- function(region, population_model,
                            removal_pr = 1,
                            removal_pr_type = c("individual",
                                                "population"),
                            remove_always = FALSE,
                            detected_only = FALSE,
                            removal_cost = NULL,
                            radius = NULL,
                            stages = NULL, schedule = NULL, ...) {

  # Build via base class
  self <- Actions(region = region,
                  population_model = population_model,
                  type = "removal",
                  stages = stages,
                  schedule = schedule,
                  class = "Removals")

  # Validate removal probability
  if (is.null(removal_pr) ||
      (!is.null(removal_pr) &&
       (!is.numeric(removal_pr) ||
        any(removal_pr < 0) || any(removal_pr > 1) ||
        !length(removal_pr) %in% c(1, region$get_locations())))) {
    stop(paste("Removal probability should be a vector with a value 0-1 for",
               "each region location."), call. = FALSE)
  }
  if (length(removal_pr) == 1) {
    removal_pr <- rep(removal_pr, region$get_locations())
  }

  # Check removal probability type
  removal_pr_type <- match.arg(removal_pr_type)

  # Validate radius
  if (!is.null(radius) && (!is.numeric(radius) || radius < 0)) {
    stop("The radius (m) parameter must be numeric and >= 0.", call. = FALSE)
  }

  # Notify if radius is provided when detected only
  if (is.numeric(radius) && detected_only) {
    message("Radius is not used when only detected individuals are removed.")
  }

  # Check and process removal cost
  if (!is.null(removal_cost)) {
    if (!is.numeric(removal_cost) ||
        !length(removal_cost) %in% c(1, region$get_locations())) {
      stop(paste("The removal cost parameter must be a numeric vector with",
                 "values for each location."), call. = FALSE)
    }
    cost_unit <- attr(removal_cost, "unit")
    if (length(removal_cost) == 1) {
      removal_cost <- removal_cost*(removal_pr > 0)
    }
    attr(removal_cost, "unit") <- cost_unit
  }

  # Get results label
  self$get_label <- function(include_id = TRUE) {
    if (!is.null(self$get_id()) && include_id) {
      return(paste0(self$get_id(), "_", "removed"))
    } else {
      return("removed")
    }
  }

  # Get the removal probability type
  self$get_removal_pr_type <- function() {
    return(removal_pr_type)
  }

  # Does the removal cost parameter have a value?
  self$include_cost <- function() {
    return(is.numeric(removal_cost))
  }

  # Get removal cost results label
  self$get_cost_label <- function(include_id = TRUE) {
    if (!is.null(self$get_id()) && include_id) {
      return(paste0(self$get_id(), "_", "removal_cost"))
    } else {
      return("removal_cost")
    }
  }

  # Get the unit of removal cost
  self$get_cost_unit <- function() {
    return(attr(removal_cost, "unit"))
  }

  # Get attached attributes
  self$get_attributes <- function(n) {
    attr_list <- list()
    if (self$get_label() %in% names(attributes(n))) {
      attr_list[[self$get_label()]] <- attr(n, self$get_label())
    }
    if (self$get_cost_label() %in% names(attributes(n))) {
      attr_list[[self$get_cost_label()]] <- attr(n, self$get_cost_label())
    }
    return(attr_list)
  }

  # Clear attached attributes
  self$clear_attributes <- function(n) {
    attr(n, self$get_label()) <- NULL
    attr(n, self$get_cost_label()) <- NULL
    return(n)
  }

  # Removal apply method
  self$apply <- function(n, tm) {

    # Initial zero removals
    removed <- list(detected = as.numeric(n*0), undetected = as.numeric(n*0))
    if (population_model$get_type() == "stage_structured") {
      removed <- lapply(removed, function(removed_i) {
        removed_i <- array(removed_i, dim(n))
        colnames(removed_i) <- attr(population_model$get_growth(), "labels")
        return(removed_i)
      })
    }

    # Initially no removal cost locations
    if (!is.null(removal_cost)) {
      cost_apply <- rep(FALSE, region$get_locations())
    }

    # Scheduled time step?
    if (is.null(schedule) || tm %in% schedule) {

      # Detection-based removal
      if (!remove_always && "undetected" %in% names(attributes(n))) {

        # Individuals to which to apply removal
        if (detected_only) {
          n_apply <- list(detected = as.numeric(n - attr(n, "undetected")),
                          undetected = as.numeric(n*0))
        } else {
          n_apply <- list(detected = as.numeric(n - attr(n, "undetected")),
                          undetected = as.numeric(attr(n, "undetected")))
        }

        # Removal locations
        if (population_model$get_type() == "stage_structured") {
          idx <-
            which(rowSums((n - attr(n, "undetected"))[,self$get_stages(),
                                                      drop = FALSE]) > 0)
        } else {
          idx <- which(n_apply$detected > 0)
        }

      } else {

        # Apply to all individuals
        if ("undetected" %in% names(attributes(n))) {
          n_apply <- list(detected = as.numeric(n - attr(n, "undetected")),
                          undetected = as.numeric(attr(n, "undetected")))
        } else {
          n_apply <- list(detected = as.numeric(n*0),
                          undetected = as.numeric(n))
        }

        # Removal locations
        if (!remove_always && detected_only) {
          idx <- c() # none detected
        } else {
          if (population_model$get_type() == "stage_structured") {
            idx <- which(rowSums(n[,self$get_stages(), drop = FALSE]) > 0)
          } else {
            idx <- which(n > 0)
          }
        }
      }

      # Expand removal locations via radius
      if ("undetected" %in% names(attributes(n)) && is.numeric(radius) &&
          length(idx) > 0 && region$get_type() %in% c("grid", "patch")) {
        idx <- region$get_nearby(idx, radius)
        idx <- idx[which(rowSums(
          as.matrix(n)[,self$get_stages(), drop = FALSE])[idx] > 0)]
      }

      # Sample and apply removals
      if (length(idx) > 0) {

        # Removal probabilities at locations
        removal_pr_idx <- removal_pr[idx]

        # Transform population level effectiveness values
        if (removal_pr_type == "population" &&
            population_model$get_type() %in% c("unstructured",
                                               "stage_structured")) {

          # Removable occupied population sizes
          if (population_model$get_type() == "stage_structured") {
            n_idx <- rowSums(n[idx, self$get_stages(), drop = FALSE])
          } else {
            n_idx <- n[idx]
          }

          # Transform overall to individual probabilities
          removal_pr_idx <- removal_pr_idx^(1/n_idx)
        }

        # Sample and apply
        if (population_model$get_type() == "stage_structured") {
          n_apply <- lapply(n_apply, function(a) array(a, dim(n)))
          for (i in c("detected", "undetected")) {
            for (j in self$get_stages()) {
              removed[[i]][idx, j] <-
                stats::rbinom(length(idx), size = n_apply[[i]][idx, j],
                              prob = removal_pr_idx)
            }
          }
          if ("undetected" %in% names(attributes(n))) { # update
            attr(n, "undetected")[idx,] <-
              attr(n, "undetected")[idx,] - removed$undetected[idx,]
          }
          removed <- removed$detected + removed$undetected
          n[idx,] <- n[idx,] - removed[idx,]
        } else {
          for (i in c("detected", "undetected")) {
            removed[[i]][idx] <-
              stats::rbinom(length(idx), size = n_apply[[i]][idx],
                            prob = removal_pr_idx)
          }
          if (population_model$get_type() == "presence_only") {
            removed <- lapply(removed, as.logical)
            if ("undetected" %in% names(attributes(n))) { # update
              attr(n, "undetected")[idx] <-
                attr(n, "undetected")[idx] & !removed$undetected[idx]
            }
            removed <- removed$detected | removed$undetected
            n[idx] <- n[idx] & !removed[idx]
          } else {
            if ("undetected" %in% names(attributes(n))) { # update
              attr(n, "undetected")[idx] <-
                attr(n, "undetected")[idx] - removed$undetected[idx]
            }
            removed <- removed$detected + removed$undetected
            n[idx] <- n[idx] - removed[idx]
          }
        }
      } else {
        removed <- removed[[1]] # zeros
      }

      # Set removal cost locations
      if (!is.null(removal_cost)) {
        if (length(idx) > 0) {
          cost_apply[idx] <- TRUE
        }
      }
    } else {
      removed <- removed[[1]] # zeros
    }

    # Attach removed as an attribute
    if (population_model$get_type() == "presence_only") {
      attr(n, self$get_label()) <- as.logical(removed)
    } else {
      if (population_model$get_type() == "stage_structured") {
        attr(n, self$get_label()) <-
          (rowSums(removed) > 0 &
             rowSums(n[,self$get_stages(), drop = FALSE]) == 0)
      } else {
        attr(n, self$get_label()) <- (removed > 0 & n == 0)
      }
      if (removal_pr_type == "individual") {
        attr(attr(n, self$get_label()), "number") <- removed
      }
    }

    # Attach removal costs as an attribute
    if (!is.null(removal_cost)) {
      if (population_model$get_region()$spatially_implicit()) { # cost/m2
        if (is.numeric(attr(n, "diffusion_radius"))) {
          total_area <- pi*(attr(n, "diffusion_radius"))^2
        } else if (is.numeric(attr(n, "spread_area"))) {
          total_area <- attr(n, "spread_area")
        } else {
          stop(paste("Cannot calculate spatially implicit removal cost",
                     "without area occupied."), call. = FALSE)
        }
        attr(n, self$get_cost_label()) <- removal_cost*total_area*cost_apply
      } else {
        attr(n, self$get_cost_label()) <- removal_cost*cost_apply
      }
    }

    return(n)
  }

  return(self)
}
