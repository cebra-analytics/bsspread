#' Detection class builder
#'
#' Builds a class for simulating the application of detection or surveillance
#' applied in the management of an invasive species.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the population spread model simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the population spread model
#'   simulations.
#' @param sensitivity Numeric vector of distributed surveillance sensitivities,
#'    or probability of detection at each location.
#' @param sensitivity_type Configures how the surveillance detection
#'   sensitivity values are specified. May be specified at either: the
#'   \code{"individual"} level (default), whereby sensitivity values denote the
#'   probability of detection for each (detectable) invasive individual; the
#'   \code{"population"} size (threshold) level, whereby sensitivity values
#'   denote the probability of detecting (local) populations having at
#'   least at certain size (of detectable individuals), specified via a
#'   threshold value (below), and the probability of detection for (local)
#'   population sizes below the threshold value is reduced proportionally (i.e.
#'   scaled by size/threshold); the \code{"presence"}/absence level, whereby
#'   sensitivity values denote the probability of detection any (local)
#'   presence of (detectable individuals within) the invasive species
#'   population (equivalent to threshold = 1).
#' @param sensitivity_threshold The threshold (minimum) population size (of
#'   detectable individuals) for (locally) applying surveillance detection
#'   probability (sensitivity) values when \code{sensitivity_type} is
#'   \code{"population"}. The probability of detection for (local) population
#'   sizes below the threshold value is reduced proportionally (i.e. scaled by
#'   size/threshold). Default is \code{NULL} when \code{sensitivity_type} is
#'   \code{"population"}, or is set to \code{1} when \code{sensitivity_type} is
#'   \code{"presence"}.
#' @param surv_cost Numeric vector of distributed surveillance costs (combined
#'   resource and fixed costs) or a single cost value for each location where
#'   surveillance is applied. For spatially-implicit area-based regions, cost
#'   should be specified as cost for the entire spatially-implicit area. Costs
#'   are accumulated for each application of the surveillance at each
#'   (scheduled) simulation time step. The cost unit may be added as an
#'   attribute (\code{attr(surv_cost, "unit")}). Default is \code{NULL} when
#'   costs are unavailable.
#' @param stages Numeric vector of population stages (indices) which the
#'   surveillance detects. Default is all stages (when set to \code{NULL}).
#' @param schedule Vector of discrete simulation time steps (t = 0, 1, 2, ...)
#'   in which to apply surveillance. Default is all time steps (when set to
#'   \code{NULL}).
#' @param ... Additional parameters.
#' @return A \code{Detection} class object (list) containing a function
#'   for accessing attributes and applying simulated detection:
#'   \describe{
#'     \item{\code{get_type()}}{Get the type of action ("detection").}
#'     \item{\code{get_id()}}{Get the actions identifier.}
#'     \item{\code{set_id(id)}}{Set the actions identifier.}
#'     \item{\code{get_label(include_id = TRUE)}}{Get the actions label used
#'       in simulation results ("<id>_detected" or "detected").
#'       Set \code{include_id} to include the action \code{id} as a label
#'       prefix (default is \code{TRUE}).}
#'     \item{\code{get_sensitivity()}}{Get the surveillance sensitivity.}
#'     \item{\code{get_sensitivity_type()}}{Get the surveillance sensitivity
#'       type.}
#'     \item{\code{get_stages()}}{Get the population stages which the
#'       surveillance detects.}
#'     \item{\code{get_schedule()}}{Get the scheduled simulation time steps in
#'       which surveillance is applied.}
#'     \item{\code{include_cost()}}{Logical indication of a cost parameter
#'       having a value.}
#'     \item{\code{get_cost_label(include_id = TRUE)}}{Get the surveillance
#'       cost label used in simulation results ("<id>_surv_cost" or
#'       "surv_cost"). Set \code{include_id} to include the action \code{id}
#'       as a label prefix (default is \code{TRUE}).}
#'     \item{\code{get_cost_unit()}}{Get the unit of surveillance cost.}
#'     \item{\code{clear_attributes(n)}}{Clear attached attributes associated
#'       with this action from a simulated population vector or matrix
#'       \code{n}, and return \code{n} without the attached attributes.}
#'     \item{\code{apply(n, tm)}}{Apply surveillance detection to a simulated
#'       population vector or matrix \code{n}, potentially with attached
#'       attributes relating to previously applied actions, providing the time
#'       step \code{tm} is in the \code{schedule}, and return the resulting
#'       population \code{n} along with attached attributes relating to the
#'       newly applied detection/surveillance.}
#'   }
#' @export
Detection <- function(region,
                      population_model,
                      sensitivity,
                      sensitivity_type = c("individual",
                                           "population",
                                           "presence"),
                      sensitivity_threshold = NULL,
                      surv_cost = NULL,
                      stages = NULL,
                      schedule = NULL, ...) {
  UseMethod("Detection")
}

#' @name Detection
#' @export
Detection.Region <- function(region,
                             population_model,
                             sensitivity,
                             sensitivity_type = c("individual",
                                                  "population",
                                                  "presence"),
                             sensitivity_threshold = NULL,
                             surv_cost = NULL,
                             stages = NULL,
                             schedule = NULL, ...) {

  # Build via base class
  self <- Actions(region = region,
                  population_model = population_model,
                  type = "detection",
                  stages = stages,
                  schedule = schedule,
                  class = "Detection")

  # Check and process the sensitivity
  if (!is.numeric(sensitivity) ||
      !length(sensitivity) %in% c(1, region$get_locations())) {
    stop(paste("The surveillance sensitivity parameter must be a numeric",
               "vector with values for each location."), call. = FALSE)
  }
  if (any(sensitivity < 0) || any(sensitivity > 1)) {
    stop("Surveillance sensitivity values should be <= 0 and <= 1.",
         call. = FALSE)
  }
  if (length(sensitivity) == 1) {
    sensitivity <- rep(sensitivity, region$get_locations())
  }

  # Check and process sensitivity type and threshold
  sensitivity_type <- match.arg(sensitivity_type)
  if (!is.null(sensitivity_threshold) &&
      (!is.numeric(sensitivity_threshold) || sensitivity_threshold < 1)) {
    stop(paste("The sensitivity population size threshold parameter should be",
               "a numeric value > 0."),
         call. = FALSE)
  }
  if (sensitivity_type == "individual" && is.numeric(sensitivity_threshold)) {
    message(paste("Ignoring the sensitivity population size threshold value,",
                  "it is not used for a individual level sensitivity type."))
  } else if (sensitivity_type == "population" &&
             population_model$get_type() %in% c("unstructured",
                                                "stage_structured") &&
             !is.numeric(sensitivity_threshold)) {
    stop(paste("A sensitivity population size threshold is required for a",
               "population level sensitivity type."),
         call. = FALSE)
  } else if (sensitivity_type == "presence") {
    if (is.numeric(sensitivity_threshold) && sensitivity_threshold != 1) {
      message(paste("The sensitivity population size threshold value of",
                    sensitivity_threshold, "has been set to 1 for a presence",
                    "level sensitivity type."))
    }
    sensitivity_threshold <- 1
  }

  # Check and process surveillance cost
  if (!is.null(surv_cost)) {
    if (!is.numeric(surv_cost) ||
        !length(surv_cost) %in% c(1, region$get_locations())) {
      stop(paste("The surveillance cost parameter must be a numeric vector",
                 "with values for each location."), call. = FALSE)
    }
    cost_unit <- attr(surv_cost, "unit")
    if (length(surv_cost) == 1) {
      surv_cost <- surv_cost*(sensitivity > 0)
    }
    attr(surv_cost, "unit") <- cost_unit
  }

  # Get results label
  self$get_label <- function(include_id = TRUE) {
    if (!is.null(self$get_id()) && include_id) {
      return(paste0(self$get_id(), "_", "detected"))
    } else {
      return("detected")
    }
  }

  # Get the surveillance sensitivity
  self$get_sensitivity <- function() {
    return(sensitivity)
  }

  # Get the surveillance sensitivity type
  self$get_sensitivity_type <- function() {
    return(sensitivity_type)
  }

  # Does the surveillance cost parameter have a value?
  self$include_cost <- function() {
    return(is.numeric(surv_cost))
  }

  # Get surveillance cost results label
  self$get_cost_label <- function(include_id = TRUE) {
    if (!is.null(self$get_id()) && include_id) {
      return(paste0(self$get_id(), "_", "surv_cost"))
    } else {
      return("surv_cost")
    }
  }

  # Get the unit of surveillance cost
  self$get_cost_unit <- function() {
    return(attr(surv_cost, "unit"))
  }

  # Clear attached attributes
  self$clear_attributes <- function(n) {
    attr(n, self$get_label()) <- NULL
    attr(n, "undetected") <- NULL
    attr(n, self$get_cost_label()) <- NULL
    return(n)
  }

  # Detection/surveillance apply method
  self$apply <- function(n, tm) {

    # Initial zero detections
    detected <- as.numeric(n)*0
    if (population_model$get_type() == "stage_structured") {
      detected <- array(detected, dim(n))
      colnames(detected) <- attr(population_model$get_growth(), "labels")
    }

    # Only survey undetected population
    if (!is.null(attr(n, "undetected"))) {
      undetected <- attr(n, "undetected")
    } else {
      undetected <- detected # zeros
      undetected[] <- n[]
    }

    # Scheduled time step?
    if (is.null(schedule) || tm %in% schedule) {

      # Occupied locations
      idx <- which(rowSums(as.matrix(undetected)) > 0)
      if (length(idx) > 0) {

        # Get detection sensitivity (probability)
        detect_pr <- sensitivity[idx]

        # Transform population/presence level sensitivities
        if (sensitivity_type %in% c("population", "presence") &&
            population_model$get_type() %in% c("unstructured",
                                               "stage_structured")) {

          # Detectable occupied population sizes
          if (population_model$get_type() == "stage_structured") {
            n_idx <- rowSums(n[idx, self$get_stages(), drop = FALSE])
          } else {
            n_idx <- n[idx]
          }

          # Apply sensitivity threshold
          detect_pr <- detect_pr*(pmin(n_idx, sensitivity_threshold)/
                                    sensitivity_threshold)

          # Transform overall to individual probabilities
          detect_pr <- 1 - (1 - detect_pr)^(1/n_idx)
        }

        # Sample detections
        if (population_model$get_type() == "stage_structured") {
          for (i in self$get_stages()) {
            detected[idx,i] <- stats::rbinom(length(idx),
                                             size = undetected[idx,i],
                                             prob = detect_pr)
          }
          undetected[idx,] <- undetected[idx,] - detected[idx,]
        } else {
          detected[idx] <- stats::rbinom(length(idx),
                                         size = undetected[idx],
                                         prob = detect_pr)
          undetected[idx] <- undetected[idx] - detected[idx]
        }
      }

      # Attach surveillance costs as an attribute
      if (!is.null(surv_cost)) {
        attr(n, self$get_cost_label()) <- surv_cost
      }

    } else {

      # Attach zero surveillance costs as an attribute
      if (!is.null(surv_cost)) {
        attr(n, self$get_cost_label()) <- surv_cost*0
      }
    }

    # Attach detected and undetected as attributes
    if (population_model$get_type() == "stage_structured") {
      attr(n, self$get_label()) <- as.logical(rowSums(detected))
    } else {
      attr(n, self$get_label()) <- as.logical(detected)
    }
    if (population_model$get_type() == "presence_only") {
      attr(n, "undetected") <- as.logical(undetected)
    } else {
      attr(n, "undetected") <- undetected
      if (sensitivity_type == "individual") {
        attr(attr(n, self$get_label()), "number") <- detected
      }
    }

    return(n)
  }

  return(self)
}
