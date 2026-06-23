#' Impacts class builder
#'
#' Builds a class for calculating incursion impacts during population spread
#' model simulations. Encapsulates \code{bsimpact} impact analysis objects.
#'
#' @param impacts A \code{bsimpact::ImpactAnalysis} or inherited class object
#'   specifying how environment, social, and/or economic impacts are
#'   calculated, classified, and/or combined based on asset values,
#'   classifications, and/or incursion loss rates.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the population spread model
#'   simulations.
#' @param impact_stages Numeric vector of population stages (indices) that
#'   contribute to impacts. Default is all stages (when set to \code{NULL}).
#' @param calc_total Logical indicator for whether or not total impacts may
#'   be (sensibly) calculated by summing across region locations. Default is
#'   \code{NULL}, whereby totals are only calculated when the impact (context)
#'   valuation type is \code{"monetary"}. Set to \code{TRUE} if it makes sense
#'   to sum total \code{"non-monetary"} or \code{"ranking"} impact valuations
#'   across region locations (providing multiple aspects can be combined).
#' @param dynamic_links Dynamic impacts (when \code{impacts} is a
#'   \code{bsimpact::ValueImpacts} object with \code{is_dynamic = TRUE}) may
#'   optionally be linked to apply (proportional) dynamic losses to threat
#'   suitability, carrying capacity, and/or dispersal attraction when
#'   applicable, such as when a threat utilises a resource or asset (e.g. as a
#'   food source). Configure dynamic links via a vector of one or more strings
#'   from \code{"suitability"}, \code{"capacity"}, and/or \code{"attractors"}.
#'   Default \code{NULL} assumes no dynamic links or impacts are not dynamic.
#' @param recovery_delay Number of simulation time steps that incursion impacts
#'   continue to be in effect at previously occupied locations before the
#'   asset value at the locations recover. Only available for spatially
#'   explicit grid or network models. Default \code{NULL} assumes no delay.
#' @param ... Additional parameters.
#' @return A \code{Impacts} class object (list) containing functions for
#'   getting the impact context and performing impact calculations:
#'   \describe{
#'     \item{\code{get_impacts()}}{Get \code{bsimpact::ImpactAnalysis} or
#'       inherited class object.}
#'     \item{\code{get_context()}}{Get \code{bsimpact::Context} object.}
#'     \item{\code{get_calc_total()}}{Get calculate total indicator.}
#'     \item{\code{includes_combined()}}{Logical indicator for when impacts are
#'       combined.}
#'     \item{\code{update_recovery_delay(n)}}{Update the recovery delay
#'       attribute attached to the population vector.}
#'     \item{\code{calculate(n, tm)}}{Perform impact calculations resulting
#'       from incursion population vector or matrix \code{n} at time step
#'       \code{tm}, and return \code{n} with a list of impact values attached,
#'       including combined impacts when applicable.}
#'   }
#' @export
Impacts <- function(impacts, population_model,
                          impact_stages = NULL,
                          calc_total = NULL,
                          dynamic_links = NULL,
                          recovery_delay = NULL, ...) {

  # Check the impacts object
  if (!is.null(impacts) && !inherits(impacts, "ImpactAnalysis")) {
    stop("Impacts must be a 'ImpactAnalysis' or inherited class object.",
         call. = FALSE)
  }

  # Check the population model
  if (!is.null(population_model) &&
      !inherits(population_model, "Population")) {
    stop("Population model must be a 'Population' or inherited class object.",
         call. = FALSE)
  }

  # Check the impact stages
  population_type <- population_model$get_type()
  if (population_type == "presence_only" ||
      population_type == "unstructured") {
    impact_stages <- 1
  } else if (population_type == "stage_structured") {
    if (is.null(impact_stages)) {
      impact_stages <- 1:population_model$get_stages()
    } else if (!is.numeric(impact_stages) ||
               !all(impact_stages %in% 1:population_model$get_stages())) {
      stop("Impact stages must be a vector of stage indices consistent ",
           "with the population model.", call. = FALSE)
    }
  }

  # Check and resolve the calculate total indicator
  if (is.null(calc_total)) {
    calc_total <- (impacts$get_context()$get_valuation_type() == "monetary")
  } else if (!is.logical(calc_total)) {
    stop("Calculate total indicator should be logical.", call. = FALSE)
  }

  # Check population capacity is available to calculate density
  if (impacts$get_incursion()$get_type() == "density" &&
      is.null(population_model$get_capacity())) {
    stop("Population capacity is required for density-based impacts.",
         call. = FALSE)
  }

  # Check dynamic links
  if (!is.null(dynamic_links)) {
    if (!is.function(impacts$get_is_dynamic) || !impacts$get_is_dynamic()) {
      stop(paste("Dynamic links are only applicable for dynamic impacts (i.e.",
                 "a bsimpact::ValueImpacts object with is_dynamic = TRUE)."),
           call. = FALSE)
    }
    if (!is.vector(dynamic_links) ||
        !all(dynamic_links %in% c("suitability", "capacity", "attractors"))) {
      stop(paste("Dynamic links should be a vector of one or more of",
                 "'suitability', 'capacity', and/or 'attractors'."),
           call. = FALSE)
    }
  }

  # Check recovery delay
  if (!is.null(recovery_delay) &&
      (!is.numeric(recovery_delay) || recovery_delay < 0)) {
    stop("Recover delay should a number >= 0.", call. = FALSE)
  }

  # Create a class structure
  self <- structure(list(), class = "Impacts")

  # Get impacts object
  self$get_impacts <- function() {
    return(impacts)
  }

  # Get context
  self$get_context <- function() {
    return(impacts$get_context())
  }

  # Get calculate total indicator
  self$get_calc_total <- function() {
    return(calc_total)
  }

  # Impacts combined?
  self$includes_combined <- function() {
    return("combined_impacts" %in% names(impacts))
  }

  # Calculate density-based incursion (internal)
  calculate_density_incursion <- function(n, tm = NULL) {
    n <- rowSums(as.matrix(n)[,impact_stages, drop = FALSE])
    n_density <- n*0
    idx <- which(population_model$get_capacity() > 0)
    n_density[idx] <-
      pmin(n[idx]/population_model$get_capacity(tm = tm)[idx], 1)
    return(n_density)
  }

  # Calculate spatially implicit area-based incursion (internal)
  calculate_area_incursion <- function(n) {

    # Get total area occupied
    if (is.numeric(attr(n, "diffusion_radius"))) {
      total_area <- pi*(attr(n, "diffusion_radius"))^2
    } else if (is.numeric(attr(n, "spread_area"))) {
      total_area <- attr(n, "spread_area")
    } else {
      stop(paste("Cannot calculate spatially implicit impacts without area",
                 "occupied."), call. = FALSE)
    }

    # Return total area occupied when population is present, else zero
    return((sum(n[impact_stages]) > 0)*total_area)
  }

  # Update recovery delay
  self$update_recovery_delay <- function(n) {
    if (is.numeric(recovery_delay)) {
      id <- impacts$get_id()
      if (is.list(attr(n, "recovery_delay")) &&
          length(attr(n, "recovery_delay")) >= id &&
          is.numeric(attr(n, "recovery_delay")[[id]])) {

        # Update for presence incursions
        if (impacts$get_incursion()$get_type() == "presence" ||
            (impacts$get_is_dynamic() &&
             impacts$get_incursion()$get_type() == "density")) { # decremented

          # Occurrence and recovery delay locations
          if (impacts$get_incursion()$get_type() == "density") {
            x <- calculate_density_incursion(n)
          } else { # presence
            x <- rowSums(as.matrix(n)[,impact_stages, drop = FALSE])
          }
          idx <- which(x > 0)
          delay_idx <- which(attr(n, "recovery_delay")[[id]] > 0)

          # Decrement delay where population removed or extirpated
          if (any(!delay_idx %in% idx)) {
            decr_idx <- delay_idx[which(!delay_idx %in% idx)]
            attr(n, "recovery_delay")[[id]][decr_idx] <-
              attr(n, "recovery_delay")[[id]][decr_idx] - 1
          }

          # Set recovery for new occurrences
          if (any(!idx %in% delay_idx)) {
            new_idx <- idx[which(!idx %in% delay_idx)]
            attr(n, "recovery_delay")[[id]][new_idx] <- recovery_delay
          }

        } else if (impacts$get_incursion()$get_type() == "density") {

          # Push current density to front of list for first impact only
          if (attr(attr(n, "recovery_delay"), "first") == id) {
            attr(attr(n, "recovery_delay"), "incursions") <-
              c(list(calculate_density_incursion(n)),
                attr(attr(n, "recovery_delay"), "incursions"))
            length(attr(attr(n, "recovery_delay"), "incursions")) <-
              min(length(attr(attr(n, "recovery_delay"), "incursions")),
                  attr(attr(n, "recovery_delay"), "max"))
          }

        } else if (population_model$get_region()$spatially_implicit() &&
                   impacts$get_incursion()$get_type() == "area") {

          # Add current area to front of list for first impact only
          if (attr(attr(n, "recovery_delay"), "first") == id) {
            attr(attr(n, "recovery_delay"), "incursions") <-
              c(calculate_area_incursion(n),
                attr(attr(n, "recovery_delay"), "incursions"))
          }

          # Add dynamic multipliers to front of each existing list
          if (impacts$get_is_dynamic() && is.list(attr(n, "dynamic_mult")) &&
              length(attr(n, "dynamic_mult")) >= id &&
              is.list(attr(n, "dynamic_mult")[[id]]) &&
              is.list(attr(attr(n, "recovery_delay")[[id]], "dynamic_mult")) &&
              (length(attr(attr(n, "recovery_delay")[[id]],
                          "dynamic_mult")) ==
               length(attr(n, "dynamic_mult")[[id]]))) {
            for (i in 1:length(attr(n, "dynamic_mult")[[id]])){
              attr(attr(n, "recovery_delay")[[id]], "dynamic_mult")[[i]] <-
                c(attr(n, "dynamic_mult")[[id]][[i]],
                  attr(attr(n, "recovery_delay")[[id]], "dynamic_mult")[[i]])
            }
          }
        }

      } else {

        # Initialise recovery delay
        if (!is.list(attr(n, "recovery_delay"))) {
          attr(n, "recovery_delay") <- list()
        }
        if (impacts$get_incursion()$get_type() == "presence" ||
            (impacts$get_is_dynamic() &&
             impacts$get_incursion()$get_type() == "density")) { # decremented
          if (impacts$get_incursion()$get_type() == "density") {
            x <- calculate_density_incursion(n)
          } else { # presence
            x <- rowSums(as.matrix(n)[,impact_stages, drop = FALSE])
          }
          attr(n, "recovery_delay")[[id]] <- (x > 0)*recovery_delay
        } else { # constant
          attr(n, "recovery_delay")[[id]] <- recovery_delay
          attr(attr(n, "recovery_delay"), "max") <-
            max(attr(attr(n, "recovery_delay"), "max"), recovery_delay)
          if (is.null(attr(attr(n, "recovery_delay"), "first"))) {
            attr(attr(n, "recovery_delay"), "first") <- id
            if (impacts$get_incursion()$get_type() == "density") {
              attr(attr(n, "recovery_delay"), "incursions") <-
                list(calculate_density_incursion(n))
            } else if (population_model$get_region()$spatially_implicit() &&
                       impacts$get_incursion()$get_type() == "area") {
              attr(attr(n, "recovery_delay"), "incursions") <-
                calculate_area_incursion(n)
            }
          }
          if (population_model$get_region()$spatially_implicit() &&
              impacts$get_incursion()$get_type() == "area" &&
              impacts$get_is_dynamic() && is.list(attr(n, "dynamic_mult")) &&
              length(attr(n, "dynamic_mult")) >= id &&
              is.list(attr(n, "dynamic_mult")[[id]])) {
            attr(attr(n, "recovery_delay")[[id]], "dynamic_mult") <-
              attr(n, "dynamic_mult")[[id]]
            attr(attr(attr(n, "recovery_delay")[[id]], "dynamic_mult"),
                 "incursion") <- NULL
            attr(attr(attr(n, "recovery_delay")[[id]], "dynamic_mult"),
                 "links") <- NULL
          }
        }
      }
    }
    return(n)
  }

  # Calculate impacts
  self$calculate <- function(n, tm) {

    # Calculate incursion values
    if (impacts$get_incursion()$get_type() == "density") {
      x <- calculate_density_incursion(n, tm = tm)
    } else if (population_model$get_region()$spatially_implicit() &&
               impacts$get_incursion()$get_type() == "area") {
      x <- calculate_area_incursion(n)
    } else if (population_model$get_type() == "stage_structured") {
      x <- rowSums(n[,impact_stages, drop = FALSE])
    } else {
      x <- as.numeric(n)
    }

    # Attach recovery delay and dynamic multipliers
    attr(x, "recovery_delay") <- attr(n, "recovery_delay")
    attr(x, "dynamic_mult") <- attr(n, "dynamic_mult")

    # Set incursion values within impact object
    impacts$get_incursion()$set_values(x)

    # Get incursion impacts
    impact_list <- impacts$incursion_impacts(raw = TRUE, time_int = tm)

    # Append combined impacts when present
    if ("combined_impacts" %in% names(impacts)) {
      impact_list$combined <- impacts$combined_impacts(raw = TRUE)
    }

    # Move attached dynamic multipliers
    if (impacts$get_is_dynamic() &&
        is.list(attr(impact_list, "dynamic_mult"))) {
      if (!is.list(attr(n, "dynamic_mult"))) {
        attr(n, "dynamic_mult") <- list()
      }
      id <- impacts$get_id()
      attr(n, "dynamic_mult")[[id]] <- attr(impact_list, "dynamic_mult")
      if (!is.null(dynamic_links)) {
        attr(attr(n, "dynamic_mult")[[id]], "links") <- dynamic_links
      }
      attr(impact_list, "dynamic_mult") <- NULL
    }

    # Attach calculated impacts to population
    attr(n, "impacts") <- impact_list

    # Update recovery delay
    n <- self$update_recovery_delay(n)

    return(n)
  }

  return(self)
}
