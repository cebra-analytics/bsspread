#' Impacts class builder
#'
#' Builds a class for calculating incursion impacts during population spread
#' model simulations. Impact values may be monetary, non-monetary, or dynamic.
#' Monetary value impacts specify periodic assets (e.g. crop yields) that loose
#' value at each simulated time step (e.g. crop season) whilst a threat is
#' present or until recovered. Non-monetary value impacts specify ongoing
#' assets (e.g. habitat condition) that reduce in value by a fixed amount
#' whilst a threat is present or until recovered. Dynamic value impacts also
#' specify ongoing assets or resources (e.g. plantation biomass) that
#' increasingly loose value whilst a threat is present and optionally remain
#' impacted until recovered. Dynamic losses in asset value may optionally be
#' linked to apply (proportional) dynamic losses to threat suitability
#' (establishment likelihood), carrying capacity, and/or dispersal attraction
#' when applicable, such as when a threat utilises a resource or asset (e.g.
#' as a food source).
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the population spread model simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the population spread model
#'   simulations.
#' @param impact_type One of \code{"presence"} (default), \code{"density"}, or
#'   \code{"area"}, to indicate if the impact of a threat incursion is
#'   calculated based on presence, density, or area of occupancy respectively.
#'   Density-based impacts may only be used when the threat is modelled via
#'   unstructured or staged-based populations with carrying capacity specified.
#'   Area-based impacts are only applicable to spatially implicit models.
#' @param valuation_type One of \code{"monetary"} (default),
#'   \code{"non-monetary"}, or \code{"dynamic"} to indicate the type of
#'   valuation associated with the impacted asset or resource.
#' @param asset_name A descriptive name of the asset or resource impacted by
#'   the threat incursion.
#' @param asset_value The periodic monetary value per simulation time step
#'   (e.g. crop yield), or the non-monetary value (e.g. habitat condition),
#'   of the asset or resource impacted by a threat incursion. This may be a
#'   spatial layer (\code{terra::SpatRaster} or \code{raster::RasterLayer}),
#'   or a vector of location values, for a grid or network-based model region.
#'   For spatially implicit models, asset value should be expressed as value
#'   per metres squared. Asset value for dynamic value impacts (e.g. plantation
#'   biomass) represents the total value or asset quantity prior to threat
#'   presence, which increasingly looses value when impacted by the threat,
#'   and may be implicitly monetary (with ongoing not periodic value)
#'   or non-monetary, as indicated by the \code{value_unit}.
#' @param value_unit A descriptive unit for asset value quantities. Default is
#'   \code{NULL}, which is changed to "$" for monetary assets.
#' @param loss_rate A fractional loss of asset or resource value resulting from
#'   a threat incursion or presence. For monetary value impacts the fractional
#'   loss is applied to the periodic asset value at every time step that the
#'   threat is present or until recovered. For non-monetary value impacts the
#'   fractional loss to the asset value is maintained during threat presence or
#'   until recovered. For dynamic value impacts the fractional loss is
#'   repeatedly applied to the asset value or resource quantity, thus
#'   increasingly reducing it, whilst the threat is present.
#' @param discount_rate An optional discount rate (per time interval) for
#'   monetary asset values to estimate future values that account
#'   for inflation. Typically the discounting uses market interest rates.
#'   Discounted impact values are calculated by dividing the asset values by
#'   \code{(1 + discount_rate)^time_int} for a specified time interval.
#'   Default is \code{NULL}.
#' @param stages Numeric vector of population stages (indices) that
#'   contribute to impacts when applicable. Default is all stages (when set to
#'   \code{NULL}).
#' @param recovery_delay Number of simulation time steps that incursion impacts
#'   continue to be in effect at previously occupied locations before the
#'   asset value at the locations recover. Only available for spatially
#'   explicit grid or network models. Default \code{NULL} assumes no delay.
#' @param dynamic_links Dynamic value impacts may optionally be linked to apply
#'   (proportional) dynamic losses to threat suitability, carrying capacity,
#'   and/or dispersal attraction when applicable, such as when a threat
#'   utilises a resource or asset (e.g. as a food source). Configure dynamic
#'   links via a vector of one or more strings from \code{"suitability"},
#'   \code{"capacity"}, and/or \code{"attractors"}. Default \code{NULL}
#'   assumes no dynamic links or that impact values are not dynamic.
#' @param ... Additional parameters.
#' @return A \code{Impacts} class object (list) containing functions for
#'   performing impact calculations:
#'   \describe{
#'     \item{\code{get_id()}}{Get the impacts numeric identifier.}
#'     \item{\code{set_id(id)}}{Set the impacts numeric identifier. Should be
#'       an integer >= 1).}
#'     \item{\code{get_impact_type()}}{Get the impact type.}
#'     \item{\code{get_valuation_type()}}{Get the valuation type.}
#'     \item{\code{get_value_unit()}}{Get the asset value unit.}
#'     \item{\code{update_recovery_delay(n)}}{Update the recovery delay
#'       attribute attached to the population vector or matrix \code{n}.}
#'     \item{\code{calculate(n, tm)}}{Perform impact calculations resulting
#'       from incursion population vector or matrix \code{n} at time step
#'       \code{tm}, and return \code{n} with impact values attached.}
#'   }
#' @export
Impacts <- function(region, population_model,
                    impact_type = c("presence", "density", "area"),
                    valuation_type = c("monetary", "non-monetary", "dynamic"),
                    asset_name,
                    asset_value,
                    value_unit = NULL,
                    loss_rate,
                    discount_rate = NULL,
                    stages = NULL,
                    recovery_delay = NULL,
                    dynamic_links = NULL, ...) {

  # Check region and population model objects
  if (!inherits(region, "Region")) {
    stop("Region model must be a 'Region' or inherited class object.",
         call. = FALSE)
  }
  if (!is.null(population_model) &&
      !inherits(population_model, "Population")) {
    stop("Population model must be a 'Population' or inherited class object.",
         call. = FALSE)
  }

  # Match and check/process impact and valuation types
  impact_type <- match.arg(impact_type)
  if (impact_type == "area" && !region$spatially_implicit()) {
    stop(paste("Area-based impacts are only applicable for spatially implicit",
               "regions."), call. = FALSE)
  }
  if (impact_type == "density" &&
      (population_model$get_type() == "presence_only" ||
       is.null(population_model$get_capacity()))) {
    stop(paste("Density-based impacts are only available for unstructured or",
               "stage-based populations with carrying capacity."),
         call. = FALSE)
  }
  valuation_type <- match.arg(valuation_type)

  # Check asset name, value and unit
  if (!is.character(asset_name)) {
    stop("Asset name should be a character string.", call. = FALSE)
  }
  if (!(class(asset_value) %in% c("SpatRaster", "RasterLayer") ||
        is.numeric(asset_value)) ||
      (class(asset_value) %in% c("SpatRaster", "RasterLayer") &&
       (region$get_type() != "grid" ||
        region$get_type() == "grid" && !region$is_compatible(asset_value))) ||
      (is.numeric(asset_value) &&
       region$get_locations() != length(asset_value))) {
    stop(paste("Asset value should be either a spatial layer or numeric",
               "vector compatible with the defined region."), call. = FALSE)
  }
  if (!is.null(value_unit)) {
    if (!is.character(value_unit)) {
      stop("Asset value unit should be a character string.", call. = FALSE)
    }
  } else if (valuation_type == "monetary") {
    value_unit <- "$"
  }

  # Check loss and discount rates
  if (!is.numeric(loss_rate) ||
      (is.numeric(loss_rate) && (loss_rate < 0 || loss_rate > 1))) {
    stop("Loss rate should be numeric, >= 0, and <= 1.", call. = FALSE)
  }
  if (!is.null(discount_rate)) {
    if (!is.numeric(discount_rate) ||
        (is.numeric(discount_rate) &&
         (discount_rate < 0 || discount_rate > 1))) {
      stop("Discount rate should be numeric, >= 0, and <= 1.", call. = FALSE)
    }
    if (valuation_type != "monetary") {
      stop("Discount rate is only applicable for monetary value impacts.",
           call. = FALSE)
    }
  }

  # Check the impact stages
  if (population_model$get_type() %in% c("presence_only", "unstructured")) {
    stages <- 1
  } else if (population_model$get_type() == "stage_structured") {
    if (is.null(stages)) {
      stages <- 1:population_model$get_stages()
    } else if (!is.numeric(stages) ||
               !all(stages %in% 1:population_model$get_stages())) {
      stop(paste("Impact stages must be a vector of stage indices consistent",
                 "with the population model."), call. = FALSE)
    }
  }

  # Check recovery delay
  if (!is.null(recovery_delay) &&
      (!is.numeric(recovery_delay) || recovery_delay < 0)) {
    stop("Recover delay should a number >= 0.", call. = FALSE)
  }

  # Check dynamic links
  if (!is.null(dynamic_links)) {
    if (valuation_type != "dynamic") {
      stop("Dynamic links are only applicable for dynamic value impacts.",
           call. = FALSE)
    }
    if (!is.vector(dynamic_links) ||
        !all(dynamic_links %in% c("suitability", "capacity", "attractors"))) {
      stop(paste("Dynamic links should be a vector of one or more of",
                 "'suitability', 'capacity', and/or 'attractors'."),
           call. = FALSE)
    }
  }

  # Create a class structure
  self <- structure(list(), class = "Impacts")

  # Id for tracking multiple impacts
  id <- 1

  # Get the impacts id
  self$get_id <- function() {
    return(id)
  }

  # Set the impacts id
  self$set_id <- function(id) {
    if (!is.numeric(id) || trunc(id) != id || id < 1) {
      stop("Impacts id should be an integer >= 1.", call. = FALSE)
    }
    id <<- id
  }

  # Get the impact type
  self$get_impact_type <- function() {
    return(impact_type)
  }

  # Get the valuation type
  self$get_valuation_type <- function() {
    return(valuation_type)
  }

  # Get the asset value unit
  self$get_value_unit <- function() {
    return(value_unit)
  }

  # Calculate density-based incursion (internal)
  calculate_density_incursion <- function(n, tm = NULL) {
    n <- rowSums(as.matrix(n)[,stages, drop = FALSE])
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
    return((sum(n[stages]) > 0)*total_area)
  }

  # Update recovery delay
  self$update_recovery_delay <- function(n) {
    if (is.numeric(recovery_delay)) {
      if (is.list(attr(n, "recovery_delay")) &&
          length(attr(n, "recovery_delay")) >= id &&
          is.numeric(attr(n, "recovery_delay")[[id]])) {

        # Update dependent on impact and asset types
        if (impact_type == "presence" ||
            (valuation_type == "dynamic" && impact_type == "density")) {

          # Occurrence and recovery delay locations
          if (impact_type == "density") {
            x <- calculate_density_incursion(n)
          } else { # presence
            x <- rowSums(as.matrix(n)[,stages, drop = FALSE])
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
          delay_idx <- which(attr(n, "recovery_delay")[[id]] == recovery_delay)
          if (any(!idx %in% delay_idx)) {
            new_idx <- idx[which(!idx %in% delay_idx)]
            attr(n, "recovery_delay")[[id]][new_idx] <- recovery_delay
          }

        } else if (impact_type == "density") {

          # Push current density to front of list for first impact only
          if (attr(attr(n, "recovery_delay"), "first") == id) {
            attr(attr(n, "recovery_delay"), "incursions") <-
              c(list(calculate_density_incursion(n)),
                attr(attr(n, "recovery_delay"), "incursions"))
            length(attr(attr(n, "recovery_delay"), "incursions")) <-
              min(length(attr(attr(n, "recovery_delay"), "incursions")),
                  attr(attr(n, "recovery_delay"), "max"))
          }

        } else if (region$spatially_implicit() && impact_type == "area") {

          # Add current area to front of list for first impact only
          if (attr(attr(n, "recovery_delay"), "first") == id) {
            attr(attr(n, "recovery_delay"), "incursions") <-
              c(calculate_area_incursion(n),
                attr(attr(n, "recovery_delay"), "incursions"))
          }

          # Add dynamic multipliers to front of each existing list
          if (valuation_type == "dynamic" &&
              is.list(attr(n, "dynamic_mult")) &&
              length(attr(n, "dynamic_mult")) >= id &&
              is.numeric(attr(n, "dynamic_mult")[[id]]) &&
              is.numeric(attr(attr(n, "recovery_delay")[[id]],
                              "dynamic_mult"))) {
            attr(attr(n, "recovery_delay")[[id]], "dynamic_mult") <-
              c(attr(n, "dynamic_mult")[[id]],
                attr(attr(n, "recovery_delay")[[id]], "dynamic_mult"))
          }
        }

      } else {

        # Initialise recovery delay
        if (!is.list(attr(n, "recovery_delay"))) {
          attr(n, "recovery_delay") <- list()
        }
        if (impact_type == "presence" ||
            (valuation_type == "dynamic" && impact_type == "density")) {
          if (impact_type == "density") {
            x <- calculate_density_incursion(n)
          } else { # presence
            x <- rowSums(as.matrix(n)[,stages, drop = FALSE])
          }
          attr(n, "recovery_delay")[[id]] <- (x > 0)*recovery_delay
        } else { # constant
          attr(n, "recovery_delay")[[id]] <- recovery_delay
          attr(attr(n, "recovery_delay"), "max") <-
            max(attr(attr(n, "recovery_delay"), "max"), recovery_delay)
          if (is.null(attr(attr(n, "recovery_delay"), "first"))) {
            attr(attr(n, "recovery_delay"), "first") <- id
            if (impact_type == "density") {
              attr(attr(n, "recovery_delay"), "incursions") <-
                list(calculate_density_incursion(n))
            } else if (region$spatially_implicit() && impact_type == "area") {
              attr(attr(n, "recovery_delay"), "incursions") <-
                calculate_area_incursion(n)
            }
          }
          if (region$spatially_implicit() && impact_type == "area" &&
              valuation_type == "dynamic" &&
              is.list(attr(n, "dynamic_mult")) &&
              length(attr(n, "dynamic_mult")) >= id &&
              is.numeric(attr(n, "dynamic_mult")[[id]])) {
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
    if (impact_type == "density") {
      x <- calculate_density_incursion(n, tm = tm)
    } else if (region$spatially_implicit() && impact_type == "area") {
      x <- calculate_area_incursion(n)
    } else if (population_model$get_type() == "stage_structured") {
      x <- rowSums(n[,stages, drop = FALSE])
    } else {
      x <- as.numeric(n)
    }

    # Truncate to 1 unless area
    if (impact_type %in% c("presence", "density")) {
      x[which(x > 1)] <- 1
    }

    # Calculate dynamic multipliers
    if (valuation_type == "dynamic") {
      if (is.list(attr(n, "dynamic_mult")) &&
          length(attr(n, "dynamic_mult")) >= id &&
          is.numeric(attr(n, "dynamic_mult")[[id]])) {
        dynamic_mult <- attr(n, "dynamic_mult")[[id]]
      } else {
        dynamic_mult <- 1
      }
      if (impact_type == "area") {
        prev_area <- attr(dynamic_mult, "incursion")
        if (x > 0) {
          if (is.numeric(prev_area) && prev_area < x) {
            prev_loss <- prev_area*(1 - dynamic_mult)
            dynamic_mult <- (x - prev_loss)*(1 - loss_rate)/x
          } else {
            dynamic_mult <- dynamic_mult*(1 - loss_rate)
          }
        }
        attr(dynamic_mult, "incursion") <- x
      } else {
        dynamic_mult <- dynamic_mult*(1 - x*loss_rate)
      }
      dynamic_mult <- as.numeric(unname(dynamic_mult))
    }

    # Recovery delays prolong impacts
    if (impact_type %in% c("presence", "density", "area") &&
        is.list(attr(n, "recovery_delay"))) {
      if (length(attr(n, "recovery_delay")) >= id &&
          is.numeric(attr(n, "recovery_delay")[[id]])) {
        if (impact_type %in% c("presence", "density") &&
            valuation_type == "dynamic") { # decremented delays
          delay <- attr(n, "recovery_delay")[[id]]
          dynamic_mult <- ((x > 0 | delay > 0)*dynamic_mult +
                             (x == 0 & delay == 0))
        } else if (impact_type == "presence") { # decremented delays
          x <- +(x > 0 | attr(n, "recovery_delay")[[id]] > 0)
        } else if (impact_type == "density") {
          delay <- attr(n, "recovery_delay")[[id]]
          prev_incursions <- attr(attr(n, "recovery_delay"), "incursions")
          if (delay > 0 && length(prev_incursions) > 0) {
            for (i in 1:min(length(prev_incursions), delay)) {
              x <- pmax(x, prev_incursions[[i]])
            }
          }
        } else if (impact_type == "area") {
          delay <- attr(n, "recovery_delay")[[id]]
          prev_incursions <- attr(attr(n, "recovery_delay"), "incursions")
          x <- as.numeric(x)
          if (valuation_type == "dynamic" &&
              is.numeric(attr(delay, "dynamic_mult"))) {
            dynamic_incursion <- (1 - (x > 0)*dynamic_mult + (x == 0))*x
            if (delay > 0 && length(prev_incursions) > 0) {
              prev_mult <- attr(delay, "dynamic_mult")
              dynamic_incursion <- ((1 - prev_mult)*prev_incursions)[1:delay]
              if (dynamic_incursion == 0) {
                dynamic_mult <- 1
              }
            }
          } else if (delay > 0 && length(prev_incursions) > 0) {
            x <- max(x, prev_incursions[1:delay], na.rm = TRUE)
          }
        }
      } else {
        x <- as.numeric(x)
        if (impact_type == "presence") {
          x <- +(x > 0)
        }
      }
    } else if (valuation_type == "dynamic" && # no delay on recovery
               is.null(attr(n, "recovery_delay"))) {
      dynamic_mult <- unname((x > 0)*dynamic_mult + (x == 0))
      if (impact_type == "area") {
        dynamic_incursion <- (1 - dynamic_mult)*x
      }
    }

    # Extract spatial raster asset layer values
    if (class(asset_value) %in% c("SpatRaster", "RasterLayer")) {
      asset_value <- asset_value[region$get_indices()][,1]
    }

    # Calculate discount multiplier
    disc_mult <- 1
    if (is.numeric(discount_rate) && is.numeric(tm)) {
      disc_mult <- 1/(1 + discount_rate)^tm
    }

    # Calculate incursion impact
    if (valuation_type == "dynamic") {
      if (impact_type == "area") {
        incursion_impacts <- asset_value*disc_mult*dynamic_incursion
      } else {
        incursion_impacts <- asset_value*disc_mult*(1 - dynamic_mult)
      }
    } else {
      incursion_impacts <- asset_value*loss_rate*x*disc_mult
    }

    # Attach dynamic multipliers
    if (valuation_type == "dynamic") {
      if (valuation_type == "dynamic" && impact_type == "area") {
        attr(dynamic_mult, "incursion") <- x
      }
      attr(n, "dynamic_mult")[[id]] <- dynamic_mult
      if (!is.null(dynamic_links)) {
        attr(attr(n, "dynamic_mult")[[id]], "links") <- dynamic_links
      }
    }

    # Attach calculated impacts to population
    if (!is.list(attr(n, "impacts"))) {
      attr(n, "impacts") <- list()
    }
    attr(n, "impacts")[[id]] <- incursion_impacts

    return(n)
  }

  return(self)
}
