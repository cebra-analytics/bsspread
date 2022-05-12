#' Generic dispersal model class builder
#'
#' Builds a generic class for representing dispersal with functionality for
#' spread simulations. Dispersal may be simulated for presence-only,
#' unstructured or stage-based \code{populations}. For unstructured or
#' stage-based populations, a specified \code{proportion} of the population at
#' each occupied location (cell or patch) is selected (sampled) for dispersal
#' at each simulation time step. Presence-only populations may disperse via
#' specifying a (mean) number of dispersal \code{events}. Dispersal events are
#' generated for each occupied location. Unstructured and stage-based population
#' dispersal may also utilize specified events, otherwise an event is assigned
#' to each dispersing individual. For each dispersal event a destination
#' location is selected via stochastic sampling using (relative) probabilities
#' of dispersal from each occupied location to reachable destinations. These
#' (relative) probabilities are calculated via \code{distance} and/or
#' \code{direction} functions, optionally combined with \code{attractor}
#' (layer) values. Paths to reachable destinations are calculated via the
#' \code{region} class object. Paths may be derived via \code{permeability} or
#' constraint grid-based layers or patch-based data, which are used to adjust
#' (effective) path distances and/or omit unreachable paths. Dispersal paths
#' may also be limited to a \code{maximum} (adjusted) distance. An optional
#' establishment likelihood (layer), which is configured via the population
#' model, may be applied to each dispersal event, resulting in potential
#' "deaths" of individuals or unsuccessful presence-only dispersal events.
#' Presence-only population dispersal may also be configured without
#' \code{events} when \code{proportion} is specified to represent a scaling
#' multiplier for (presumed) actual (rather than relative) dispersal
#' probabilities. Under these circumstances destination locations are sampled
#' from all reachable destinations using the scaled probabilities. The
#' dispersal functionality utilizes a wrapped population list of separate
#' values for the \code{original}, \code{remaining}, and \code{relocated}
#' populations. This separation enables multiple dispersal models, representing
#' different dispersal vectors, to be run in sequence.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   grid-based spatial region (template) or patch-based network for spread
#'   simulations. The region object contains functionality for calculating path
#'   distances and directions, permeability graphs, and structures to
#'   facilitate two-tier grid-based dispersal.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the spread simulations.
#' @param dispersal_stages Numeric vector of population stages (indices) that
#'   disperse. Default is all stages (when set to \code{NULL}).
#' @param proportion The proportion of the (unstructured or staged) population
#'   that disperses at each time step. This parameter is also used to scale the
#'   the number of dispersal destinations selected when the population is
#'   presence-only and the number of dispersal \code{events} is not defined.
#'   Default is \code{NULL} (producing no dispersal unless the population is
#'   presence-only and \code{events} is defined).
#' @param events The mean number of dispersal events generated via a Poisson
#'   distribution for each location at each time step. A dispersal destination
#'   (location) is selected for each dispersal event. Default is \code{NULL}
#'   (resulting in destinations being selected for each individual within
#'   unstructured or staged populations, or stochastic sampling of destinations
#'   for presence-only populations).
#' @param distance_function A function (or kernel) in the form
#'   \code{function(distances)}, that calculates the (relative) probability of
#'   dispersal for each distance (in m) specified as a numeric vector. Default
#'   is none.
#' @param distance_adjust Logical indication of whether the (relative)
#'   probabilities returned by \code{distance_function} should be distributed
#'   across the (approximate) number of grid cells at each distance. When not
#'   specified (\code{NULL}), the value will resolve to \code{TRUE} for
#'   grid-based regions only.
#' @param direction_function A function (or kernel) in the form
#'   \code{function(directions)}, that calculates the (relative) probability of
#'   dispersal for each direction (0-360 degrees) specified as an integer
#'   vector. Default is none.
#' @param attractors List containing \code{Attractor} (or inherited) class
#'   objects for spatially weighted dispersal, and/or containing a named
#'   numeric multiplier \code{source_density = <numeric>}, which (when present)
#'   is used to dynamically create an "attractor" based on the population
#'   density (number/capacity) at each dispersal location source (origin) at
#'   each simulation time step. Default is empty.
#' @param permeability A \code{Permeability} class (or inherited) class object
#'   for representing spatial permeability or constraints. Default is none.
#' @param max_distance The maximum dispersal distance (in m). Default is
#'   \code{NULL} (resulting in no explicit distance limit).
#' @param class Character class name for inherited classes. Default is
#'   \code{NULL}.
#' @param ... Additional parameters.
#' @return A \code{Dispersal} class object (list) containing functions for
#'   accessing attributes (of the function environment) and performing
#'   dispersal:
#'   \describe{
#'     \item{\code{set_cores(cores)}}{Set the number of cores available for
#'       parallel processing.}
#'     \item{\code{pack(n)}}{Packs a population vector or matrix \code{n} into
#'       a list containing a vector of occupied location (cell or patch)
#'       \code{indices}, the \code{original} population values at the occupied
#'       locations only, the \code{remaining} occupied values (initially a
#'       duplicate of \code{original}), and a vector or matrix for the
#'       \code{relocated} population values at all locations (initially all
#'       zero).}
#'     \item{\code{unpack(n)}}{Unpacks a population list by combining the
#'       \code{remaining} and \code{relocated} population values to form a
#'       new post-dispersal population vector or matrix.}
#'     \item{\code{disperse(n)}}{Perform location dispersal on a list \code{n}
#'       of vectors or matrices, representing the occupied location (cell or
#'       patch) \code{indices}, the \code{original} occupied populations, the
#'       \code{remaining} occupied populations, and the \code{relocated}
#'       populations (at all region locations), and return the transformed list
#'       of vectors or matrices. The separation of original, remaining and
#'       relocated populations enables multiple models for different dispersal
#'       vectors to run in sequence.}
#'   }
#' @include Region.R
#' @export
Dispersal <- function(region, population_model,
                      dispersal_stages = NULL,
                      proportion = NULL,
                      events = NULL,
                      distance_function = NULL,
                      distance_adjust = NULL,
                      direction_function = NULL,
                      attractors = list(),
                      permeability = NULL,
                      max_distance = NULL,
                      class = character(), ...) {
  UseMethod("Dispersal")
}

#' @name Dispersal
#' @export
Dispersal.Region <- function(region, population_model,
                             dispersal_stages = NULL,
                             proportion = NULL,
                             events = NULL,
                             distance_function = NULL,
                             distance_adjust = NULL,
                             direction_function = NULL,
                             attractors = list(),
                             permeability = NULL,
                             max_distance = NULL,
                             class = character(), ...) {

  # Check the population model
  if (!is.null(population_model) &&
      !inherits(population_model, "Population")) {
    stop("Population model must be a 'Population' or inherited class object.",
         call. = FALSE)
  }

  # Variables that depend on population type
  population_type <- population_model$get_type()
  if (population_type == "presence_only") {
    dispersal_stages <- 1
  } else if (population_type == "unstructured") {
    capacity_stages <- 1
    dispersal_stages <- 1
  } else if (population_type == "stage_structured") {
    capacity_stages <- population_model$get_capacity_stages()
    if (is.null(dispersal_stages)) {
      dispersal_stages <- 1:population_model$get_stages()
    } else if (!is.numeric(dispersal_stages) ||
               !all(dispersal_stages %in% 1:population_model$get_stages())) {
      stop("Dispersal stages must be a vector of stage indices consistent ",
           "with the population model.", call. = FALSE)
    }
  }

  # Check proportion, events, distance function, and direction function
  if (!is.null(proportion) && (!is.numeric(proportion) || proportion <= 0)) {
    stop("The proportion parameter must be numeric and > 0.", call. = FALSE)
  }
  if (!is.null(events) && (!is.numeric(events) || events <= 0)) {
    stop("The events parameter must be numeric and > 0.", call. = FALSE)
  }
  if (!is.null(distance_function) && !is.function(distance_function)) {
    stop("The distance function must be a function.", call. = FALSE)
  }
  if (!is.null(distance_adjust) && !is.logical(distance_adjust)) {
    stop("The distance adjust must be logical.", call. = FALSE)
  }
  if (is.null(distance_adjust)) {
    distance_adjust <- region$get_type() == "grid"
  }
  if (!is.null(direction_function) && !is.function(direction_function)) {
    stop("The direction function must be a function.", call. = FALSE)
  }

  # Check the attractors
  if (!is.list(attractors) ||
      length(attractors) && !all(sapply(1:length(attractors), function(i) {
        (inherits(attractors[[i]], "Attractor") ||
         (is.numeric(attractors[[i]]) & attractors[[i]] > 0 &&
          names(attractors)[i] == "source_density"))
      }))) {
    stop("Attractors must be a list containing zero or more 'Attractor' or ",
         "inherited class objects, and/or a numeric attractor (> 0) named ",
         "'source_density'.", call. = FALSE)
  }

  # Check permeability
  if (!is.null(permeability) &&
      !inherits(permeability, "Permeability")) {
    stop("Permeability parameter must be a 'Permeability' or inherited class ",
         "object.", call. = FALSE)
  }

  # Check maximum distance parameter
  if (!is.null(max_distance) &&
      (!is.numeric(max_distance) || max_distance <= 0)) {
    stop("The maximum distance parameter must be numeric and > 0.",
         call. = FALSE)
  }

  # Configure paths via region (sets permeability id when present)
  region$configure_paths(directions = is.function(direction_function),
                         max_distance = max_distance,
                         permeability = permeability)

  # Get the permeability identifier
  perm_id <- NULL
  if (inherits(permeability, "Permeability")) {
    perm_id <- permeability$get_id()
  }

  # Create a class structure
  self <- structure(list(), class = c(class, "Dispersal"))

  # Set the number of cores available for parallel processing
  parallel_cores <- 1
  self$set_cores <- function(cores) {
    parallel_cores <<- cores
  }

  # Pack population vector/matrix into a list for dispersal
  self$pack <- function(n) {
    indices = which(rowSums(as.matrix(n)) > 0)
    n <- list(indices = indices,
              original = as.matrix(n)[indices,, drop = FALSE],
              remaining = as.matrix(n)[indices,, drop = FALSE],
              relocated = n)
    n$relocated[] <- FALSE
    return(n)
  }

  # Unpack post-dispersal population vector/matrix from packed list
  self$unpack <- function(n) {
    if (population_type == "presence_only") {
      n$relocated[n$indices] <- n$remaining
    } else if (population_type == "unstructured") {
      n$relocated[n$indices] <- n$relocated[n$indices] + n$remaining
    } else if (population_type == "stage_structured") {
      n$relocated[n$indices,] <- n$relocated[n$indices,] + n$remaining
    }
    return(n$relocated)
  }

  # Generic disperse method (may be overridden in subclass)
  self$disperse <- function(n) {

    # Ensure cell characters are stored without scientific notation
    options(scipen = 999)

    # Calculate paths
    region$calculate_paths(n$indices)

    # Limit to dispersal-ready when presence-only with spread delay
    if (population_type == "presence_only" &&
        "spread_delay" %in% names(attributes(n$relocated))) {
      dispersal_ready <-
        which(attr(n$relocated, "spread_delay")[n$indices] == 0)
    } else {
      dispersal_ready <- 1:length(n$indices)
    }

    # Perform dispersal for each (dispersal-ready) occupied location
    doParallel::registerDoParallel(cores = min(parallel_cores,
                                               length(dispersal_ready)))
    dispersal_list <- foreach(
      i = dispersal_ready,
      .errorhandling = c("stop"),
      .packages = c(),
      .export = c(),
      .noexport = c()) %dopar% {

        ## Map path lists use location indices as characters
        loc_i <- n$indices[i]
        loc_char <- as.character(loc_i)

        ## Get paths for location
        paths <- region$get_paths(loc_i,
                                  directions = is.function(direction_function),
                                  max_distance = max_distance,
                                  perm_id = perm_id)
        if (region$get_type() == "patch") { # for code reuse
          paths <- lapply(paths, function(p)
            lapply(p, function(l) list(cell = l)))
        }

        ## Relative dispersal probabilities for source and destination
        source_p <- 1
        destination_p <- list(cell = rep(1,
                                         length(paths$idx[[loc_char]]$cell)))
        if (region$two_tier()) { # multiple cells
          destination_p$aggr <- region$get_aggr()$rast[
            region$get_aggr()$indices[paths$idx[[loc_char]]$aggr]][,1]
        }

        ## Adjust destination (relative) probabilities based on distance
        if (is.function(distance_function)) {

          # Additional adjustment given (approx.) cells at each distance
          dist_p_adj <- list(cell = 1, aggr = 1)
          if (distance_adjust) {
            dist_p_adj$cell <- (region$get_res()/
                                  (2*pi*paths$distances[[loc_char]]$cell))
            if (region$two_tier()) {
              dist_p_adj$aggr <- (region$get_res()/
                                    (2*pi*paths$distances[[loc_char]]$aggr))
            }
          }

          # Adjust via distance function using appropriate distances
          if (is.numeric(perm_id) && "perm_dist" %in% names(paths)) {
            destination_p$cell <-
              (distance_function(paths$perm_dist[[loc_char]]$cell)*
                 dist_p_adj$cell*destination_p$cell)
            if (region$two_tier()) {
              destination_p$aggr <-
                (distance_function(paths$perm_dist[[loc_char]]$aggr)*
                   dist_p_adj$aggr*destination_p$aggr)
            }
          } else {
            destination_p$cell <-
              (distance_function(paths$distances[[loc_char]]$cell)*
                 dist_p_adj$cell*destination_p$cell)
            if (region$two_tier()) {
              destination_p$aggr <-
                (distance_function(paths$distances[[loc_char]]$aggr)*
                   dist_p_adj$aggr*destination_p$aggr)
            }
          }
          rm(dist_p_adj)
        }

        ## Adjust destination (relative) probabilities based on direction
        if (is.function(direction_function)) {
          destination_p$cell <-
            (destination_p$cell*
               direction_function(paths$directions[[loc_char]]$cell))
          if (region$two_tier()) {
            destination_p$aggr <-
              (destination_p$aggr*
                 direction_function(paths$directions[[loc_char]]$aggr))
          }
        }

        ## Adjust (relative) probabilities via attractors
        if (length(attractors)) {

          # Source density attractor
          if ("source_density" %in% names(attractors)) {
            cell_k <- population_model$get_capacity(loc_i)
            if (is.numeric(cell_k)) {
              if (cell_k > 0) {
                source_p <- min((attractors$source_density*source_p*
                                   sum(n$original[i, capacity_stages])/cell_k),
                                1)
              } else {
                source_p <- 0
              }
            }
          }

          # Other attractors
          for (attractor in attractors) {
            if (inherits(attractor, "Attractor")) {

              # Source attractors
              if (attractor$get_type() %in% c("source", "both")) {
                source_p <- (source_p*attractor$get_values(loc_i))
              }

              # Destination attractors
              if (attractor$get_type() %in% c("destination", "both")) {
                destination_p$cell <-
                  (destination_p$cell*
                     attractor$get_values(paths$idx[[loc_char]]$cell))
                if (region$two_tier()) {
                  destination_p$aggr <-
                    (destination_p$aggr*
                       attractor$get_aggr_values(paths$idx[[loc_char]]$aggr))
                }
              }
            }
          }
        }

        ## Calculate dispersers from (original) cell population
        dispersers <- FALSE
        if (population_type == "presence_only") {
          dispersers <- TRUE
        } else if (is.numeric(proportion) && proportion  > 0) {

          # Generate dispersers
          dispersers <- stats::rbinom(length(dispersal_stages),
                                      size = n$original[i, dispersal_stages],
                                      prob = proportion*source_p)

          # Ensure there are sufficient remaining (after other dispersals)
          dispersers <- pmin(dispersers, n$remaining[i, dispersal_stages])
        }

        ## Distribute dispersers across dispersal events
        dispersals <- FALSE
        if (sum(dispersers) > 0) {

          if (is.numeric(events)) {

            # Generate the number of events via Poisson distribution
            dispersals <- stats::rpois(1, events*source_p)

            # Distribute dispersers
            if (dispersals > 0) {

              if (population_type == "presence_only") {

                # Distribute logical presence
                dispersers <- as.matrix(rep(TRUE, dispersals))

              } else {

                # Distribute dispersers via a multinomial generation
                dispersers <- t(array(dispersers,
                                      c(length(dispersers), dispersals)))
                for (stage in dispersal_stages) {
                  dispersers[, stage] <-
                    stats::rmultinom(1, size = dispersers[1, stage],
                                     prob = rep(1, dispersals))
                }

                # Remove any empty events
                dispersers <- dispersers[which(rowSums(dispersers) > 0),,
                                         drop = FALSE]
                dispersals <- nrow(dispersers)
              }
            }

          } else {

            if (population_type == "presence_only") {

              # Probabilistic dispersal to each reachable destination
              if (is.numeric(proportion) && proportion  > 0) {
                dispersals <- TRUE
              }

            } else if (population_type == "unstructured") {

              # Dispersal event per individual
              dispersals <- dispersers
              dispersers <- as.matrix(rep(1, dispersals))

            } else if (population_type == "stage_structured") {

              # Dispersal event per individual
              dispersals <- sum(dispersers)
              dispersers <- +(
                array(unlist(sapply(1:length(dispersers),
                                    function(j) rep(j, dispersers[j]))),
                      c(dispersals, length(dispersers))) ==
                  matrix(rep(1:length(dispersers), each = dispersals),
                         nrow = dispersals))
            }
          }
        }

        ## Sample dispersal destinations
        if (dispersals) {

          if (population_type == "presence_only" && is.null(events)) {

            # Probabilistic dispersal to each reachable destination
            destinations <- which(as.logical(stats::rbinom(
              length(destination_p$cell), size = 1,
              prob = destination_p$cell*proportion*source_p)))

            if (region$two_tier()) {

              # Number or (non-NA) region cells in each aggregate cell
              aggr_cells <- region$get_aggr()$rast[
                region$get_aggr()$indices[paths$idx[[loc_char]]$aggr]][,1]

              # Sample region cell counts for aggregate cells
              aggr_cells <- stats::rbinom(
                length(destination_p$aggr), size = aggr_cells,
                prob = destination_p$aggr*proportion*source_p/aggr_cells)

              # Replicate aggregate destination cells using cell counts
              aggr_dest <- which(as.logical(aggr_cells))
              aggr_dest <- unlist(sapply(aggr_dest,
                                         function(d) rep(d, aggr_cells[d])))

              # Append aggregate destinations
              destinations <- c(destinations,
                                aggr_dest + length(destination_p$cell))
            }

          } else {

            # Sample specified number of dispersal event destinations
            destinations <- sample(1:length(unlist(destination_p)),
                                   size = dispersals, replace = TRUE,
                                   prob = unlist(destination_p))
          }
        }

        ## Resolve dispersal destinations as region cell indices
        if (dispersals) {
          if (region$two_tier()) {

            # Identify cell and aggregate cell destinations
            cell_dest <- which(destinations <= length(destination_p$cell))
            aggr_dest <- which(destinations > length(destination_p$cell))

            # Map cell destinations
            destinations[cell_dest] <-
              paths$idx[[loc_char]]$cell[destinations[cell_dest]]

            # Map aggregate destinations
            destinations[aggr_dest] <-
              paths$idx[[loc_char]]$aggr[destinations[aggr_dest] -
                                           length(destination_p$cell)]

            # Sample region cell(s) within each aggregate destination cell
            for (ad_i in aggr_dest) {

              # Get region cell raster indices
              aggr_cells <- region$get_aggr()$get_cells(destinations[ad_i])

              # Perform a weighted sample via attractors when present
              aggr_p <- rep(1, length(aggr_cells))
              for (attractor in attractors) {
                if (inherits(attractor, "Attractor")) {
                  if (attractor$get_type() %in% c("destination", "both")) {
                    aggr_p <- aggr_p*attractor$get_values(aggr_cells)
                  }
                }
              }
              aggr_i <- aggr_cells[sample(1:length(aggr_cells), size = 1,
                                          prob = aggr_p)]

              # Substitute region cell for aggregate destination
              destinations[ad_i] <- aggr_i
            }

          } else {
            destinations <- paths$idx[[loc_char]]$cell[destinations]
          }
        }

        ## Perform dispersal to cell destinations
        if (dispersals) {

          # Calculate remaining population
          if (population_type %in% c("unstructured", "stage_structured")) {
            remaining <- n$remaining[i, dispersal_stages] - colSums(dispersers)
          } else {
            remaining <- n$remaining[i, dispersal_stages]
          }

          # Apply establishment survival to dispersers (some deaths)
          establish_p <- population_model$get_establish_pr(destinations)
          if (population_type == "presence_only" && is.null(events)) {
            destinations <- destinations[which(as.logical(
              stats::rbinom(length(destinations), size = 1,
                            prob = establish_p)))]
          } else {
            for (stage in dispersal_stages) {
              dispersers[, stage] <- stats::rbinom(dispersals,
                                                   size = dispersers[, stage],
                                                   prob = establish_p)
            }
            destinations <- destinations[which(rowSums(dispersers) > 0)]
            dispersers <- dispersers[which(rowSums(dispersers) > 0),,
                                     drop = FALSE]
            dispersals <- nrow(dispersers)
          }

          # Return population relocation components
          list(i = i, dispersals = dispersals, remaining = remaining,
               destinations = destinations, dispersers = dispersers)

        } else {

          # Return empty dispersals
          list(i = i, dispersals = dispersals)
        }
      }
    doParallel::stopImplicitCluster()

    # Perform dispersal to cell destinations
    for (d in dispersal_list) {

      if (d$dispersals) {

        # Update remaining population
        n$remaining[d$i, dispersal_stages] <- d$remaining

        # Update relocated population
        if (population_type == "presence_only") {
          n$relocated[d$destinations] <- TRUE
        } else if (population_type == "unstructured") {
          n$relocated[d$destinations] <-
            n$relocated[d$destinations] + d$dispersers
        } else if (population_type == "stage_structured") {
          n$relocated[d$destinations, dispersal_stages] <-
            n$relocated[d$destinations, dispersal_stages] + d$dispersers
        }
      }
    }

    # Reinstate default scientific notation
    options(scipen = 0)

    return(n)
  }

  return(self)
}
