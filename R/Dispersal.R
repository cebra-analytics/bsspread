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
#'   that disperses from each occupied location at each time step. It may be
#'   specified as a single numeric value or, when applicable, with spatial
#'   and/or temporal variation via a matrix of spatial (rows) and/or temporal
#'   (columns). Spatial values should be specified via a row for each location,
#'   else a single row may specify temporal variation only. Likewise, a single
#'   column may specify spatial variation only. The number of columns for
#'   temporal variation should either coincide with the number of simulation
#'   time steps, or be a cyclic pattern (e.g. 12 columns for seasonal variation
#'   with monthly time steps). This parameter may also be used to scale the
#'   number of dispersal destinations selected when the population is
#'   presence-only and the number of dispersal \code{events} is not defined.
#'   Default is \code{NULL} (producing no dispersal unless the population is
#'   presence-only and \code{events} is defined).
#' @param events The mean number of dispersal events generated via a Poisson
#'   distribution for each location at each time step. It may be specified as
#'   a single numeric value or, when applicable, with spatial and/or temporal
#'   variation via a matrix of spatial (rows) and/or temporal (columns).
#'   Spatial values should be specified via a row for each location, else a
#'   single row may specify temporal variation only. Likewise, a single column
#'   may specify spatial variation only. The number of columns for temporal
#'   variation should either coincide with the number of simulation time steps,
#'   or be a cyclic pattern (e.g. 12 columns for seasonal variation with
#'   monthly time steps). A dispersal destination (location) is selected for
#'   each dispersal event. Default is \code{NULL} (resulting in destinations
#'   being selected for each individual within unstructured or staged
#'   populations, or stochastic sampling of destinations for presence-only
#'   populations).
#' @param density_dependent Logical to indicate that dispersal is density
#'   dependent, whereby the proportion dispersing and/or the number of
#'   dispersal events generated is scaled by the (unstructured or staged)
#'   population density (number/capacity) at each occupied location at each
#'   simulation time step. Default is \code{FALSE} for no density dependence.
#' @param distance_function A function (or kernel) in the form
#'   \code{function(distances)}, that calculates the (relative) probability of
#'   dispersal for each distance (in m) specified as a numeric vector. Default
#'   is none.
#' @param direction_function A function (or kernel) in the form
#'   \code{function(directions)}, that calculates the (relative) probability of
#'   dispersal for each direction (0-360 degrees) specified as an integer
#'   vector. Default is none.
#' @param combined_function A function (or kernel) in the form
#'   \code{function(distances, directions)}, that calculates the (relative)
#'   probability of dispersal for distances, a list of equal-length numeric
#'   vectors for (1) actual and (2) optional permeable distances (in m), and
#'   a corresponding equal-length numeric vector of directions (0-360 degrees).
#'   Default is none.
#' @param distance_adjust Logical indication of whether the (relative)
#'   probabilities returned by \code{distance_function} or
#'   \code{combined_function} should be distributed across the (approximate)
#'   number of grid cells at each distance. When not specified (\code{NULL}),
#'   the value will resolve to \code{TRUE} for grid-based regions only.
#' @param attractors List containing \code{Attractor} (or inherited) class
#'   objects for spatially weighted dispersal to destination locations. Default
#'   is empty.
#' @param permeability A \code{Permeability} class (or inherited) class object
#'   for representing spatial permeability or constraints. Default is none.
#' @param max_distance The maximum dispersal distance (in m) in each time step.
#'   Default is \code{NULL} (resulting in no explicit distance limit).
#' @param class Character class name for inherited classes. Default is
#'   \code{NULL}.
#' @param ... Additional parameters.
#' @return A \code{Dispersal} class object (list) containing functions for
#'   accessing attributes (of the function environment) and performing
#'   dispersal:
#'   \describe{
#'     \item{\code{set_cores(cores)}}{Set the number of cores available for
#'       parallel processing and thus enable parallel processing for
#'       calculating path distances and directions.}
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
#'     \item{\code{disperse(n, tm)}}{Perform location dispersal at simulation
#'       time step \code{tm} on a list \code{n} of vectors or matrices,
#'       representing the occupied location (cell or patch) \code{indices},
#'       the \code{original} occupied populations, the \code{remaining}
#'       occupied populations, and the \code{relocated} populations (at all
#'       region locations), and return the transformed list of vectors or
#'       matrices. The separation of original, remaining and relocated
#'       populations enables multiple models for different dispersal vectors
#'       to run in sequence. Spread and establishment control (suppression) may
#'       also be processed when passed via attributes (see
#'       \code{bsmanage::ManageControls}).}
#'   }
#' @references
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
#'   Robinet, C., Kehlenbeck, H., Kriticos, D. J., Baker, R. H. A.,
#'   Battisti, A., Brunel, S., Dupin, M., Eyre, D., Faccoli, M., Ilieva, Z.,
#'   Kenis, M., Knight, J., Reynaud, P., Yart, A., & van der Werf, W. (2012).
#'   A Suite of Models to Support the Quantitative Assessment of Spread in Pest
#'   Risk Analysis. \emph{PLoS ONE}, 7(10), 1–18.
#'   \doi{10.1371/journal.pone.0043366}
#' @include Region.R
#' @export
Dispersal <- function(region, population_model,
                      dispersal_stages = NULL,
                      proportion = NULL,
                      events = NULL,
                      density_dependent = FALSE,
                      distance_function = NULL,
                      direction_function = NULL,
                      combined_function = NULL,
                      distance_adjust = NULL,
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
                             density_dependent = FALSE,
                             distance_function = NULL,
                             direction_function = NULL,
                             combined_function = NULL,
                             distance_adjust = NULL,
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

  # Check proportion, events, density_dependent, distance function, and
  # direction function
  if (!is.null(proportion)) {
    if ((!is.matrix(proportion) &&
         !(length(proportion) %in% c(1, region$get_locations()))) ||
        (is.matrix(proportion) &&
         !(nrow(proportion) %in% c(1, region$get_locations())))) {
      stop(paste("The proportion parameter should be a single value or",
                 "a matrix with a single row or a row for each region",
                 "location."), call. = FALSE)
    } else if (!is.numeric(proportion) ||  any(proportion < 0)) {
      stop("The proportion parameter must be numeric >= 0.", call. = FALSE)
    }
    proportion <- as.matrix(proportion)
  }
  if (!is.null(events)) {
    if ((!is.matrix(events) &&
         !(length(events) %in% c(1, region$get_locations()))) ||
        (is.matrix(events) &&
         !(nrow(events) %in% c(1, region$get_locations())))) {
      stop(paste("The events parameter should be a single value or",
                 "a matrix with a single row or a row for each region",
                 "location."), call. = FALSE)
    } else if (!is.numeric(events) ||  any(events < 0)) {
      stop("The events parameter must be numeric >= 0.", call. = FALSE)
    }
    events <- as.matrix(events)
  }
  if (!is.null(density_dependent) && !is.logical(density_dependent)) {
    stop("The density dependent parameter must be logical.", call. = FALSE)
  }
  if (!is.null(distance_function) && !is.function(distance_function)) {
    stop("The distance function must be a function.", call. = FALSE)
  }
  if (!is.null(direction_function) && !is.function(direction_function)) {
    stop("The direction function must be a function.", call. = FALSE)
  }
  if (!is.null(combined_function) && !is.function(combined_function)) {
    stop("The combined function must be a function.", call. = FALSE)
  }
  if (!is.null(distance_adjust) && !is.logical(distance_adjust)) {
    stop("The distance adjust parameter must be logical.", call. = FALSE)
  }
  if (is.null(distance_adjust)) {
    distance_adjust <- region$get_type() == "grid"
  }

  # Check the attractors
  if (!is.list(attractors) ||
      length(attractors) && !all(sapply(1:length(attractors), function(i) {
        (inherits(attractors[[i]], "Attractor"))
      }))) {
    stop("Attractors must be a list containing zero or more 'Attractor' or ",
         "inherited class objects.", call. = FALSE)
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
  region$configure_paths(directions = (is.function(direction_function) ||
                                         is.function(combined_function)),
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
  parallel_cores <- NULL
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
  self$disperse <- function(n, tm) {

    # Ensure cell characters are stored without scientific notation
    options(scipen = 999)

    # Calculate paths
    region$calculate_paths(n$indices)

    # Limit to dispersal-ready when presence-only with spread delay
    if (population_type == "presence_only" &&
        "spread_delay" %in% names(attributes(n$relocated))) {
      dispersal_ready <-
        which(attr(n$relocated, "spread_delay")[n$indices] == 0)
    } else if (length(n$indices)){
      dispersal_ready <- 1:length(n$indices)
    } else {
      dispersal_ready <- n$indices
    }

    # Function to calculate dispersals from a single occupied location (index)
    calculate_dispersals <- function(i) {

      ## Map path lists use location indices as characters
      loc_i <- n$indices[i]
      loc_char <- as.character(loc_i)

      ## Get paths for location
      paths <- region$get_paths(
        loc_i,
        directions = (is.function(direction_function) ||
                        is.function(combined_function)),
        max_distance = max_distance,
        perm_id = perm_id)
      if (region$get_type() == "patch") { # for code reuse
        paths <- lapply(paths, function(p)
          lapply(p, function(l) list(cell = l)))
      }

      ## Check that dispersal paths are present
      if (length(unlist(paths$idx[[loc_char]])) == 0) {
        return(list(i = i, dispersals = FALSE))
      }
      aggr_paths_present <- region$two_tier()
      if (length(paths$idx[[loc_char]]$aggr) == 0) {
        aggr_paths_present <- FALSE
      }

      ## Relative dispersal probabilities for destinations
      destination_p <- list(cell = rep(1,
                                       length(paths$idx[[loc_char]]$cell)))
      if (aggr_paths_present) {
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
          if (aggr_paths_present) {
            dist_p_adj$aggr <- (region$get_res()/
                                  (2*pi*paths$distances[[loc_char]]$aggr))
          }
        }

        # Adjust via distance function using appropriate distances
        if (is.numeric(perm_id) && "perm_dist" %in% names(paths)) {
          destination_p$cell <-
            (distance_function(paths$perm_dist[[loc_char]]$cell)*
               dist_p_adj$cell*destination_p$cell)
          if (aggr_paths_present) {
            destination_p$aggr <-
              (distance_function(paths$perm_dist[[loc_char]]$aggr)*
                 dist_p_adj$aggr*destination_p$aggr)
          }
        } else {
          destination_p$cell <-
            (distance_function(paths$distances[[loc_char]]$cell)*
               dist_p_adj$cell*destination_p$cell)
          if (aggr_paths_present) {
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
        if (aggr_paths_present) {
          destination_p$aggr <-
            (destination_p$aggr*
               direction_function(paths$directions[[loc_char]]$aggr))
        }
      }

      ## Adjust destination (relative) probabilities based on combined distance
      ## and direction
      if (is.function(combined_function)) {

        # Additional adjustment given (approx.) cells at each distance
        dist_p_adj <- list(cell = 1, aggr = 1)
        if (distance_adjust) {
          dist_p_adj$cell <- (region$get_res()/
                                (2*pi*paths$distances[[loc_char]]$cell))
          if (aggr_paths_present) {
            dist_p_adj$aggr <- (region$get_res()/
                                  (2*pi*paths$distances[[loc_char]]$aggr))
          }
        }

        # Adjust via combined function using appropriate distances/directions
        if (is.numeric(perm_id) && "perm_dist" %in% names(paths)) {
          destination_p$cell <-
            (combined_function(list(paths$distances[[loc_char]]$cell,
                                    paths$perm_dist[[loc_char]]$cell),
                               paths$directions[[loc_char]]$cell)*
               dist_p_adj$cell*destination_p$cell)
          if (aggr_paths_present) {
            destination_p$aggr <-
              (combined_function(list(paths$distances[[loc_char]]$aggr,
                                      paths$perm_dist[[loc_char]]$aggr),
                                 paths$directions[[loc_char]]$aggr)*
                 dist_p_adj$aggr*destination_p$aggr)
          }
        } else {
          destination_p$cell <-
            (combined_function(list(paths$distances[[loc_char]]$cell),
                               paths$directions[[loc_char]]$cell)*
               dist_p_adj$cell*destination_p$cell)
          if (aggr_paths_present) {
            destination_p$aggr <-
              (combined_function(list(paths$distances[[loc_char]]$aggr),
                                 paths$directions[[loc_char]]$aggr)*
                 dist_p_adj$aggr*destination_p$aggr)
          }
        }
        rm(dist_p_adj)
      }

      ## Adjust (relative) probabilities via attractors
      if (length(attractors)) {
        for (attractor in attractors) {
          if (inherits(attractor, "Attractor")) {
            destination_p$cell <-
              (destination_p$cell*
                 attractor$get_values(paths$idx[[loc_char]]$cell))
            if (aggr_paths_present) {
              destination_p$aggr <-
                (destination_p$aggr*
                   attractor$get_aggr_values(paths$idx[[loc_char]]$aggr))
            }
          }
        }
      }

      ## Apply density dependence
      source_p <- 1
      if (density_dependent) {
        cell_k <- population_model$get_capacity(cells = loc_i, tm = tm)
        if (is.numeric(cell_k)) {
          if (cell_k > 0) {
            source_p <- min(sum(n$original[i, capacity_stages])/cell_k, 1)
          } else {
            source_p <- 0
          }
        }
      }

      ## Extract proportion and/or events for location / time step
      if (is.numeric(proportion)) {
        if (!is.numeric(tm) || (is.numeric(tm) && tm == 0)) {
          tm_i = 1
        } else { # wrap
          tm_i <- ((tm + (ncol(proportion) - 1)) %% ncol(proportion)) + 1
        }
        if (nrow(proportion) == region$get_locations()) {
          proportion_i <- proportion[loc_i, tm_i]
        } else {
          proportion_i <- proportion[, tm_i]
        }
      } else {
        proportion_i <- proportion
      }
      if (is.numeric(events)) {
        if (!is.numeric(tm) || (is.numeric(tm) && tm == 0)) {
          tm_i = 1
        } else { # wrap
          tm_i <- ((tm + (ncol(events) - 1)) %% ncol(events)) + 1
        }
        if (nrow(events) == region$get_locations()) {
          events_i <- events[loc_i, tm_i]
        } else {
          events_i <- events[, tm_i]
        }
      } else {
        events_i <- events
      }

      ## Process spread control/suppression
      if (!is.null(attr(n$relocated, "control_spread"))) {
        if (length(attr(n$relocated, "control_spread")) == 1) {
          control_i <- attr(n$relocated, "control_spread")
        } else {
          control_i <- attr(n$relocated, "control_spread")[loc_i]
        }
        if (is.numeric(proportion_i) && proportion_i  > 0) {
          proportion_i <- proportion_i*control_i
        }
        if (is.numeric(events_i) && events_i > 0) {
          events_i <- events_i*control_i
        }
      }

      ## Calculate dispersers from (original) cell population
      dispersers <- FALSE
      if (population_type == "presence_only") {
        dispersers <- TRUE
      } else if (is.numeric(proportion_i) && proportion_i  > 0) {

        # Generate dispersers
        dispersers <- stats::rbinom(length(dispersal_stages),
                                    size = n$original[i, dispersal_stages],
                                    prob = pmin(proportion_i*source_p, 1))

        # Ensure there are sufficient remaining (after other dispersals)
        dispersers <- pmin(dispersers, n$remaining[i, dispersal_stages])
      }

      ## Distribute dispersers across dispersal events
      dispersals <- FALSE
      if (sum(dispersers) > 0) {

        if (is.numeric(events_i) && events_i > 0) {

          # Generate the number of events via Poisson distribution
          dispersals <- stats::rpois(1, events_i*source_p)

          # Distribute dispersers
          if (dispersals > 0) {

            if (population_type == "presence_only") {

              # Distribute logical presence
              dispersers <- as.matrix(rep(TRUE, dispersals))

            } else {

              # Distribute dispersers via a multinomial generation
              dispersers <- t(array(dispersers,
                                    c(length(dispersers), dispersals)))
              for (ds_i in 1:length(dispersal_stages)) {
                dispersers[, ds_i] <-
                  stats::rmultinom(1, size = dispersers[1, ds_i],
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
            if (is.numeric(proportion_i) && proportion_i  > 0) {
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
            prob = pmin(destination_p$cell*proportion_i, 1))))

          if (aggr_paths_present) {

            # Number or (non-NA) region cells in each aggregate cell
            aggr_cells <- region$get_aggr()$rast[
              region$get_aggr()$indices[paths$idx[[loc_char]]$aggr]][,1]

            # Sample region cell counts for aggregate cells
            aggr_cells <- stats::rbinom(
              length(destination_p$aggr), size = aggr_cells,
              prob = pmin(destination_p$aggr*proportion_i/aggr_cells, 1))

            # Replicate aggregate destination cells using cell counts
            aggr_dest <- which(as.logical(aggr_cells))
            aggr_dest <- unlist(sapply(aggr_dest,
                                       function(d) rep(d, aggr_cells[d])))

            # Append aggregate destinations
            destinations <- c(destinations,
                              aggr_dest + length(destination_p$cell))
          }

          # Any destinations?
          if (length(destinations) == 0) {
            dispersals <- FALSE
          }

        } else {

          # Sample specified number of dispersal event destinations
          if (any(unlist(destination_p) > 0)) {
            destinations <- sample(1:length(unlist(destination_p)),
                                   size = dispersals, replace = TRUE,
                                   prob = unlist(destination_p))
          } else {
            dispersals <- FALSE
          }
        }
      }

      ## Resolve dispersal destinations as region cell indices
      if (dispersals) {
        if (aggr_paths_present) {

          # Identify cell and aggregate cell destinations
          cell_dest <- which(destinations <= length(destination_p$cell))
          aggr_dest <- which(destinations > length(destination_p$cell))

          # Map cell destinations
          destinations[cell_dest] <-
            paths$idx[[loc_char]]$cell[destinations[cell_dest]]

          if (length(aggr_dest)) {

            # Map aggregate destinations
            destinations[aggr_dest] <-
              paths$idx[[loc_char]]$aggr[destinations[aggr_dest] -
                                           length(destination_p$cell)]

            # Collect consecutive repeated values in a list
            aggr_dest_rep <- list(destinations[aggr_dest][1])
            for (ag_i in destinations[aggr_dest][-1]) {
              if (ag_i ==
                  unlist(aggr_dest_rep)[length(unlist(aggr_dest_rep))]) {
                aggr_dest_rep[[length(aggr_dest_rep)]] <-
                  c(aggr_dest_rep[[length(aggr_dest_rep)]], ag_i)
              } else {
                aggr_dest_rep <- c(aggr_dest_rep, list(ag_i))
              }
            }

            # Sample region cell(s) within each aggregate destination cell
            aggr_sample_replace <- !(population_type == "presence_only" &&
                                       is.null(events))
            for (rep_i in 1:length(aggr_dest_rep)) {

              # Get region cell raster indices
              aggr_cells <-
                region$get_aggr()$get_cells(aggr_dest_rep[[rep_i]][1])

              # Perform a weighted sample via attractors when present
              aggr_p <- rep(1, length(aggr_cells))
              for (attractor in attractors) {
                if (inherits(attractor, "Attractor")) {
                  aggr_p <- aggr_p*attractor$get_values(aggr_cells)
                }
              }
              aggr_i <- aggr_cells[sample(
                1:length(aggr_cells),
                size = length(aggr_dest_rep[[rep_i]]),
                replace = aggr_sample_replace,
                prob = aggr_p)]

              # Substitute region cells for repeated aggregate destinations
              aggr_dest_rep[[rep_i]] <- aggr_i
            }

            # Substitute region cells for all aggregate destinations
            destinations[aggr_dest] <- unlist(aggr_dest_rep)
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

        # Get establishment probability for destination cells at time step
        establish_p <- population_model$get_establish_pr(cells = destinations,
                                                         tm = tm)

        # Process establishment control/suppression
        if (!is.null(attr(n$relocated, "control_establishment")) &&
            is.numeric(establish_p)) {
          if (length(attr(n$relocated, "control_establishment")) == 1) {
            establish_p <-
              establish_p*attr(n$relocated, "control_establishment")
          } else {
            establish_p <-
              establish_p*attr(n$relocated,
                               "control_establishment")[destinations]
          }
        }

        # Apply establishment survival to dispersers (some deaths)
        if (is.numeric(establish_p)) {
          if (population_type == "presence_only" && is.null(events)) {
            destinations <- destinations[which(as.logical(
              stats::rbinom(length(destinations), size = 1,
                            prob = establish_p)))]
          } else {
            for (sd_i in 1:length(dispersal_stages)) {
              dispersers[, sd_i] <- stats::rbinom(dispersals,
                                                   size = dispersers[, sd_i],
                                                   prob = establish_p)
            }
            destinations <- destinations[which(rowSums(dispersers) > 0)]
            dispersers <- dispersers[which(rowSums(dispersers) > 0),,
                                     drop = FALSE]
            dispersals <- nrow(dispersers)
          }
        }

        # Return population relocation components
        list(i = i, dispersals = dispersals, remaining = remaining,
             destinations = destinations, dispersers = dispersers)

      } else {

        # Return empty dispersals
        list(i = i, dispersals = FALSE)
      }
    }

    # Calculate dispersal for each (dispersal-ready) occupied location
    if (is.numeric(parallel_cores) &&
        min(parallel_cores, length(dispersal_ready)) > 1) {

      # Calculate and collect in parallel
      doParallel::registerDoParallel(cores = min(parallel_cores,
                                                 length(dispersal_ready)))
      dispersal_list <- foreach(
        i = dispersal_ready,
        .errorhandling = c("stop")) %dopar% {
          calculate_dispersals(i)
        }
      doParallel::stopImplicitCluster()

      # Recover from parallel memory failures via serial calculations
      if (any(sapply(dispersal_list, is.null))) {
        message("Parallel dispersal memory failures detected - trying serial")
        for (i in which(sapply(dispersal_list, is.null))) {
          dispersal_list[[i]] <- calculate_dispersals(dispersal_ready[i])
        }
      }

    } else {

      # Calculate and collect in serial
      dispersal_list <- list()
      for (i in dispersal_ready) {
        dispersal_list[[length(dispersal_list) + 1]] <- calculate_dispersals(i)
      }
    }

    # Perform dispersal to cell destinations
    for (d in dispersal_list) {

      if (d$dispersals) {

        # Update remaining population
        n$remaining[d$i, dispersal_stages] <- d$remaining

        # Update relocated population
        if (population_type == "presence_only") {
          n$relocated[unique(d$destinations)] <- TRUE
        } else if (population_type == "unstructured") {
          destinations <- unique(d$destinations)
          dispersers <- sapply(destinations, function(di)
            sum(d$dispersers[which(d$destinations == di)]))
          n$relocated[destinations] <-
            n$relocated[destinations] + dispersers
        } else if (population_type == "stage_structured") {
          destinations <- unique(d$destinations)
          dispersers <- sapply(destinations, function(di)
            colSums(d$dispersers[which(d$destinations == di),,drop = FALSE]))
          if (is.matrix(dispersers)) {
            dispersers <- t(dispersers)
          }
          n$relocated[destinations, dispersal_stages] <-
            n$relocated[destinations, dispersal_stages] + dispersers
        }
      }
    }

    # Reinstate default scientific notation
    options(scipen = 0)

    return(n)
  }

  return(self)
}
