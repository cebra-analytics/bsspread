#' Results class builder
#'
#' Builds a class for encapsulating, calculating, and collating population
#' spread modelling simulation results, including the population at each
#' location at each collation time step, the total population size and area
#' occupied at each time step, as well as spatio-temporal impacts and
#' management action quantities applied (e.g. removals), plus costs when
#' configured. When simulations are replicated, summary results (means and
#' standard deviations) are produced.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the spread simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation and growth functionality for the
#'   spread simulations.
#' @param impacts A list of \code{Impacts} class objects specifying various
#'   impacts of the simulated population.
#' @param actions A list of \code{Actions} or inherited class objects
#'   specifying simulated management actions, such as detection, control, and
#'   removal.
#' @param time_steps The number of discrete time steps to simulate. Default is
#'   1.
#' @param step_duration The duration of the simulation time steps in units
#'   specified by \code{step_units}. Default is 1.
#' @param step_units The units for the simulation step duration
#'   (\code{step_duration}) as a character string. Default is "years".
#' @param collation_steps The interval in time steps for collating results.
#'   Default is 1, that is, results are collated at every time step.
#' @param replicates The number of replicate or repeated simulations to be run.
#'   Default is 1. Note that replicate simulations results are collated as
#'   summary statistics across simulations.
#' @param combine_stages Optionally combine (sum) specified stages (a vector of
#'   stage indices) of stage-based population results. The default
#'   (\code{NULL}) maintains results for each stage.
#' @param class Character class name for inherited classes. Default is
#'   \code{NULL}.
#' @param ... Additional parameters.
#' @return A \code{Results} class object (list) containing functions for
#'   calculating and collating results, as well as accessing lists of results
#'   and the simulation parameters used to produce them:
#'   \describe{
#'     \item{\code{collate(r, tm, n)}}{Collate the results at simulation
#'       replicate \code{r} and time step \code{tm} using the current vector
#'       or array \code{n} representing the population at each location, as
#'       well as any attachments relating to impacts or actions.}
#'     \item{\code{finalize()}}{Finalize the results collation (summary
#'       calculations).}
#'     \item{\code{get_list()}}{Return the results as a list (collated
#'       populations and/or occupancy, totals, area occupied, impacts, and
#'       actions).}
#'     \item{\code{get_params()}}{Get the simulation parameters used.}
#'     \item{\code{save_rasters(...)}}{Save the collated results as raster TIF
#'       files when the region is grid-based. \code{Terra} raster write options
#'       may be passed to the function. Returns a list of saved \code{Terra}
#'       raster layers with attached \code{"metadata"} attribute list, which
#'       includes: plot category; impact or action type & name; impact
#'       incursion type; population stage (when applicable); an indication of
#'       whether or not layers represent costs and/or are cumulative outputs;
#'       multi-replicate summary (when applicable); a descriptive label; units;
#'       scale type; and an indication of whether or not each layer contains
#'       non-zero values.}
#'     \item{\code{save_csv()}}{Save the collated results as comma-separated
#'       values (CSV) files when the region is patch-based. Also saves the
#'       population totals, area occupied, as well as impact and action totals
#'       (where applicable) to CSV files for both grid and patch-based region
#'       types.}
#'     \item{\code{save_plots(width = 480, height = 480)}}{Save plots of the
#'       population (staged) totals, the area occupied, total impacts (where
#'       applicable), and total actions as PNG files having specified
#'       \code{width} and \code{height} in pixels.}
#'   }
#' @references
#'   Baker, C. M., Bower, S., Tartaglia, E., Bode, M., Bower, H., & Pressey,
#'   R. L. (2018). Modelling the spread and control of cherry guava on Lord
#'   Howe Island. \emph{Biological Conservation}, 227, 252–258.
#'   \doi{10.1016/j.biocon.2018.09.017}
#'
#'   Bradhurst, R., Spring, D., Stanaway, M., Milner, J., & Kompas, T. (2021).
#'   A generalised and scalable framework for modelling incursions,
#'   surveillance and control of plant and environmental pests.
#'   \emph{Environmental Modelling & Software}, 139, N.PAG.
#'   \doi{10.1016/j.envsoft.2021.105004}
#'
#'   Cacho, O. J., & Hester, S. M. (2022). Modelling biocontrol of invasive
#'   insects: An application to European Wasp (Vespula germanica) in Australia.
#'   \emph{Ecological Modelling}, 467. \doi{10.1016/j.ecolmodel.2022.109939}
#'
#'   García Adeva, J. J., Botha, J. H., & Reynolds, M. (2012). A simulation
#'   modelling approach to forecast establishment and spread of Bactrocera
#'   fruit flies. \emph{Ecological Modelling}, 227, 93–108.
#'   \doi{10.1016/j.ecolmodel.2011.11.026}
#'
#'   Gormley, A. M., Holland, E. P., Barron, M. C., Anderson, D. P., & Nugent,
#'   G. (2016). A modelling framework for predicting the optimal balance
#'   between control and surveillance effort in the local eradication of
#'   tuberculosis in New Zealand wildlife.
#'   \emph{Preventive Veterinary Medicine}, 125, 10–18.
#'   \doi{10.1016/j.prevetmed.2016.01.007}
#'
#'   Krug, R. M., Roura-Pascual, N., & Richardson, D. M. (2010). Clearing of
#'   invasive alien plants under different budget scenarios: using a
#'   simulation model to test efficiency. \emph{Biological Invasions}, 12(12),
#'   4099–4112. \doi{10.1007/s10530-010-9827-3}
#'
#'   Lustig, A., James, A., Anderson, D., & Plank, M. (2019). Pest control at a
#'   regional scale: Identifying key criteria using a spatially explicit,
#'   agent‐based model. \emph{Journal of Applied Ecology}, 56(7 pp.1515–1527),
#'   1527–1515. \doi{10.1111/1365-2664.13387}
#'
#'   Rout, T. M., Moore, J. L., & McCarthy, M. A. (2014). Prevent, search or
#'   destroy? A partially observable model for invasive species management.
#'   \emph{Journal of Applied Ecology}, 51(3), 804–813.
#'   \doi{10.1111/1365-2664.12234}
#'
#'   Spring, D., Croft, L., & Kompas, T. (2017). Look before you treat:
#'   increasing the cost effectiveness of eradication programs with aerial
#'   surveillance. \emph{Biological Invasions}, 19(2), 521.
#'   \doi{10.1007/s10530-016-1292-1}
#'
#'   Wadsworth, R. A., Collingham, Y. C., Willis, S. G., Huntley, B., & Hulme,
#'   P. E. (2000). Simulating the Spread and Management of Alien Riparian
#'   Weeds: Are They Out of Control? \emph{Journal of Applied Ecology}, 37,
#'   28–38. \doi{10.1046/j.1365-2664.2000.00551.x}
#'
#'   Warburton, B., & Gormley, A. M. (2015). Optimising the Application of
#'   Multiple-Capture Traps for Invasive Species Management Using Spatial
#'   Simulation. \emph{PLoS ONE}, 10(3), 1–14.
#'   \doi{10.1371/journal.pone.0120373}
#'
#'   Zub, K., García-Díaz, P., Sankey, S., Eisler, R., & Lambin, X. (2022).
#'   Using a Modeling Approach to Inform Progress Towards Stoat Eradication
#'   From the Orkney Islands. \emph{Frontiers in Conservation Science}, 2.
#'   \doi{10.3389/fcosc.2021.780102}
#' @include Region.R
#' @export
Results <- function(region, population_model,
                    impacts = list(),
                    actions = list(),
                    time_steps = 1,
                    step_duration = 1,
                    step_units = "years",
                    collation_steps = 1,
                    replicates = 1,
                    combine_stages = NULL,
                    class = character(), ...) {
  UseMethod("Results")
}

#' @name Results
#' @export
Results.Region <- function(region, population_model,
                           impacts = list(),
                           actions = list(),
                           time_steps = 1,
                           step_duration = 1,
                           step_units = "years",
                           collation_steps = 1,
                           replicates = 1,
                           combine_stages = NULL,
                           class = character(), ...) {

  # Validate population model
  if (is.null(population_model) || !inherits(population_model, "Population")) {
    stop("Population model must be a 'Population' or inherited class object.",
         call. = FALSE)
  }

  # Validate impact objects
  if (length(impacts) > 0 &&
      (!is.list(impacts) ||
       !all(sapply(impacts, function(i) inherits(i, "Impacts"))))) {
    stop(paste("Impacts must be a list of 'Impacts' or inherited class",
               "objects."), call. = FALSE)
  }

  # Validate action objects
  if (length(actions) > 0 &&
      (!is.list(actions) ||
       !all(sapply(actions, function(i) inherits(i, "Actions"))))) {
    stop(paste("Actions must be a list of 'Actions' or inherited class",
               "objects."), call. = FALSE)
  }

  # Population stages (NULL or number of stages)
  stages <- population_model$get_stages()

  # Stage labels
  if (is.numeric(stages) && is.null(combine_stages)) {
    stage_labels <- attr(population_model$get_growth(), "labels")
  } else if (is.numeric(combine_stages)) {
    stage_labels <- "combined"
  }

  # Include collated and total?
  include_collated <- (region$get_locations() > 1)

  # Include population?
  include_population <- (population_model$get_type() %in%
                           c("unstructured", "stage_structured"))

  # Direct action?
  direct_action <- function(a) {
    return(a$get_label(include_id = FALSE) %in%
             c("detected", "control_search_destroy", "removed"))
  }

  # Individual type action probabilities?
  indiv_type_action <- function(a) {
    action_label <- a$get_label(include_id = FALSE)
    return(population_model$get_type() %in%
             c("unstructured", "stage_structured") &&
             ((action_label == "detected" &&
                 a$get_sensitivity_type() == "individual") ||
                (action_label == "control_search_destroy" &&
                   a$get_manage_pr_type() == "individual") ||
                (action_label == "removed" &&
                   a$get_removal_pr_type() == "individual")))
  }

  # Initialize result lists
  if (include_population) {
    results <- list(population = list())
    if (include_collated) {
      results$total_pop <- list()
    }
    results$occupancy <- list()
  } else {
    results <- list(occupancy = list())
  }
  if (include_collated) {
    results$total_occup <- list()
  }
  results$area <- list()
  if (region$get_type() == "grid") {
    cell_areas <-
      terra::cellSize(region$get_template())[region$get_indices()][,1]
    attr(results$area, "units") <- "square metres"
  } else if (region$spatially_implicit()) {
    attr(results$area, "units") <- "square metres"
  } else {
    attr(results$area, "units") <- "patches"
  }
  if (include_population) {
    if (is.numeric(stages) && is.null(combine_stages)) {
      zeros <- list(population = population_model$make(initial = 0L))
      if (include_collated) {
        zeros$total_pop <- zeros$population[1,, drop = FALSE]
      }
    } else if (is.numeric(stages) && is.numeric(combine_stages)) {
      zeros <- list(population = array(0L, c(region$get_locations(), 1)))
      colnames(zeros$population) <- stage_labels
      if (include_collated) {
        zeros$total_pop <- zeros$population[1,, drop = FALSE]
      }
    } else {
      zeros <- list(population = rep(0L, region$get_locations()))
      if (include_collated) {
        zeros$total_pop <- 0L
      }
    }
    zeros$occupancy <- rep(0L, region$get_locations())
  } else {
    zeros <- list(occupancy = rep(0L, region$get_locations()))
  }
  if (include_collated) {
    zeros$total_occup <- 0L
  }
  zeros$area <- 0L
  if (replicates > 1) { # summaries
    if (include_population) {
      zeros$population <- list(mean = zeros$population, sd = zeros$population)
      if (include_collated) {
        zeros$total_pop <- list(mean = zeros$total_pop, sd = zeros$total_pop)
      }
    }
    zeros$occupancy <- list(mean = zeros$occupancy)
    if (include_collated) {
      zeros$total_occup <- list(mean = zeros$total_occup,
                                sd = zeros$total_occup)
    }
    zeros$area <- list(mean = zeros$area, sd = zeros$area)
  }
  if (include_collated) {
    for (tmc in as.character(c(0, seq(collation_steps, time_steps,
                                      by = collation_steps)))) {
      if (include_population) {
        results$population[[tmc]] <- zeros$population
      }
      results$occupancy[[tmc]] <- zeros$occupancy
    }
    for (tmc in as.character(0:time_steps)) {
      if (include_population) {
        results$total_pop[[tmc]] <- zeros$total_pop
      }
      results$total_occup[[tmc]] <- zeros$total_occup
      results$area[[tmc]] <- zeros$area
    }
  } else {
    for (tmc in as.character(0:time_steps)) {
      if (include_population) {
        results$population[[tmc]] <- zeros$population
      }
      results$occupancy[[tmc]] <- zeros$occupancy
      results$area[[tmc]] <- zeros$area
    }
  }
  rm(zeros)

  # Impact results
  if (length(impacts) > 0) {
    zeros <- list()
    zeros$impact <- rep(0L, region$get_locations())
    if (include_collated) {
      zeros$total_impact <- 0L
    }
    if (replicates > 1) { # summaries
      zeros$impact <- list(mean = zeros$impact, sd = zeros$impact)
      if (include_collated) {
        zeros$total_impact <- list(mean = zeros$total_impact,
                                   sd = zeros$total_impact)
      }
    }
    zeros$impact_steps <- list()
    zeros$total_impact_steps <- list()
    if (include_collated) {
      for (tmc in as.character(c(0, seq(collation_steps, time_steps,
                                        by = collation_steps)))) {
        zeros$impact_steps[[tmc]] <- zeros$impact
      }
      for (tmc in as.character(0:time_steps)) {
        zeros$total_impact_steps[[tmc]] <- zeros$total_impact
      }
    } else {
      for (tmc in as.character(0:time_steps)) {
        zeros$impact_steps[[tmc]] <- zeros$impact
      }
    }
    # Nested impacts result list based on valuation types
    results$impacts <- list(idx = list())
    valuation_types <-
      sapply(impacts, function(impacts_i) impacts_i$get_valuation_type())
    for (vt in unique(valuation_types)) {
      results$impacts[[vt]] <- list()
      results$impacts$idx[[vt]] <- which(vt == valuation_types)
      names(results$impacts$idx[[vt]]) <-
        sapply(results$impacts$idx[[vt]],
               function(i) impacts[[i]]$get_asset_name())
      if (include_collated) {
        results$impacts[[vt]]$total <- list()
      }
    }
    for (i in 1:length(impacts)) {
      vt <- impacts[[i]]$get_valuation_type()
      a <- impacts[[i]]$get_asset_name()
      value_unit <- impacts[[i]]$get_value_unit()
      results$impacts[[vt]][[a]] <- zeros$impact_steps
      attr(results$impacts[[vt]][[a]], "unit") <- value_unit
      if (include_collated) {
        results$impacts[[vt]]$total[[a]] <- zeros$total_impact_steps
        attr(results$impacts[[vt]]$total[[a]], "unit") <- value_unit
      }
    }
    for (vt in names(results$impacts)[-1]) {
      idx <- unname(results$impacts$idx[[vt]])
      value_units <- unlist(
        sapply(impacts[idx], function(impacts_i) impacts_i$get_value_unit()))
      if (length(idx) > 1 && # multiple & same unit
          (all(is.character(value_units)) && !any(value_units == "") &&
           length(unique(value_units)) == 1)) {
        results$impacts[[vt]]$combined <- zeros$impact_steps
        attr(results$impacts[[vt]]$combined, "unit") <- value_units[1]
        if (include_collated) {
          results$impacts[[vt]]$total$combined <- zeros$total_impact_steps
          attr(results$impacts[[vt]]$total$combined, "unit") <- value_units[1]
        }
      }
      if (vt == "monetary") {
        results$impacts[[vt]]$cumulative <- results$impacts[[vt]]
      }
    }
    rm(zeros)
  }

  # Action results
  if (length(actions) > 0) {

    # Individual action results
    results$actions <- lapply(actions, function(actions_i) {

      # Build list of initial zeros
      zeros <- list()

      # Binary indicator of action success
      zeros$action <- rep(FALSE, region$get_locations())
      if (include_collated) {
        zeros$total_action <- 0L
      }

      # Number of individuals when applicable
      include_indiv <- indiv_type_action(actions_i)
      if (include_indiv) {
        if (is.numeric(stages)) {
          if (is.numeric(combine_stages)) {
            zeros$action_num <- array(0L, c(region$get_locations(), 1))
            colnames(zeros$action_num) <- stage_labels
            if (include_collated) {
              zeros$total_action_num <- zeros$action_num[1,, drop = FALSE]
            }
          } else {
            zeros$action_num <- population_model$make(initial = 0L)
            if (include_collated) {
              zeros$total_action_num <- zeros$action_num[1,, drop = FALSE]
            }
          }
        } else {
          zeros$action_num <- rep(0L, region$get_locations())
          if (include_collated) {
            zeros$total_action_num <- 0L
          }
        }
      }

      # Action costs
      if (actions_i$include_cost()) {
        zeros$action_cost <- rep(0L, region$get_locations())
        if (include_collated) {
          zeros$total_action_cost <- 0L
        }
      }

      # Summary statistics when applicable
      if (replicates > 1) {
        zeros$action <- list(mean = +zeros$action)
        if (include_collated) {
          zeros$total_action <- list(mean = zeros$total_action,
                                     sd = zeros$total_action)
        }
        if (include_indiv) {
          zeros$action_num <- list(mean = zeros$action_num,
                                   sd = zeros$action_num)
          if (include_collated) {
            zeros$total_action_num <- list(mean = zeros$total_action_num,
                                           sd = zeros$total_action_num)
          }
        }
        if (actions_i$include_cost()) {
          zeros$action_cost <- list(mean = zeros$action_cost,
                                    sd = zeros$action_cost)
          if (include_collated) {
            zeros$total_action_cost <- list(mean = zeros$total_action_cost,
                                            sd = zeros$total_action_cost)
          }
        }
      }

      # Build initial action list
      actions_list <- list()
      action_label <- actions_i$get_label(include_id = FALSE)
      actions_list[[action_label]] <- list()
      if (include_collated) {
        for (tm in as.character(c(0, seq(collation_steps, time_steps,
                                         by = collation_steps)))) {
          actions_list[[action_label]][[tm]] <- zeros$action
        }
        actions_list$total <- list()
        for (tm in as.character(0:time_steps)) {
          actions_list$total[[tm]] <- zeros$total_action
        }
      } else {
        for (tm in as.character(0:time_steps)) {
          actions_list[[action_label]][[tm]] <- zeros$action
        }
      }
      if (include_indiv) {
        actions_list$number <- list()
        actions_list$number[[action_label]] <- list()
        if (include_collated) {
          for (tm in as.character(c(0, seq(collation_steps, time_steps,
                                           by = collation_steps)))) {
            actions_list$number[[action_label]][[tm]] <- zeros$action_num
          }
          actions_list$number$total <- list()
          for (tm in as.character(0:time_steps)) {
            actions_list$number$total[[tm]] <- zeros$total_action_num
          }
        } else {
          for (tm in as.character(0:time_steps)) {
            actions_list$number[[action_label]][[tm]] <- zeros$action_num
          }
        }
      }
      if (actions_i$include_cost()) {
        actions_list$cost <- list()
        actions_list$cost[[action_label]] <- list()
        if (include_collated) {
          for (tm in as.character(c(0, seq(collation_steps, time_steps,
                                           by = collation_steps)))) {
            actions_list$cost[[action_label]][[tm]] <- zeros$action_cost
          }
          actions_list$cost$total <- list()
          for (tm in as.character(0:time_steps)) {
            actions_list$cost$total[[tm]] <- zeros$total_action_cost
          }
        } else {
          for (tm in as.character(0:time_steps)) {
            actions_list$cost[[action_label]][[tm]] <- zeros$action_cost
          }
        }
        actions_list$cost$cumulative <- actions_list$cost
        attr(actions_list$cost, "unit") <- actions_i$get_cost_unit()
      }
      actions_list
    })

    # Combined cost for multiple actions
    if (length(actions) > 1 &&
        all(sapply(actions, function(a) a$include_cost())) &&
        all(sapply(actions, function(a) a$get_cost_unit()) ==
            actions[[1]]$get_cost_unit())) {
      results$actions$cost <- list(combined = results$actions[[1]]$cost[[1]])
      if (include_collated) {
        results$actions$cost$total <- results$actions[[1]]$cost$total
      }
      results$actions$cost$cumulative <-
        list(combined = results$actions[[1]]$cost$cumulative[[1]])
      attr(results$actions$cost, "unit") <- actions[[1]]$get_cost_unit()
      if (include_collated) {
        results$actions$cost$cumulative$total <-
          results$actions[[1]]$cost$cumulative$total
      }
    }
  }

  # Combined monetary impacts and action costs
  if (length(impacts) > 0 &&
      all(sapply(impacts, function(i) {
        i$get_valuation_type() == "monetary" })) &&
      length(actions) > 0 && is.list(results$actions$cost) &&
      all(sapply(impacts, function(i) i$get_value_unit()) ==
          actions[[1]]$get_cost_unit())) {
    results$cost <- results$actions$cost
  }

  # Create a class structure
  self <- structure(list(), class = c(class, "Results"))

  # Collate results
  self$collate <- function(r, tm, n) {

    # Detach attributes related to impacts and actions
    if (length(impacts)) {
      calc_impacts <- attr(n, "impacts")
      attr(n, "impacts") <- NULL
      # impacts_attr <- attributes(n)[impacts[[1]]$get_attr_names(n)]
      attributes(n)[impacts[[1]]$get_attr_names(n)] <- NULL
      # for (a in names(impacts_attr)) {
      #   attr(n, a) <- impacts_attr[[a]]
      # }
    }

    # Use character list index (allows initial time = 0 and collated times)
    tmc <- as.character(tm)

    # Calculate total area occupied
    occupied_idx <- which(rowSums(as.matrix(n)) > 0)
    if (region$get_type() == "grid") {
      total_area <- sum(cell_areas[occupied_idx])
    } else if (region$spatially_implicit() &&
               is.numeric(attr(n, "diffusion_radius"))) {
      total_area <- pi*(attr(n, "diffusion_radius"))^2
    } else if (region$spatially_implicit() &&
               is.numeric(attr(n, "spread_area"))) {
      total_area <- attr(n, "spread_area")
    } else { # patches
      total_area <- length(occupied_idx)
    }

    # Clean n (remove attributes etc.)
    if (population_model$get_type() == "presence_only") {
      n <- as.integer(n)
    }

    # Combine stages when required
    if (population_model$get_type() == "stage_structured" &&
        is.numeric(combine_stages)) {
      n <- matrix(rowSums(n[,combine_stages, drop = FALSE]), ncol = 1)
      colnames(n) <- stage_labels
    }

    # Calculate occupancy
    occupancy <- +(rowSums(as.matrix(n)) > 0)
    if (include_collated) {
      total_occup <- sum(occupancy)
    }

    # Shape total when population is staged
    if (include_population && include_collated) {
      if (is.numeric(stages)) {
        total_n <- array(colSums(n), c(1, ncol(n)))
        colnames(total_n) <- stage_labels
      } else {
        total_n <- sum(n)
      }
    }

    if (replicates > 1) { # summaries

      # Calculates running mean and standard deviation (note: variance*r is
      # stored as SD and transformed at the final replicate and time step)

      # Population summaries at each location
      if (include_population &&
          (!include_collated || (tm %% collation_steps == 0))) {
        previous_mean <- results$population[[tmc]]$mean
        results$population[[tmc]]$mean <<-
          previous_mean + (n - previous_mean)/r
        previous_sd <- results$population[[tmc]]$sd
        results$population[[tmc]]$sd <<-
          (previous_sd + ((n - previous_mean)*
                            (n - results$population[[tmc]]$mean)))
      }

      # Total population size summaries
      if (include_population && include_collated) {
        previous_mean <- results$total_pop[[tmc]]$mean
        results$total_pop[[tmc]]$mean <<-
          previous_mean + (total_n - previous_mean)/r
        previous_sd <- results$total_pop[[tmc]]$sd
        results$total_pop[[tmc]]$sd <<-
          (previous_sd + ((total_n - previous_mean)*
                            (total_n - results$total_pop[[tmc]]$mean)))
      }

      # Occupancy summaries at each location
      if (!include_collated || (tm %% collation_steps == 0)) {
        previous_mean <- results$occupancy[[tmc]]$mean
        results$occupancy[[tmc]]$mean <<-
          previous_mean + (occupancy - previous_mean)/r
      }

      # Total occupancy summaries
      if (include_collated) {
        previous_mean <- results$total_occup[[tmc]]$mean
        results$total_occup[[tmc]]$mean <<- (previous_mean +
                                               (total_occup - previous_mean)/r)
        previous_sd <- results$total_occup[[tmc]]$sd
        results$total_occup[[tmc]]$sd <<-
          (previous_sd + ((total_occup - previous_mean)*
                            (total_occup - results$total_occup[[tmc]]$mean)))
      }

      # Total area occupied summaries
      previous_mean <- results$area[[tmc]]$mean
      results$area[[tmc]]$mean <<- (previous_mean +
                                      (total_area - previous_mean)/r)
      previous_sd <- results$area[[tmc]]$sd
      results$area[[tmc]]$sd <<-
        (previous_sd + ((total_area - previous_mean)*
                          (total_area - results$area[[tmc]]$mean)))

    } else {

      # Population at each location
      if (include_population &&
          (!include_collated || (tm %% collation_steps == 0))) {
        results$population[[tmc]] <<- n
      }

      # Total population size
      if (include_population && include_collated) {
        results$total_pop[[tmc]] <<- total_n
      }

      # Occupancy at each location
      if (!include_collated || (tm %% collation_steps == 0)) {
        results$occupancy[[tmc]] <<- occupancy
      }

      # Total occupancy
      if (include_collated) {
        results$total_occup[[tmc]] <<- total_occup
      }

      # Total area occupied
      results$area[[tmc]] <<- total_area
    }

    # Collate impacts
    if (length(impacts) > 0) {

      # Place calculated impacts in existing results structure
      for (vt in names(results$impacts)[-1]) { # valuation types

        # Initialise combined impacts
        if ("combined" %in% names(results$impacts[[vt]])) {
          combined_impact <- 0
          if ("total" %in% names(results$impacts[[vt]])) {
            total_combined_impact <- 0
          }
          if (tm == 0 && "cumulative" %in% names(results$impacts[[vt]])) {
            results$impacts[[vt]]$cumulative$combined$current <<- 0
            if ("total" %in% names(results$impacts[[vt]])) {
              results$impacts[[vt]]$cumulative$total$combined$current <<- 0
            }
          }
        }

        # Impact assets
        for (i in unname(results$impacts$idx[[vt]])) {

          # Asset name
          a <- impacts[[i]]$get_asset_name()

          # Current total impact for asset
          if ("total" %in% names(results$impacts[[vt]])) {
            total_impact <- sum(calc_impacts[[i]])
          }

          # Update combined impacts
          if ("combined" %in% names(results$impacts[[vt]])) {
            combined_impact <- combined_impact + calc_impacts[[i]]
            if ("total" %in% names(results$impacts[[vt]])) {
              total_combined_impact <- total_combined_impact + total_impact
            }
          }

          # Update current cumulative impacts
          if ("cumulative" %in% names(results$impacts[[vt]])) {
            if (tm == 0) {
              results$impacts[[vt]]$cumulative[[a]]$current <<- 0
            }
            results$impacts[[vt]]$cumulative[[a]]$current <<-
              (results$impacts[[vt]]$cumulative[[a]]$current +
                 calc_impacts[[i]])
            if ("combined" %in% names(results$impacts[[vt]]$cumulative)) {
              results$impacts[[vt]]$cumulative$combined$current <<-
                (results$impacts[[vt]]$cumulative$combined$current +
                   calc_impacts[[i]])
            }
            if ("total" %in% names(results$impacts[[vt]]$cumulative)) {
              if (tm == 0) {
                results$impacts[[vt]]$cumulative$total[[a]]$current <<- 0
              }
              results$impacts[[vt]]$cumulative$total[[a]]$current <<-
                (results$impacts[[vt]]$cumulative$total[[a]]$current +
                   total_impact)
              if ("combined" %in%
                  names(results$impacts[[vt]]$cumulative$total)) {
                results$impacts[[vt]]$cumulative$total$combined$current <<-
                  (results$impacts[[vt]]$cumulative$total$combined$current +
                     total_impact)
              }
            }
          }

          if (replicates > 1) { # summaries

            # Calculates running mean and standard deviation
            # (note: variance*r is stored as SD and transformed at the final
            #  replicate and time step)

            # Calculated impacts recorded in specified time steps
            if (!include_collated || tm %% collation_steps == 0) {
              previous_mean <- results$impacts[[vt]][[a]][[tmc]]$mean
              results$impacts[[vt]][[a]][[tmc]]$mean <<-
                previous_mean + (calc_impacts[[i]] - previous_mean)/r
              previous_sd <- results$impacts[[vt]][[a]][[tmc]]$sd
              results$impacts[[vt]][[a]][[tmc]]$sd <<-
                (previous_sd + ((calc_impacts[[i]] - previous_mean)*
                                  (calc_impacts[[i]] -
                                     results$impacts[[vt]][[a]][[tmc]]$mean)))
            }

            # Total impacts at every time step
            if ("total" %in% names(results$impacts[[vt]])) {
              previous_mean <- results$impacts[[vt]]$total[[a]][[tmc]]$mean
              results$impacts[[vt]]$total[[a]][[tmc]]$mean <<-
                previous_mean + (total_impact - previous_mean)/r
              previous_sd <- results$impacts[[vt]]$total[[a]][[tmc]]$sd
              results$impacts[[vt]]$total[[a]][[tmc]]$sd <<-
                (previous_sd +
                   ((total_impact - previous_mean)*
                      (total_impact -
                         results$impacts[[vt]]$total[[a]][[tmc]]$mean)))
            }

            # Cumulative impacts
            if ("cumulative" %in% names(results$impacts[[vt]])) {

              # Collated cumulative impacts
              if (!include_collated || tm %% collation_steps == 0) {
                previous_mean <-
                  results$impacts[[vt]]$cumulative[[a]][[tmc]]$mean
                results$impacts[[vt]]$cumulative[[a]][[tmc]]$mean <<-
                  (previous_mean +
                     (results$impacts[[vt]]$cumulative[[a]]$current
                      - previous_mean)/r)
                previous_sd <- results$impacts[[vt]]$cumulative[[a]][[tmc]]$sd
                results$impacts[[vt]]$cumulative[[a]][[tmc]]$sd <<-
                  (previous_sd +
                     ((results$impacts[[vt]]$cumulative[[a]]$current -
                         previous_mean)*
                        (results$impacts[[vt]]$cumulative[[a]]$current -
                           results$impacts[[vt]]$cumulative[[a]][[tmc]]$mean)))
              }

              # Cumulative totals at every time step
              if ("total" %in% names(results$impacts[[vt]]$cumulative)) {
                previous_mean <-
                  results$impacts[[vt]]$cumulative$total[[a]][[tmc]]$mean
                results$impacts[[vt]]$cumulative$total[[a]][[tmc]]$mean <<-
                  (previous_mean +
                     (results$impacts[[vt]]$cumulative$total[[a]]$current -
                        previous_mean)/r)
                previous_sd <-
                  results$impacts[[vt]]$cumulative$total[[a]][[tmc]]$sd
                results$impacts[[vt]]$cumulative$total[[a]][[tmc]]$sd <<-
                  (previous_sd +
                     ((results$impacts[[vt]]$cumulative$total[[a]]$current -
                         previous_mean)*
                        (results$impacts[[vt]]$cumulative$total[[a]]$current -
                           results$impacts[[vt]]$cumulative$total[[a]][[
                             tmc]]$mean)))
              }
            }

          } else {

            # Calculated impacts recorded in specified time steps
            if (!include_collated || tm %% collation_steps == 0) {
              results$impacts[[vt]][[a]][[tmc]] <<- calc_impacts[[i]]
            }

            # Total impacts at every time step
            if ("total" %in% names(results$impacts[[vt]])) {
              results$impacts[[vt]]$total[[a]][[tmc]] <<- total_impact
            }

            # Cumulative impacts
            if ("cumulative" %in% names(results$impacts[[vt]])) {

              # Collated cumulative impacts
              if (!include_collated || tm %% collation_steps == 0) {
                results$impacts[[vt]]$cumulative[[a]][[tmc]] <<-
                  results$impacts[[vt]]$cumulative[[a]]$current
              }

              # Cumulative totals at every time step
              if ("total" %in% names(results$impacts[[vt]]$cumulative)) {
                results$impacts[[vt]]$cumulative$total[[a]][[tmc]] <<-
                  results$impacts[[vt]]$cumulative$total[[a]]$current
              }
            }
          }
        }

        # Collate combined impacts
        if ("combined" %in% names(results$impacts[[vt]])) {
          if (replicates > 1) { # summaries

            # Calculates running mean and standard deviation
            # (note: variance*r is stored as SD and transformed at the final
            #  replicate and time step)

            # Combined impacts recorded in specified time steps
            if (!include_collated || tm %% collation_steps == 0) {
              previous_mean <- results$impacts[[vt]]$combined[[tmc]]$mean
              results$impacts[[vt]]$combined[[tmc]]$mean <<-
                previous_mean + (combined_impact - previous_mean)/r
              previous_sd <- results$impacts[[vt]]$combined[[tmc]]$sd
              results$impacts[[vt]]$combined[[tmc]]$sd <<-
                (previous_sd +
                   ((combined_impact - previous_mean)*
                      (combined_impact -
                         results$impacts[[vt]]$combined[[tmc]]$mean)))
            }

            # Total combined impacts at every time step
            if ("total" %in% names(results$impacts[[vt]])) {
              previous_mean <- results$impacts[[vt]]$total$combined[[tmc]]$mean
              results$impacts[[vt]]$total$combined[[tmc]]$mean <<-
                previous_mean + (total_combined_impact - previous_mean)/r
              previous_sd <- results$impacts[[vt]]$total$combined[[tmc]]$sd
              results$impacts[[vt]]$total$combined[[tmc]]$sd <<-
                (previous_sd +
                   ((total_combined_impact - previous_mean)*
                      (total_combined_impact -
                         results$impacts[[vt]]$total$combined[[tmc]]$mean)))
            }

            # Cumulative combined impacts
            if ("cumulative" %in% names(results$impacts[[vt]])) {

              # Collated cumulative combined impacts
              if (!include_collated || tm %% collation_steps == 0) {
                previous_mean <-
                  results$impacts[[vt]]$cumulative$combined[[tmc]]$mean
                results$impacts[[vt]]$cumulative$combined[[tmc]]$mean <<-
                  (previous_mean +
                     (results$impacts[[vt]]$cumulative$combined$current
                      - previous_mean)/r)
                previous_sd <-
                  results$impacts[[vt]]$cumulative$combined[[tmc]]$sd
                results$impacts[[vt]]$cumulative$combined[[tmc]]$sd <<-
                  (previous_sd +
                     ((results$impacts[[vt]]$cumulative$combined$current -
                         previous_mean)*
                        (results$impacts[[vt]]$cumulative$combined$current -
                           results$impacts[[vt]]$cumulative$combined[[
                             tmc]]$mean)))
              }

              # Cumulative combined totals at every time step
              if ("total" %in% names(results$impacts[[vt]]$cumulative)) {
                previous_mean <-
                  results$impacts[[vt]]$cumulative$total$combined[[tmc]]$mean
                results$impacts[[vt]]$cumulative$total$combined[[tmc]]$mean <<-
                  (previous_mean +
                     (results$impacts[[vt]]$cumulative$total$combined$current -
                        previous_mean)/r)
                previous_sd <-
                  results$impacts[[vt]]$cumulative$total$combined[[tmc]]$sd
                results$impacts[[vt]]$cumulative$total$combined[[tmc]]$sd <<-
                  (previous_sd +
                     ((results$impacts[[
                       vt]]$cumulative$total$combined$current - previous_mean)*
                        (results$impacts[[
                          vt]]$cumulative$total$combined$current -
                           results$impacts[[vt]]$cumulative$total$combined[[
                             tmc]]$mean)))
              }
            }

          } else {

            # Combined impacts recorded in specified time steps
            if (!include_collated || tm %% collation_steps == 0) {
              results$impacts[[vt]]$combined[[tmc]] <<- combined_impact
            }

            # Total combined impacts at every time step
            if ("total" %in% names(results$impacts[[vt]])) {
              results$impacts[[vt]]$total$combined[[tmc]] <<-
                total_combined_impact
            }

            # Cumulative combined impacts
            if ("cumulative" %in% names(results$impacts[[vt]])) {

              # Collated cumulative impacts
              if (!include_collated || tm %% collation_steps == 0) {
                results$impacts[[vt]]$cumulative$combined[[tmc]] <<-
                  results$impacts[[vt]]$cumulative$combined$current
              }

              # Cumulative combined totals at every time step
              if ("total" %in% names(results$impacts[[vt]]$cumulative)) {
                results$impacts[[vt]]$cumulative$total$combined[[tmc]] <<-
                  results$impacts[[vt]]$cumulative$total$combined$current
              }
            }
          }
        }
      }
    }

  }

  # Finalize the results collation
  self$finalize <- function() {

    if (replicates > 1) { # summaries

      # Transform population standard deviations
      if (include_population) {
        for (tmc in names(results$population)) {
          results$population[[tmc]]$sd <<-
            sqrt(results$population[[tmc]]$sd/(replicates - 1))
        }
      }

      # Transform total population standard deviations
      if (include_population && include_collated) {
        for (tmc in names(results$total_pop)) {
          results$total_pop[[tmc]]$sd <<-
            sqrt(results$total_pop[[tmc]]$sd/(replicates - 1))
        }
      }

      # Transform collated total occupancy standard deviations
      if (include_collated) {
        for (tmc in names(results$total_occup)) {
          results$total_occup[[tmc]]$sd <<-
            sqrt(results$total_occup[[tmc]]$sd/(replicates - 1))
        }
      }

      # Transform area occupied standard deviations
      for (tmc in names(results$area)) {
        results$area[[tmc]]$sd <<-
          sqrt(results$area[[tmc]]$sd/(replicates - 1))
      }
    }

    # Add labels to staged populations (again)
    if (population_model$get_type() == "stage_structured") {
      for (tmc in names(results$population)) {
        if (replicates > 1) {
          for (s in names(results$population[[tmc]])) {
            colnames(results$population[[tmc]][[s]]) <<- stage_labels
          }
          if (include_collated) {
            for (s in names(results$total_pop[[tmc]])) {
              colnames(results$total_pop[[tmc]][[s]]) <<- stage_labels
            }
          }
        } else {
          colnames(results$population[[tmc]]) <<- stage_labels
          if (include_collated) {
            colnames(results$total_pop[[tmc]]) <<- stage_labels
          }
        }
      }
    }

    # Finalize impact results
    if (length(impacts) > 0) {

      # Clear working memory for current cumulative impacts
      for (vt in names(results$impacts)[-1]) {
        if ("cumulative" %in% names(results$impacts[[vt]])) {
          for (i in unname(results$impacts$idx[[vt]])) {
            a <- impacts[[i]]$get_asset_name()
            results$impacts[[vt]]$cumulative[[a]]$current <<- NULL
            if ("total" %in% names(results$impacts[[vt]]$cumulative)) {
              results$impacts[[vt]]$cumulative$total[[a]]$current <<- NULL
            }
          }
          if ("combined" %in% names(results$impacts[[vt]]$cumulative)) {
            results$impacts[[vt]]$cumulative$combined$current <<- NULL
          }
          if ("total" %in% names(results$impacts[[vt]]$cumulative) &&
              "combined" %in% names(results$impacts[[vt]]$cumulative$total)) {
            results$impacts[[vt]]$cumulative$total$combined$current <<- NULL
          }
        }
      }

      # Transform impact standard deviations
      if (replicates > 1) { # summaries
        for (vt in names(results$impacts)[-1]) {
          for (i in unname(results$impacts$idx[[vt]])) {
            a <- impacts[[i]]$get_asset_name()
            for (tmc in names(results$impacts[[vt]][[a]])) {
              results$impacts[[vt]][[a]][[tmc]]$sd <<-
                sqrt(results$impacts[[vt]][[a]][[tmc]]$sd/(replicates - 1))
            }
            if ("total" %in% names(results$impacts[[vt]])) {
              for (tmc in names(results$impacts[[vt]]$total[[a]])) {
                results$impacts[[vt]]$total[[a]][[tmc]]$sd <<-
                  sqrt(results$impacts[[vt]]$total[[a]][[tmc]]$sd/
                         (replicates - 1))
              }
            }
          }
          if ("combined" %in% names(results$impacts[[vt]])) {
            for (tmc in names(results$impacts[[vt]]$combined)) {
              results$impacts[[vt]]$combined[[tmc]]$sd <<-
                sqrt(results$impacts[[vt]]$combined[[tmc]]$sd/(replicates - 1))
            }
            if ("total" %in% names(results$impacts[[vt]])) {
              for (tmc in names(results$impacts[[vt]]$total$combined)) {
                results$impacts[[vt]]$total$combined[[tmc]]$sd <<-
                  sqrt(results$impacts[[vt]]$total$combined[[tmc]]$sd/
                         (replicates - 1))
              }
            }
          }
          if ("cumulative" %in% names(results$impacts[[vt]])) {
            for (i in unname(results$impacts$idx[[vt]])) {
              a <- impacts[[i]]$get_asset_name()
              for (tmc in names(results$impacts[[vt]]$cumulative[[a]])) {
                results$impacts[[vt]]$cumulative[[a]][[tmc]]$sd <<-
                  sqrt(results$impacts[[vt]]$cumulative[[a]][[tmc]]$sd/
                         (replicates - 1))
              }
              if ("total" %in% names(results$impacts[[vt]]$cumulative)) {
                for (tmc in
                     names(results$impacts[[vt]]$cumulative$total[[a]])) {
                  results$impacts[[vt]]$cumulative$total[[a]][[tmc]]$sd <<-
                    sqrt(results$impacts[[vt]]$cumulative$total[[a]][[tmc]]$sd/
                           (replicates - 1))
                }
              }
            }
            if ("combined" %in% names(results$impacts[[vt]]$cumulative)) {
              for (tmc in names(results$impacts[[vt]]$cumulative$combined)) {
                results$impacts[[vt]]$cumulative$combined[[tmc]]$sd <<-
                  sqrt(results$impacts[[vt]]$cumulative$combined[[tmc]]$sd/
                         (replicates - 1))
              }
              if ("total" %in% names(results$impacts[[vt]]$cumulative)) {
                for (tmc in
                     names(results$impacts[[vt]]$cumulative$total$combined)) {
                  results$impacts[[vt]]$cumulative$total$combined[[tmc]]$sd <<-
                    sqrt(results$impacts[[vt]]$cumulative$total$combined[[
                      tmc]]$sd/(replicates - 1))
                }
              }
            }
          }
        }
      }
    }

  }

  # Get results list
  self$get_list <- function() {
    return(results)
  }

  # Get simulation parameters
  self$get_params  <- function() {
    return(list(time_steps = time_steps,
                step_duration = step_duration,
                step_units = step_units,
                collation_steps = collation_steps,
                replicates = replicates,
                stages = stages,
                combine_stages = combine_stages))
  }

  # Convert first letter of each word of string to upper-case
  title_case <- function(title_str, all = FALSE) {
    title_str <- as.character(title_str)
    if (all) {
      title_str <- unlist(strsplit(title_str, " ", fixed = TRUE))
    }
    return(paste(paste0(toupper(substr(title_str, 1, 1)),
                        substr(title_str, 2, nchar(title_str))),
                 collapse = " "))
  }

  # Post-fix and metadata for single replicate or summaries
  s_post <- list("", mean = "_mean", sd = "_sd")
  s_data <- list(NULL, mean = "mean", sd = "sd")

  # Save collated results as raster files
  if (include_collated && region$get_type() == "grid") {
    self$save_rasters  <- function(...) {

      # Output and non-zero indicator lists
      output_list <- list()
      nonzero_list <- list()

      # Save population rasters for each summary/time step
      if (include_population) {

        # Replicate summaries or single replicate
        if (replicates > 1) {
          summaries <- c("mean", "sd")
        } else {
          summaries <- ""
        }

        # Create and save a results raster per stage/summary/time step
        if (is.null(stages) || is.numeric(combine_stages)) {
          stages <- 1
        }
        for (i in 1:stages) {

          # Stage post-fix
          if (population_model$get_type() == "stage_structured" &&
              is.null(combine_stages)) {
            ic <- paste0("_stage_", i)
          } else {
            ic <- ""
          }

          for (s in summaries) {

            # Summary post-fix
            if (replicates > 1) {
              sc <- paste0("_", s)
            } else {
              sc <- s
            }

            # Add nested list to output list
            output_key <- paste0("population", ic, sc)
            output_list[[output_key]] <- list()

            # Initialise non-zero indicator
            nonzero_list[[output_key]] <- FALSE

            for (tmc in names(results$population)) {

              # Create result raster and update non-zero indicator
              if (population_model$get_type() == "stage_structured") {
                if (replicates > 1) {
                  output_rast <-
                    region$get_rast(results$population[[tmc]][[s]][,i])
                  nonzero_list[[output_key]] <-
                    (nonzero_list[[output_key]] |
                       sum(results$population[[tmc]][[s]][,i]) > 0)
                } else {
                  output_rast <- region$get_rast(results$population[[tmc]][,i])
                  nonzero_list[[output_key]] <-
                    (nonzero_list[[output_key]] |
                       sum(results$population[[tmc]][,i]) > 0)
                }
                if (is.null(combine_stages)) {
                  names(output_rast) <- stage_labels[i]
                } else if (is.numeric(combine_stages)) {
                  names(output_rast) <- "combined"
                }
              } else {
                if (replicates > 1) {
                  output_rast <-
                    region$get_rast(results$population[[tmc]][[s]])
                  nonzero_list[[output_key]] <-
                    (nonzero_list[[output_key]] |
                       sum(results$population[[tmc]][[s]]) > 0)
                } else {
                  output_rast <- region$get_rast(results$population[[tmc]])
                  nonzero_list[[output_key]] <-
                    (nonzero_list[[output_key]] |
                       sum(results$population[[tmc]]) > 0)
                }
              }

              # Write raster to file and add to output list
              filename <- sprintf(paste0("population%s_t%0",
                                         nchar(as.character(time_steps)),
                                         "d%s.tif"), ic, as.integer(tmc), sc)
              output_list[[output_key]][[tmc]] <-
                terra::writeRaster(output_rast, filename, ...)
            }

            # Add list of population metadata as an attribute
            label <- "population size"
            if (population_model$get_type() == "unstructured") {
              stage <- NULL
            } else if (population_model$get_type() == "stage_structured") {
              if (is.null(combine_stages)) {
                stage <- stage_labels[i]
                label <- paste(stage_labels[i], label)
              } else if (is.numeric(combine_stages)) {
                stage <- "combined"
                label <- paste("combined stage", label)
              }
            }
            if (s == "mean") {
              label <- paste("mean", label)
            } else if (s == "sd") {
              label <- paste(label, "std. dev.")
            }
            attr(output_list[[output_key]], "metadata") <- list(
              category = "population",
              population_type = population_model$get_type(),
              stage = stage,
              summary = s,
              label = title_case(label),
              units = NULL,
              scale_type = "continuous_log10",
              nonzero = nonzero_list[[output_key]]
            )
          }
        }
      }

      # Save occupancy rasters for each summary/time step

      # Postfix
      if (replicates > 1) {
        s <- "mean"
        sc <- "_mean"
      } else {
        s <- sc <- ""
      }

      # Add nested list to output list
      output_key <- paste0("occupancy", sc)
      output_list[[output_key]] <- list()

      # Initialise non-zero indicator
      nonzero_list[[output_key]] <- FALSE

      for (tmc in names(results$occupancy)) {

        # Copy results into a raster & update non-zero indicator
        if (replicates > 1) {
          output_rast <- region$get_rast(results$occupancy[[tmc]][["mean"]])
          nonzero_list[[output_key]] <-
            (nonzero_list[[output_key]] |
               sum(results$occupancy[[tmc]][["mean"]]) > 0)
        } else {
          output_rast <- region$get_rast(results$occupancy[[tmc]])
          nonzero_list[[output_key]] <-
            (nonzero_list[[output_key]] |
               sum(results$occupancy[[tmc]]) > 0)
        }

        # Write raster to file and add to output list
        filename <- sprintf(paste0("occupancy_t%0",
                                   nchar(as.character(time_steps)),
                                   "d%s.tif"), as.integer(tmc), sc)
        output_list[[output_key]][[tmc]] <-
          terra::writeRaster(output_rast, filename, ...)
      }

      # Add list of occupancy metadata as an attribute
      label <- "occupancy"
      if (replicates > 1) {
        label <- paste("mean", label)
        scale_type <- "percent" # 0-1
      } else {
        scale_type <- "discrete" # 0|1
      }
      attr(output_list[[output_key]], "metadata") <- list(
        category = "occupancy",
        population_type = population_model$get_type(),
        stage = NULL,
        summary = s,
        label = title_case(label),
        units = NULL,
        scale_type = scale_type,
        nonzero = nonzero_list[[output_key]]
      )

      # Save impacts
      if (length(impacts) > 0) {

        # Replicate summaries or single replicate
        if (replicates > 1) {
          summaries <- c("mean", "sd")
        } else {
          summaries <- 1
        }

        # Save rasters for each valuation type at each time step
        for (i in names(results$impacts)[-1]) {

          # Impacts post-fix and metadata label
          ic <- paste0("_", i)
          i_label <- paste0(i, " ")
          if (is.character(i)) {
            i_label <- sub("_", " ", i_label, fixed = TRUE)
          }

          # Impact aspects and their cumulative when present
          aspects <- names(results$impacts[[i]])
          aspects <- aspects[which(!aspects %in% c("total", "cumulative"))]
          unit_map <- sapply(aspects, function(a) {
            unit <- attr(results$impacts[[i]][[a]], "unit")
            if (is.character(unit)) {
              unit
            } else {
              ""
            }
          })
          if ("cumulative" %in% names(results$impacts[[i]])) {
            cum_aspects <- names(results$impacts[[i]]$cumulative)
            cum_aspects <- cum_aspects[which(cum_aspects != "total")]
            names(cum_aspects) <- paste0("cum_", cum_aspects)
            aspects <- c(aspects, names(cum_aspects))
            unit_map_cum <- unit_map[cum_aspects]
            names(unit_map_cum) <- names(cum_aspects)
            unit_map <- c(unit_map, unit_map_cum)
          } else {
            cum_aspects <- NULL
          }
          for (a in aspects) {
            for (s in summaries) {

              # Add nested list to output list
              if (a %in% names(cum_aspects)) {
                output_key <- paste0("cumulative_impacts", ic, "_",
                                     cum_aspects[a], s_post[[s]])
              } else {
                output_key <- paste0("impacts", ic, "_", a, s_post[[s]])
              }
              output_list[[output_key]] <- list()

              # Initialise non-zero indicator
              nonzero_list[[output_key]] <- FALSE

              # Collated time steps
              if (a %in% names(cum_aspects)) {
                collated_tmc <-
                  names(results$impacts[[i]]$cumulative[[cum_aspects[a]]])
              } else {
                collated_tmc <- names(results$impacts[[i]][[a]])
              }
              for (tmc in collated_tmc) {

                # Copy impacts into a raster & update non-zero indicator
                if (replicates > 1) {
                  if (a %in% names(cum_aspects)) {
                    output_rast <-
                      region$get_rast(results$impacts[[i]]$cumulative[[
                        cum_aspects[a]]][[tmc]][[s]])
                    nonzero_list[[output_key]] <-
                      (nonzero_list[[output_key]] |
                         sum(results$impacts[[i]]$cumulative[[
                           cum_aspects[a]]][[tmc]][[s]]) > 0)
                  } else {
                    output_rast <-
                      region$get_rast(results$impacts[[i]][[a]][[tmc]][[s]])
                    nonzero_list[[output_key]] <-
                      (nonzero_list[[output_key]] |
                         sum(results$impacts[[i]][[a]][[tmc]][[s]]) > 0)
                  }
                } else {
                  if (a %in% names(cum_aspects)) {
                    output_rast <-
                      region$get_rast(results$impacts[[i]]$cumulative[[
                        cum_aspects[a]]][[tmc]])
                    nonzero_list[[output_key]] <-
                      (nonzero_list[[output_key]] |
                         sum(results$impacts[[i]]$cumulative[[
                           cum_aspects[a]]][[tmc]]) > 0)
                  } else {
                    output_rast <-
                      region$get_rast(results$impacts[[i]][[a]][[tmc]])
                    nonzero_list[[output_key]] <-
                      (nonzero_list[[output_key]] |
                         sum(results$impacts[[i]][[a]][[tmc]]) > 0)
                  }
                }

                # Write raster to file
                if (a %in% names(cum_aspects)) {
                  filename <- sprintf(
                    paste0("cumulative_impacts%s_%s_t%0",
                           nchar(as.character(time_steps)), "d%s.tif"),
                    ic, cum_aspects[a], as.integer(tmc), s_post[[s]])
                } else {
                  filename <- sprintf(
                    paste0("impacts%s_%s_t%0", nchar(as.character(time_steps)),
                           "d%s.tif"),
                    ic, a, as.integer(tmc), s_post[[s]])
                }
                output_list[[output_key]][[tmc]] <-
                  terra::writeRaster(output_rast, filename, ...)
              }

              # Add list of impacts metadata as an attribute
              if (a %in% names(cum_aspects)) {
                aspect <- unname(cum_aspects[a])
                cumulative <- TRUE
                cumulative_label <- "cumulative "
              } else {
                aspect <- a
                cumulative <- FALSE
                cumulative_label <- ""
              }
              label <- paste0(cumulative_label, i_label, aspect, " impacts")
              if (s == "mean") {
                label <- paste("mean", label)
              } else if (s == "sd") {
                label <- paste(label, "std. dev.")
              }
              impact_type <- sapply(results$impacts$idx[[i]],
                                    function(j) impacts[[j]]$get_impact_type())
              if (length(impact_type) > 1 && aspect == "combined") {
                impact_type <- paste(impact_type, collapse = " & ")
              } else {
                impact_type <- unname(impact_type[aspect])
              }
              impact_type <- paste0(impact_type, "-based")
              attr(output_list[[output_key]], "metadata") <- list(
                category = "impact",
                type = i,
                name = aspect,
                impact = impact_type,
                cost = FALSE,
                cumulative = cumulative,
                summary = s_data[[s]],
                label = title_case(label),
                units = unname(unit_map[a]),
                scale_type = "continuous",
                nonzero = nonzero_list[[output_key]]
              )
            }
          }
        }
      }

      # Return output list as multi-layer rasters
      return(lapply(output_list, function(rast_list) {
        raster_layers <- terra::rast(rast_list)
        attr(raster_layers, "metadata") <- attr(rast_list, "metadata")
        raster_layers
      }))
    }
  }

  # Save collated (patch only) and summary (both) results as CSV files
  self$save_csv  <- function() {

    # Filename parts summaries or single replicate
    s_fname <- list("", mean = "_mean", sd = "_sd")

    # Time step labels
    time_steps_labels <- sprintf(
      paste0("t%0", nchar(as.character(time_steps)), "d"),
      as.integer(0:time_steps))
    collated_labels <- sprintf(
      paste0("t%0", nchar(as.character(time_steps)), "d"),
      c(0, seq(collation_steps, time_steps, by = collation_steps)))

    # Resolve result stages
    result_stages <- stages
    if (is.null(stages) || is.numeric(combine_stages)) {
      result_stages <- 1
      i_fname <- ""
    } else {
      i_fname <- paste0("_stage_", 1:result_stages)
    }

    # Results for single or multi-patch only
    if (region$get_type() == "patch") {

      # Location coordinates and labels
      if (include_collated) {
        coords <- region$get_coords(extra_cols = TRUE)
        coords <- coords[, c("lon", "lat",
                             names(which(sapply(coords, is.character))))]
      }

      # Save population CSV file(s)
      if (include_population) {
        if (include_collated) {
          if (replicates > 1) {
            summaries <- c("mean", "sd")
          } else {
            summaries <- 1
          }
          for (i in 1:result_stages) {
            output_df <- list()
            for (s in summaries) {
              if (population_model$get_type() == "stage_structured") {
                if (replicates > 1) {
                  output_df[[s]] <- lapply(results$population,
                                           function(c_tm) c_tm[[s]][,i])
                } else {
                  output_df[[s]] <- lapply(results$population,
                                           function(c_tm) c_tm[,i])
                }
              } else {
                if (replicates > 1) {
                  output_df[[s]] <- lapply(results$population,
                                           function(c_tm) c_tm[[s]])
                } else {
                  output_df[[s]] <- results$population
                }
              }
              names(output_df[[s]]) <- collated_labels
              output_df[[s]] <- cbind(coords, as.data.frame(output_df[[s]]))
              filename <- sprintf(paste0("population%s%s.csv"), i_fname[i],
                                  s_fname[[s]])
              utils::write.csv(output_df[[s]], filename, row.names = FALSE)
            }
          }
        } else {
          if (population_model$get_type() == "stage_structured") {
            if (replicates > 1) {
              for (i in 1:result_stages) {
                output_df <- sapply(results$population,
                                    function(pop) as.data.frame(
                                      lapply(pop, function(m) m[,i])))
                colnames(output_df) <- time_steps_labels
                filename <- sprintf(paste0("population%s.csv"), i_fname[i])
                utils::write.csv(output_df, filename, row.names = TRUE)
              }
            } else {
              if (is.numeric(combine_stages)) {
                output_df <- as.data.frame(results$population)
                if (length(combine_stages) == 1) {
                  rownames(output_df) <- attr(population_model$get_growth(),
                                              "labels")[combine_stages]
                } else {
                  rownames(output_df) <- sprintf(
                    "stages %s-%s", min(combine_stages), max(combine_stages))
                }
              } else {
                output_df <- sapply(results$population,
                                    function(pop) as.data.frame(pop))
              }
              colnames(output_df) <- time_steps_labels
              utils::write.csv(output_df, "population.csv", row.names = TRUE)
            }
          } else {
            if (replicates > 1) {
              output_df <- sapply(results$population,
                                  function(pop) as.data.frame(pop))
            } else {
              output_df <- as.data.frame(results$population)
              rownames(output_df) <- "population"
            }
            colnames(output_df) <- time_steps_labels
            utils::write.csv(output_df, "population.csv", row.names = TRUE)
          }
        }
      }

      # Save occupancy CSV file(s)
      if (include_collated) {
        if (replicates > 1) {
          output_df <- lapply(results$occupancy,
                              function(o_tm) o_tm[["mean"]])
          filename <- "occupancy_mean.csv"
        } else {
          output_df <- results$occupancy
          filename <- "occupancy.csv"
        }
        names(output_df) <- collated_labels
        output_df <- cbind(coords, as.data.frame(output_df))
        utils::write.csv(output_df, filename, row.names = FALSE)
      } else {
        if (replicates > 1) {
          output_df <- as.data.frame(lapply(results$occupancy,
                                            function(occ) occ[["mean"]]))
          rownames(output_df) <- "mean"
        } else {
          output_df <- as.data.frame(results$occupancy)
          rownames(output_df) <- "occupancy"
        }
        colnames(output_df) <- time_steps_labels
        utils::write.csv(output_df, "occupancy.csv", row.names = TRUE)
      }
    }

    # Save total population CSV file(s)
    if (include_population && include_collated) {
      total_pop <- list()
      if (replicates > 1 && is.numeric(stages)) {
        for (i in 1:result_stages) {
          total_pop[[i]] <-
            sapply(results$total_pop,
                   function(tot) as.data.frame(lapply(tot, function(m) m[,i])))
          colnames(total_pop[[i]]) <- time_steps_labels
        }
      } else {
        if (replicates > 1 || result_stages > 1) {
          total_pop[[1]] <- sapply(results$total_pop,
                                   function(tot) as.data.frame(tot))
        } else {
          total_pop[[1]] <- matrix(results$total_pop, ncol = time_steps + 1)
          if (population_model$get_type() == "stage_structured" &&
              is.numeric(combine_stages)) {
            if (length(combine_stages) == 1) {
              rownames(total_pop[[1]]) <-
                attr(population_model$get_growth(), "labels")[combine_stages]
            } else {
              rownames(total_pop[[1]]) <-
                sprintf("stages %s-%s", min(combine_stages),
                        max(combine_stages))
            }
          } else {
            rownames(total_pop[[1]]) <- "population"
          }
        }
        colnames(total_pop[[1]]) <- time_steps_labels
      }
      if (length(total_pop) > 1) {
        for (i in 1:length(total_pop)) {
          utils::write.csv(total_pop[[i]],
                           sprintf("total_population%s.csv", i_fname[i]),
                           row.names = TRUE)
        }
      } else {
        utils::write.csv(total_pop[[1]], "total_population.csv",
                         row.names = TRUE)
      }
    }

    # Save total occupancy CSV file
    if (include_collated) {
      if (replicates > 1) {
        total_occup <- sapply(results$total_occup,
                              function(occup) as.data.frame(occup))
      } else {
        total_occup <- matrix(results$total_occup, ncol = time_steps + 1)
        rownames(total_occup) <- "occupancy"
      }
      colnames(total_occup) <- time_steps_labels
      utils::write.csv(total_occup, "total_occupancy.csv", row.names = TRUE)
    }

    # Save (total) area occupied CSV file
    if (replicates > 1) {
      area <- sapply(results$area,
                     function(area) as.data.frame(area))
    } else {
      area <- matrix(results$area, ncol = time_steps + 1)
      rownames(area) <- "area"
    }
    colnames(area) <- time_steps_labels
    if (include_collated) {
      utils::write.csv(area, "total_area_occupied.csv", row.names = TRUE)
    } else {
      utils::write.csv(area, "area_occupied.csv", row.names = TRUE)
    }

    # Save impacts
    if (length(impacts) > 0) {

      # Results for each impact valuation type
      for (i in names(results$impacts)[-1]) {

        # Include impact name
        ic <- paste0("_", i)

        # Results for single or multi-patch only
        if (region$get_type() == "patch") {

          # Impact aspects and their cumulative when present
          aspects <- names(results$impacts[[i]])
          aspects <- aspects[which(!aspects %in% c("total", "cumulative"))]
          if ("cumulative" %in% names(results$impacts[[i]])) {
            cum_aspects <- names(results$impacts[[i]]$cumulative)
            cum_aspects <- cum_aspects[which(cum_aspects != "total")]
            names(cum_aspects) <- paste0("cum_", cum_aspects)
            aspects <- c(aspects, names(cum_aspects))
          } else {
            cum_aspects <- NULL
          }

          # Save impact CSV file(s)
          if (include_collated) { # spatial with coordinates
            for (a in aspects) {
              output_df <- list()
              if (replicates > 1) {
                summaries <- c("mean", "sd")
              } else {
                summaries <- 1
              }
              for (s in summaries) {
                if (replicates > 1) {
                  if (a %in% names(cum_aspects)) {
                    output_df[[s]] <-
                      lapply(results$impacts[[i]]$cumulative[[cum_aspects[a]]],
                             function(a_tm) a_tm[[s]])
                  } else {
                    output_df[[s]] <- lapply(results$impacts[[i]][[a]],
                                             function(a_tm) a_tm[[s]])
                  }
                } else {
                  if (a %in% names(cum_aspects)) {
                    output_df[[s]] <-
                      results$impacts[[i]]$cumulative[[cum_aspects[a]]]
                  } else {
                    output_df[[s]] <- results$impacts[[i]][[a]]
                  }
                }
                names(output_df[[s]]) <- collated_labels
                output_df[[s]] <- cbind(coords, as.data.frame(output_df[[s]]))
                if (a %in% names(cum_aspects)) {
                  filename <- sprintf("cumulative_impacts%s_%s%s.csv", ic,
                                      cum_aspects[a], s_post[[s]])
                } else {
                  filename <- sprintf("impacts%s_%s%s.csv", ic, a,
                                      s_post[[s]])
                }
                utils::write.csv(output_df[[s]], filename, row.names = FALSE)
              }
            }

          } else { # spatially implicit

            # Collect and save to CSV
            if (replicates > 1) {

              # Save mean/SD (rows) for each impact
              for (a in aspects) {
                if (a %in% names(cum_aspects)) {
                  output_df <- sapply(
                    results$impacts[[i]]$cumulative[[cum_aspects[a]]],
                    function(imp) as.data.frame(imp))
                  filename <- sprintf("cumulative_impacts%s_%s.csv", ic,
                                      cum_aspects[a])
                } else {
                  output_df <- sapply(results$impacts[[i]][[a]],
                                      function(imp) as.data.frame(imp))
                  filename <- sprintf("impacts%s_%s.csv", ic, a)
                }
                colnames(output_df) <- time_steps_labels
                utils::write.csv(output_df, filename, row.names = TRUE)
              }

            } else {

              # Saves impacts in rows
              a_i <- which(names(results$impacts[[i]]) != "cumulative")
              output_df <- t(sapply(results$impacts[[i]][a_i], function(a) a))
              colnames(output_df) <- time_steps_labels
              filename <- sprintf("impacts%s.csv", ic)
              utils::write.csv(output_df, filename, row.names = TRUE)
              if ("cumulative" %in% names(results$impacts[[i]])) {
                output_df <- t(sapply(results$impacts[[i]]$cumulative,
                                      function(a) unlist(a)))
                colnames(output_df) <- time_steps_labels
                filename <- sprintf("cumulative_impacts%s.csv", ic)
                utils::write.csv(output_df, filename, row.names = TRUE)
              }
            }
          }
        }

        # Impact totals when present
        if (include_collated) {
          if (replicates > 1) {

            # Save mean/SD (rows) impact totals for each impact
            if ("total" %in% names(results$impacts[[i]])) {
              for (a in names(results$impacts[[i]]$total)) {
                output_df <- as.data.frame(
                  sapply(results$impacts[[i]]$total[[a]],
                         function(tot) unlist(tot)))
                colnames(output_df) <- time_steps_labels
                filename <- sprintf("total_impacts%s_%s.csv", ic, a)
                utils::write.csv(output_df, filename, row.names = TRUE)
              }
            }

            # Save mean/SD (rows) cumulative impact totals for each impact
            if ("cumulative" %in% names(results$impacts[[i]]) &&
                "total" %in% names(results$impacts[[i]]$cumulative)) {
              for (a in names(results$impacts[[i]]$cumulative$total)) {
                output_df <- as.data.frame(
                  sapply(results$impacts[[i]]$cumulative$total[[a]],
                         function(tot) unlist(tot)))
                colnames(output_df) <- time_steps_labels
                filename <- sprintf("total_cumulative_impacts%s_%s.csv", ic, a)
                utils::write.csv(output_df, filename, row.names = TRUE)
              }
            }

          } else {

            # Save impact totals with impacts in rows
            if ("total" %in% names(results$impacts[[i]])) {
              output_df <- as.data.frame(t(
                sapply(results$impacts[[i]]$total,
                       function(a) unlist(a))))
              colnames(output_df) <- time_steps_labels
              filename <- sprintf("total_impacts%s.csv", ic)
              utils::write.csv(output_df, filename, row.names = TRUE)
            }

            # Save cumulative impact totals with impacts in rows
            if ("cumulative" %in% names(results$impacts[[i]]) &&
                "total" %in% names(results$impacts[[i]]$cumulative)) {
              output_df <- as.data.frame(t(
                sapply(results$impacts[[i]]$cumulative$total,
                       function(a) unlist(a))))
              colnames(output_df) <- time_steps_labels
              filename <- sprintf("total_cumulative_impacts%s.csv", ic)
              utils::write.csv(output_df, filename, row.names = TRUE)
            }
          }
        }
      }
    } # impacts

  }

  # Plot total population (per stage) and area occupied as PNG files
  self$save_plots  <- function(width = 480, height = 480) {

    # Resolve the number of (combined) stages used in the results
    result_stages <- stages
    if (is.null(stages) || is.numeric(combine_stages)) {
      result_stages <- 1
    }

    # Stage label for plot headings and files
    stage_label <- ""
    stage_file <- ""
    if (population_model$get_type() == "stage_structured") {
      if (is.numeric(stages) && is.null(combine_stages)) {
        stage_label <- paste0(stage_labels, " ")
        stage_file <- paste0("_stage_", 1:result_stages)
      } else if (is.numeric(stages) && is.numeric(combine_stages)) {
        if (length(combine_stages) == 1) {
          stage_label <- paste0(
            attr(population_model$get_growth(), "labels")[combine_stages], " ")
        } else {
          stage_label <- paste0(sprintf("stages %s-%s", min(combine_stages),
                                        max(combine_stages)), " ")
        }
      }
    }

    # Label population and totals appropriately
    pop_label <- "population size"
    if (region$spatially_implicit()) {
      tot_label <- list(text = "", file = "")
    } else {
      tot_label <- list(text = "total ", file = "total_")
    }

    # All plots have time steps on x-axis
    plot_x_label <- paste0("Time steps (", step_units, ")")

    if (replicates > 1) { # plot summary mean +/- 2 SD

      # Collect total population
      if (include_population) {
        population <- list()
        for (s in c("mean", "sd")) {
          if (include_collated) {
            population[[s]] <- array(sapply(results$total_pop,
                                            function(pop) pop[[s]]),
                                     c(result_stages, time_steps + 1))
          } else {
            population[[s]] <- array(sapply(results$population,
                                            function(pop) pop[[s]]),
                                     c(result_stages, time_steps + 1))
          }
        }
      }

      # Collect total occupancy
      occup <- list()
      if (include_collated) {
        for (s in c("mean", "sd")) {
          occup[[s]] <- sapply(results$total_occup, function(occup) occup[[s]])
        }
      } else {
        occup[["mean"]] <- sapply(results$occupancy,
                                  function(occup) occup[["mean"]])
      }

      # Collect total area
      area <- list()
      for (s in c("mean", "sd")) {
        area[[s]] <- sapply(results$area, function(area) area[[s]])
      }

      # Plot population (per result stage)
      if (include_population) {
        for (s in 1:result_stages) {
          filename <- paste0(tot_label$file, "population", stage_file[s], ".png")
          main_title <- paste0(tot_label$text, stage_label[s], pop_label,
                               " (mean +/- 2 SD)")
          y_label <- title_case(pop_label)
          grDevices::png(filename = filename, width = width, height = height)
          graphics::plot(0:time_steps, population$mean[s,], type = "l",
                         main = title_case(main_title),
                         xlab = plot_x_label, ylab = y_label,
                         ylim = c(0, 1.1*max(population$mean[s,] +
                                               2*population$sd[s,])))
          graphics::lines(0:time_steps,
                          population$mean[s,] + 2*population$sd[s,],
                          lty = "dashed")
          graphics::lines(0:time_steps,
                          pmax(0, population$mean[s,] - 2*population$sd[s,]),
                          lty = "dashed")
          invisible(grDevices::dev.off())
        }
      }

      # Plot occupancy
      if (region$spatially_implicit()) {
        filename <- "occupancy.png"
        main_title <- "occupancy (mean)"
      } else {
        filename <- "total_occupancy.png"
        main_title <- "total occupancy (mean +/- 2 SD)"
      }
      grDevices::png(filename = filename, width = width, height = height)
      if (include_collated) {
        ylim <- c(0, 1.1*max(occup$mean + 2*occup$sd))
      } else {
        ylim <- c(0, 1.1*max(occup$mean))
      }
      graphics::plot(0:time_steps, occup$mean, type = "l",
                     main = title_case(main_title),
                     xlab = plot_x_label,
                     ylab = "Occupancy",
                     ylim = ylim)
      if (include_collated) {
        graphics::lines(0:time_steps, occup$mean + 2*occup$sd, lty = "dashed")
        graphics::lines(0:time_steps, pmax(0, occup$mean - 2*occup$sd),
                        lty = "dashed")
      }
      invisible(grDevices::dev.off())

      # Plot area
      if (region$spatially_implicit()) {
        filename <- "area_occupied.png"
        main_title <- "area occupied (mean +/- 2 SD)"
      } else {
        filename <- "total_area_occupied.png"
        main_title <- "total area occupied (mean +/- 2 SD)"
      }
      grDevices::png(filename = filename, width = width, height = height)
      graphics::plot(0:time_steps, area$mean, type = "l",
                     main = title_case(main_title),
                     xlab = plot_x_label,
                     ylab = paste0("Area (", attr(results$area, "units"),
                                   ")"),
                     ylim = c(0, 1.1*max(area$mean + 2*area$sd)))
      graphics::lines(0:time_steps, area$mean + 2*area$sd, lty = "dashed")
      graphics::lines(0:time_steps, pmax(0, area$mean - 2*area$sd),
                      lty = "dashed")
      invisible(grDevices::dev.off())

    } else { # plot values

      # Collect total population
      if (include_population) {
        if (include_collated) {
          population <- array(sapply(results$total_pop, function(pop) pop),
                              c(result_stages, time_steps + 1))
        } else {
          population <- array(sapply(results$population, function(pop) pop),
                              c(result_stages, time_steps + 1))
        }
      }

      # Collect total occupancy
      if (include_collated) {
        occup <- sapply(results$total_occup, function(occup) occup)
      } else {
        occup <- sapply(results$occupancy, function(occup) occup)
      }

      # Collect total area
      area <- sapply(results$area, function(area) area)

      # Plot population (per result stage)
      if (include_population) {
        for (s in 1:result_stages) {
          filename <- paste0(tot_label$file, "population", stage_file[s], ".png")
          main_title <- paste0(tot_label$text, stage_label[s], pop_label)
          y_label <- title_case(pop_label)
          grDevices::png(filename = filename, width = width, height = height)
          graphics::plot(0:time_steps, population[s,], type = "l",
                         main = title_case(main_title),
                         xlab = plot_x_label, ylab = y_label,
                         ylim = c(0, 1.1*max(population[s,])))
          invisible(grDevices::dev.off())
        }
      }

      # Plot total occupancy
      if (region$spatially_implicit()) {
        filename <- "occupancy.png"
        main_title <- "occupancy"
      } else {
        filename <- "total_occupancy.png"
        main_title <- "total occupancy"
      }
      grDevices::png(filename = filename, width = width, height = height)
      graphics::plot(0:time_steps, occup, type = "l",
                     main = title_case(main_title),
                     xlab = plot_x_label,
                     ylab = "Occupancy",
                     ylim = c(0, 1.1*max(occup)))
      invisible(grDevices::dev.off())

      # Plot area
      if (region$spatially_implicit()) {
        filename <- "area_occupied.png"
        main_title <- "area occupied"
      } else {
        filename <- "total_area_occupied.png"
        main_title <- "total area occupied"
      }
      grDevices::png(filename = filename, width = width, height = height)
      graphics::plot(0:time_steps, area, type = "l",
                     main = title_case(main_title),
                     xlab = plot_x_label,
                     ylab = paste0("Area (", attr(results$area, "units"),
                                   ")"),
                     ylim = c(0, 1.1*max(area)))
      invisible(grDevices::dev.off())
    }

    # Plot impacts when present
    if (length(impacts) > 0) {

      # Plots for each impact valuation type
      for (i in names(results$impacts)[-1]) {

        # Include impact name
        ic <- c(paste0("_", i), paste0(i, " "))
        ic[2] <- sub("_", " ", ic[2], fixed = TRUE)

        # Impact aspects and their cumulative when present
        aspects <- names(results$impacts[[i]])
        aspects <- aspects[which(!aspects %in% c("total", "cumulative"))]
        if ("cumulative" %in% names(results$impacts[[i]])) {
          cum_aspects <- names(results$impacts[[i]]$cumulative)
          cum_aspects <- cum_aspects[which(cum_aspects != "total")]
          names(cum_aspects) <- paste0("cum_", cum_aspects)
          aspects <- c(aspects, names(cum_aspects))
        } else {
          cum_aspects <- NULL
        }

        # Plot for each impact or cumulative impact
        for (a in aspects) {

          # Get the values to be plotted and their unit
          if (include_collated) { # use totals
            if (a %in% names(cum_aspects)) {
              a_cum <- cum_aspects[a]
              values <- sapply(results$impacts[[i]]$cumulative$total[[a_cum]],
                               function(tot) unlist(tot))
              unit <- attr(results$impacts[[i]]$cumulative$total[[a_cum]],
                            "unit")
            } else {
              values <- sapply(results$impacts[[i]]$total[[a]],
                               function(tot) unlist(tot))
              unit <- attr(results$impacts[[i]]$total[[a]], "unit")
            }
          } else { # spatially-implicit
            if (a %in% names(cum_aspects)) {
              a_cum <- cum_aspects[a]
              values <- sapply(results$impacts[[i]]$cumulative[[a_cum]],
                               function(val) unlist(val))
              unit <- attr(results$impacts[[i]]$cumulative[[a_cum]], "unit")
            } else {
              values <- sapply(results$impacts[[i]][[a]],
                               function(val) unlist(val))
              unit <- attr(results$impacts[[i]][[a]], "unit")
            }
          }
          if (is.null(unit) || unit == "") {
            unit <- ""
          } else {
            unit <- paste0(" (", unit, ")")
          }

          # Impact labels
          if (a %in% names(cum_aspects)) {
            ac <- c(paste0("_", cum_aspects[a]), paste0(cum_aspects[a], " "))
            cc <- c("cumulative_", "cumulative ")
          } else {
            ac <- c(paste0("_", a), paste0(a, " "))
            cc <- c("", "")
          }

          # Plot file names and main titles
          if (include_collated) {
            filename <- sprintf("total_%simpacts%s%s.png", cc[1], ic[1], ac[1])
            main_title <- sprintf("total %s%s%simpacts", cc[2], ic[2], ac[2])
          } else {
            filename <- sprintf("%simpacts%s%s.png", cc[1], ic[1], ac[1])
            main_title <- sprintf("%s%s%simpacts", cc[2], ic[2], ac[2])
          }

          # Create and save plot
          if (replicates > 1) { # plot summary mean +/- 2 SD
            main_title <- paste(main_title, "(mean +/- 2 SD)")
            values <- as.data.frame(t(values))
            grDevices::png(filename = filename, width = width,
                           height = height)
            graphics::plot(0:time_steps, values$mean, type = "l",
                           main = title_case(main_title),
                           xlab = plot_x_label,
                           ylab = sprintf("Impact%s", unit),
                           ylim = c(0, 1.1*max(values$mean + 2*values$sd)))
            graphics::lines(0:time_steps, values$mean + 2*values$sd,
                            lty = "dashed")
            graphics::lines(0:time_steps,
                            pmax(0, values$mean - 2*values$sd),
                            lty = "dashed")
            invisible(grDevices::dev.off())
          } else {
            grDevices::png(filename = filename, width = width,
                           height = height)
            graphics::plot(0:time_steps, values, type = "l",
                           main = title_case(main_title),
                           xlab = plot_x_label,
                           ylab = sprintf("Impact%s", unit),
                           ylim = c(0, 1.1*max(values)))
            invisible(grDevices::dev.off())
          }
        }
      }
    } # impacts

  }

  return(self)
}
