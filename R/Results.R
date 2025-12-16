#' Results class builder
#'
#' Builds a class for encapsulating, calculating, and collating spread
#' simulation results, including the population at each location at each
#' collation time step, the total population size, and the area occupied at
#' each time step. When simulations are replicated, summary results (means and
#' standard deviations) are produced.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the spread simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation and growth functionality for the
#'   spread simulations.
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
#'       replicate \code{r} and time step \code{tm} using the current vector or
#'       array \code{n} representing the population at each location.}
#'     \item{\code{finalize()}}{Finalize the results collation (summary
#'       calculations).}
#'     \item{\code{as_list()}}{Return the results as a list (collated, total,
#'       area).}
#'     \item{\code{get_params()}}{Get the simulation parameters used.}
#'     \item{\code{save_rasters(...)}}{Save the collated results as raster TIF
#'       files when the region is grid-based. \code{Terra} raster write options
#'       may be passed to the function. Returns a list of saved \code{Terra}
#'       raster layers with attached \code{"metadata"} attribute list, which
#'       includes: plot category; population type & stage (when applicable);
#'       multi-replicate summary (when applicable); a descriptive label; units;
#'       scale type; and an indication of whether or not each layer contains
#'       non-zero values.}
#'     \item{\code{save_csv()}}{Save the collated results as comma-separated
#'       values (CSV) files when the region is patch-based. Also saves the
#'       population totals and area occupied to CSV files for both grid and
#'       patch-based region types.}
#'     \item{\code{save_plots(width = 480, height = 480)}}{Save plots of the
#'       population (staged) totals and the area occupied as PNG files having
#'       specified \code{width} and \code{height} in pixels.}
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
#' @include Region.R
#' @export
Results <- function(region, population_model,
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
    for (tm in as.character(c(0, seq(collation_steps, time_steps,
                                     by = collation_steps)))) {
      if (include_population) {
        results$population[[tm]] <- zeros$population
      }
      results$occupancy[[tm]] <- zeros$occupancy
    }
    for (tm in as.character(0:time_steps)) {
      if (include_population) {
        results$total_pop[[tm]] <- zeros$total_pop
      }
      results$total_occup[[tm]] <- zeros$total_occup
      results$area[[tm]] <- zeros$area
    }
  } else {
    for (tm in as.character(0:time_steps)) {
      if (include_population) {
        results$population[[tm]] <- zeros$population
      }
      results$occupancy[[tm]] <- zeros$occupancy
      results$area[[tm]] <- zeros$area
    }
  }
  rm(zeros)

  # Create a class structure
  self <- structure(list(), class = c(class, "Results"))

  # Collate results
  self$collate <- function(r, tm, n) {

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
                label <- paste("stage", as.character(stage), label)
              } else if (is.numeric(combine_stages)) {
                stage <- "combined"
                label <- paste("combined stage", label)
              }
            }
            if (s == "mean") {
              label <- paste("mean", label)
            } else if (s == "sd") {
              label <- paste(label, "standard deviation")
            }
            attr(output_list[[output_key]], "metadata") <- list(
              category = "population",
              population_type = population_model$get_type(),
              stage = stage,
              summary = s,
              label = title_case(label),
              units = NULL,
              scale_type = "continuous",
              nonzero = nonzero_list[[output_key]]
            )
          }
        }
      }

      # Save occupancy rasters for each summary/time step

      # Postfix
      if (replicates > 1) {
        sc <- "_mean"
      } else {
        sc <- ""
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
      if (s == "mean") {
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
  }

  return(self)
}
