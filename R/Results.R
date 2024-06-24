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
#'       may be passed to the function.}
#'     \item{\code{save_csv()}}{Save the collated results as comma-separated
#'       values (CSV) files when the region is patch-based. Also saves the
#'       population totals and area occupied to CSV files for both grid and
#'       patch-based region types.}
#'     \item{\code{save_plots()}}{Save plots of the population (staged) totals
#'       and the area occupied as PNG files.}
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

  # Include occupancy?
  include_occupancy <- ((population_model$get_type() %in%
                          c("unstructured", "stage_structured")) &&
                          !region$spatially_implicit())

  # Initialize result lists
  results <- list(collated = list(), total = list(), area = list())
  if (region$get_type() == "grid") {
    cell_areas <-
      terra::cellSize(region$get_template())[region$get_indices()][,1]
    attr(results$area, "units") <- "square metres"
  } else if (region$spatially_implicit()) {
    attr(results$area, "units") <- "square metres"
  } else {
    attr(results$area, "units") <- "patches"
  }
  if (is.numeric(stages) && is.null(combine_stages)) {
    zeros <- list(collated = population_model$make(initial = 0L))
    zeros$total <- zeros$collated[1,, drop = FALSE]
  } else if (is.numeric(combine_stages)) {
    zeros <- list(collated = array(0L, c(region$get_locations(), 1)))
    colnames(zeros$collated) <- stage_labels
    zeros$total <- zeros$collated[1,, drop = FALSE]
  } else {
    zeros <- list(collated = rep(0L, region$get_locations()), total = 0L)
  }
  zeros$area <- 0L
  if (include_occupancy) {
    results$occupancy <- list()
    zeros$occupancy <- rep(0L, region$get_locations())
    results$total_occup <- list()
    zeros$total_occup <- 0L
  }
  if (replicates > 1) { # summaries
    zeros$collated <- list(mean = zeros$collated, sd = zeros$collated)
    zeros$total <- list(mean = zeros$total, sd = zeros$total)
    zeros$area <- list(mean = zeros$area, sd = zeros$area)
    if (include_occupancy) {
      zeros$occupancy <- list(mean = zeros$occupancy, sd = zeros$occupancy)
      zeros$total_occup <- list(mean = zeros$total_occup,
                                sd = zeros$total_occup)
    }
  }
  for (tm in as.character(c(0, seq(collation_steps, time_steps,
                                   by = collation_steps)))) {
    results$collated[[tm]] <- zeros$collated
    if (include_occupancy) {
      results$occupancy[[tm]] <- zeros$occupancy
    }
  }
  for (tm in as.character(0:time_steps)) {
    results$total[[tm]] <- zeros$total
    results$area[[tm]] <- zeros$area
    if (include_occupancy) {
      results$total_occup[[tm]] <- zeros$total_occup
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
    if (include_occupancy) {
      occupancy <- +(rowSums(as.matrix(n)) > 0)
      total_occup <- sum(occupancy)
    }

    # Shape total when population is staged
    if (is.numeric(stages)) {
      total_n <- array(colSums(n), c(1, ncol(n)))
      colnames(total_n) <- stage_labels
    } else {
      total_n <- sum(n)
    }

    if (replicates > 1) { # summaries

      # Calculates running mean and standard deviation (note: variance*r is
      # stored as SD and transformed at the final replicate and time step)

      # Population summaries at each location
      if (tm %% collation_steps == 0) {
        previous_mean <- results$collated[[tmc]]$mean
        results$collated[[tmc]]$mean <<- previous_mean + (n - previous_mean)/r
        previous_sd <- results$collated[[tmc]]$sd
        results$collated[[tmc]]$sd <<-
          (previous_sd + ((n - previous_mean)*
                            (n - results$collated[[tmc]]$mean)))
      }

      # Total population size summaries
      previous_mean <- results$total[[tmc]]$mean
      results$total[[tmc]]$mean <<- previous_mean + (total_n - previous_mean)/r
      previous_sd <- results$total[[tmc]]$sd
      results$total[[tmc]]$sd <<-
        (previous_sd + ((total_n - previous_mean)*
                          (total_n - results$total[[tmc]]$mean)))

      # Total area occupied summaries
      previous_mean <- results$area[[tmc]]$mean
      results$area[[tmc]]$mean <<- (previous_mean +
                                      (total_area - previous_mean)/r)
      previous_sd <- results$area[[tmc]]$sd
      results$area[[tmc]]$sd <<-
        (previous_sd + ((total_area - previous_mean)*
                          (total_area - results$area[[tmc]]$mean)))

      # Occupancy summaries at each location
      if (include_occupancy && (tm %% collation_steps == 0)) {
        previous_mean <- results$occupancy[[tmc]]$mean
        results$occupancy[[tmc]]$mean <<-
          previous_mean + (occupancy - previous_mean)/r
        previous_sd <- results$occupancy[[tmc]]$sd
        results$occupancy[[tmc]]$sd <<-
          (previous_sd + ((occupancy - previous_mean)*
                            (occupancy - results$occupancy[[tmc]]$mean)))
      }

      # Total occupancy summaries
      if (include_occupancy) {
        previous_mean <- results$total_occup[[tmc]]$mean
        results$total_occup[[tmc]]$mean <<- (previous_mean +
                                               (total_occup - previous_mean)/r)
        previous_sd <- results$total_occup[[tmc]]$sd
        results$total_occup[[tmc]]$sd <<-
          (previous_sd + ((total_occup - previous_mean)*
                            (total_occup - results$total_occup[[tmc]]$mean)))
      }

    } else {

      # Population at each location
      if (tm %% collation_steps == 0) {
        results$collated[[tmc]] <<- n
      }

      # Total population size
      results$total[[tmc]] <<- total_n

      # Total area occupied
      results$area[[tmc]] <<- total_area

      # Occupancy at each location
      if (include_occupancy && (tm %% collation_steps == 0)) {
        results$occupancy[[tmc]] <<- occupancy
      }

      # Total occupancy
      if (include_occupancy) {
        results$total_occup[[tmc]] <<- total_occup
      }
    }
  }

  # Finalize the results collation
  self$finalize <- function() {

    if (replicates > 1) { # summaries

      # Transform collated population standard deviations
      for (tmc in names(results$collated)) {
        results$collated[[tmc]]$sd <<-
          sqrt(results$collated[[tmc]]$sd/(replicates - 1))
      }

      # Transform total population standard deviations
      for (tmc in names(results$total)) {
        results$total[[tmc]]$sd <<-
          sqrt(results$total[[tmc]]$sd/(replicates - 1))
      }

      # Transform area occupied standard deviations
      for (tmc in names(results$area)) {
        results$area[[tmc]]$sd <<-
          sqrt(results$area[[tmc]]$sd/(replicates - 1))
      }

      # Transform collated occupancy standard deviations
      if (include_occupancy) {
        for (tmc in names(results$occupancy)) {
          results$occupancy[[tmc]]$sd <<-
            sqrt(results$occupancy[[tmc]]$sd/(replicates - 1))
        }
      }

      # Transform collated total occupancy standard deviations
      if (include_occupancy) {
        for (tmc in names(results$total_occup)) {
          results$total_occup[[tmc]]$sd <<-
            sqrt(results$total_occup[[tmc]]$sd/(replicates - 1))
        }
      }
    }

    # Add labels to staged populations (again)
    if (population_model$get_type() == "stage_structured") {
      if (replicates > 1) {
        summaries <- c("mean", "sd")
      } else {
        summaries <- ""
      }
      for (tmc in names(results$collated)) {
        for (s in summaries) {
          if (replicates > 1) {
            colnames(results$collated[[tmc]][[s]]) <<- stage_labels
            colnames(results$total[[tmc]][[s]]) <<- stage_labels
          } else {
            colnames(results$collated[[tmc]]) <<- stage_labels
            colnames(results$total[[tmc]]) <<- stage_labels
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

  # Save collated results as raster files
  if (region$get_type() == "grid") {
    self$save_rasters  <- function(...) {

      # Replicate summaries or single replicate
      if (replicates > 1) {
        summaries <- c("mean", "sd")
      } else {
        summaries <- ""
      }

      # Label population appropriately
      if (population_model$get_type() == "presence_only") {
        pop_label <- "occupancy"
      } else {
        pop_label <- "population"
      }

      # Save rasters for each time step
      for (tmc in names(results$collated)) {
        for (s in summaries) {

          # Create and save a results raster per stage
          if (is.null(stages) || is.numeric(combine_stages)) {
            stages <- 1
          }
          for (i in 1:stages) {

            # Create result raster
            if (population_model$get_type() == "stage_structured") {
              if (replicates > 1) {
                output_rast <-
                  region$get_rast(results$collated[[tmc]][[s]][,i])
              } else {
                output_rast <- region$get_rast(results$collated[[tmc]][,i])
              }
              if (is.null(combine_stages)) {
                names(output_rast) <- stage_labels[i]
                i <- paste0("_stage_", i)
              } else if (is.numeric(combine_stages)) {
                names(output_rast) <- "combined"
                i <- ""
              }
            } else {
              if (replicates > 1) {
                output_rast <- region$get_rast(results$collated[[tmc]][[s]])
              } else {
                output_rast <- region$get_rast(results$collated[[tmc]])
              }
              i <- ""
            }

            # Write raster to file
            if (replicates > 1) {
              sc <- paste0("_", s)
            } else {
              sc <- s
            }
            filename <- sprintf(paste0(pop_label, "%s_t%0",
                                       nchar(as.character(time_steps)),
                                       "d%s.tif"), i, as.integer(tmc), sc)
            terra::writeRaster(output_rast, filename, ...)
          }
        }
      }

      # Save occupancy rasters for each time step
      if (include_occupancy) {
        for (tmc in names(results$occupancy)) {
          for (s in summaries) {

            # Copy results into a raster
            if (replicates > 1) {
              output_rast <- region$get_rast(results$occupancy[[tmc]][[s]])
              s <- paste0("_", s)
            } else {
              output_rast <- region$get_rast(results$occupancy[[tmc]])
            }

            # Write raster to file
            filename <- sprintf(paste0("occupancy_t%0",
                                       nchar(as.character(time_steps)),
                                       "d%s.tif"), as.integer(tmc), s)
            terra::writeRaster(output_rast, filename, ...)
          }
        }
      }
    }
  }

  # Save collated (patch only) and summary (both) results as CSV files
  self$save_csv  <- function() {

    # Replicate summaries or single replicate
    if (replicates > 1) {
      summaries <- c("mean", "sd")
    } else {
      summaries <- 1
    }
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

    # Label population appropriately
    if (population_model$get_type() == "presence_only") {
      pop_label <- "occupancy"
    } else {
      pop_label <- "population"
    }

    # Collated results for patch only
    if (region$spatially_implicit()) {

      # Collect values
      output_df <- list()
      for (s in summaries) {
        if (replicates > 1) {
          output_df[[s]] <- array(sapply(results$collated,
                                   function(c_tm) c_tm[[s]]),
                                  c(result_stages, length(results$collated)))
        } else {
          output_df[[s]] <- array(sapply(results$collated,
                                         function(c_tm) c_tm),
                                  c(result_stages, length(results$collated)))
        }
        colnames(output_df[[s]]) <- collated_labels
        if (population_model$get_type() == "stage_structured") {
          rownames(output_df[[s]]) <- stage_labels
        } else {
          rownames(output_df[[s]]) <- pop_label
        }
      }

      # Write to CSV files
      for (s in summaries) {
        filename <- sprintf(paste0(pop_label, "%s.csv"), s_fname[[s]])
        utils::write.csv(output_df[[s]], filename, row.names = TRUE)
      }

    } else if (region$get_type() == "patch") {

      # Location coordinates and labels
      coords <- region$get_coords(extra_cols = TRUE)
      coords <- coords[, c("lon", "lat",
                           names(which(sapply(coords, is.character))))]

      # Save population CSV file(s)
      for (i in 1:result_stages) {

        # Combine coordinates and population values
        output_df <- list()
        for (s in summaries) {
          if (population_model$get_type() == "stage_structured") {
            if (replicates > 1) {
              output_df[[s]] <- lapply(results$collated,
                                       function(c_tm) c_tm[[s]][,i])
            } else {
              output_df[[s]] <- lapply(results$collated,
                                       function(c_tm) c_tm[,i])
            }
          } else {
            if (replicates > 1) {
              output_df[[s]] <- lapply(results$collated,
                                       function(c_tm) c_tm[[s]])
            } else {
              output_df[[s]] <- results$collated
            }
          }
          names(output_df[[s]]) <- collated_labels
          output_df[[s]] <- cbind(coords, as.data.frame(output_df[[s]]))
        }

        # Write to CSV files
        for (s in summaries) {
          filename <- sprintf(paste0(pop_label, "%s%s.csv"), i_fname[i],
                              s_fname[[s]])
          utils::write.csv(output_df[[s]], filename, row.names = FALSE)
        }
      }

      # Save occupancy CSV file(s)
      if (include_occupancy) {

        # Combine coordinates and occupancy values
        output_df <- list()
        for (s in summaries) {
          if (replicates > 1) {
            output_df[[s]] <- lapply(results$occupancy,
                                     function(o_tm) o_tm[[s]])
          } else {
            output_df[[s]] <- results$occupancy
          }
          names(output_df[[s]]) <- collated_labels
          output_df[[s]] <- cbind(coords, as.data.frame(output_df[[s]]))
        }

        # Write to CSV files
        for (s in summaries) {
          filename <- sprintf("occupancy%s.csv", s_fname[[s]])
          utils::write.csv(output_df[[s]], filename, row.names = FALSE)
        }
      }
    }

    # Population totals and area occupied
    for (s in summaries) {

      # Collect total population, area occupied, and occupancy
      if (replicates > 1) {
        totals <- array(sapply(results$total, function(tot) tot[[s]]),
                        c(result_stages, time_steps + 1))
        areas <- array(sapply(results$area, function(area) area[[s]]),
                       c(1, time_steps + 1))
        if (include_occupancy) {
          total_occup <- array(sapply(results$total_occup,
                                      function(occup) occup[[s]]),
                         c(1, time_steps + 1))
        }
      } else {
        totals <- array(sapply(results$total, function(tot) tot),
                        c(result_stages, time_steps + 1))
        areas <- array(sapply(results$area, function(area) area),
                       c(1, time_steps + 1))
        if (include_occupancy) {
          total_occup <- array(sapply(results$total_occup,
                                      function(occup) occup),
                         c(1, time_steps + 1))
        }
      }

      # Label columns and rows
      colnames(totals) <- time_steps_labels
      colnames(areas) <- time_steps_labels
      if (population_model$get_type() == "stage_structured") {
        rownames(totals) <- stage_labels
      } else {
        rownames(totals) <- pop_label
      }
      rownames(areas) <- attr(results$area, "units")
      if (include_occupancy) {
        colnames(total_occup) <- time_steps_labels
        rownames(total_occup) <- "occupancy"
      }

      # Write to CSV files
      if (!region$spatially_implicit()) {
        filename <- sprintf(paste0("total_", pop_label, "%s.csv"),
                                   s_fname[[s]])
        utils::write.csv(totals, filename, row.names = TRUE)
      }
      if (region$spatially_implicit()) {
        filename <- sprintf("area_occupied%s.csv", s_fname[[s]])
      } else {
        filename <- sprintf("total_area_occupied%s.csv", s_fname[[s]])
      }
      utils::write.csv(areas, filename, row.names = TRUE)
      if (include_occupancy) {
        filename <- sprintf("total_occupancy%s.csv", s_fname[[s]])
        utils::write.csv(total_occup, filename, row.names = TRUE)
      }
    }
  }

  # Plot total population (per stage) and area occupied as PNG files
  self$save_plots  <- function() {

    # Resolve the number of (combined) stages used in the results
    result_stages <- stages
    if (is.null(stages) || is.numeric(combine_stages)) {
      result_stages <- 1
    }

    # Stage label for plot headings
    stage_label <- ""
    stage_file <- ""
    if (population_model$get_type() == "stage_structured") {
      if (is.numeric(stages) && is.null(combine_stages)) {
        stage_label <- paste0(stage_labels, " ")
        stage_file <- paste0("_stage_", 1:result_stages)
      } else if (is.numeric(combine_stages)) {
        stage_label <- "combined "
      }
    }

    # Label population and totals appropriately
    if (population_model$get_type() == "presence_only") {
      pop_label <- "occupancy"
    } else {
      pop_label <- "population"
    }
    if (region$spatially_implicit()) {
      tot_label <- list(text = "", file = "")
    } else {
      tot_label <- list(text = "total ", file = "total_")
    }

    # All plots have time steps on x-axis
    plot_x_label <- paste0("Time steps (", step_units, ")")

    if (replicates > 1) { # plot summary mean +/- 2 SD

      # Collect totals and areas
      totals <- list()
      areas <- list()
      for (s in c("mean", "sd")) {
        totals[[s]] <- array(sapply(results$total, function(tot) tot[[s]]),
                             c(result_stages, time_steps + 1))
        areas[[s]] <- sapply(results$area, function(area) area[[s]])
      }

      # Collect total occupancy
      if (include_occupancy) {
        occup <- list()
        for (s in c("mean", "sd")) {
          occup[[s]] <- sapply(results$total_occup, function(occup) occup[[s]])
        }
      }

      # Plot totals (per result stage)
      for (s in 1:result_stages) {
        filename <- paste0(tot_label$file, pop_label, stage_file[s], ".png")
        main_title <- paste0(tot_label$text, stage_label[s], pop_label,
                             " (mean +/- 2 SD)")
        main_title <- paste0(toupper(substr(main_title, 1, 1)),
                             substr(main_title, 2, nchar(main_title)))
        y_label <- paste0(toupper(substr(pop_label, 1, 1)),
                             substr(pop_label, 2, nchar(pop_label)))
        grDevices::png(filename = filename)
        graphics::plot(0:time_steps, totals$mean[s,], type = "l",
                       main = main_title,
                       xlab = plot_x_label, ylab = y_label,
                       ylim = c(0, 1.1*max(totals$mean[s,] +
                                             2*totals$sd[s,])))
        graphics::lines(0:time_steps, totals$mean[s,] + 2*totals$sd[s,],
                        lty = "dashed")
        graphics::lines(0:time_steps,
                        pmax(0, totals$mean[s,] - 2*totals$sd[s,]),
                        lty = "dashed")
        invisible(grDevices::dev.off())
      }

      # Plot areas
      if (region$spatially_implicit()) {
        filename <- "area_occupied.png"
        main_title <- "Area occupied (mean +/- 2 SD)"
      } else {
        filename <- "total_area_occupied.png"
        main_title <- "Total area occupied (mean +/- 2 SD)"
      }
      grDevices::png(filename = filename)
      graphics::plot(0:time_steps, areas$mean, type = "l",
                     main = main_title,
                     xlab = plot_x_label,
                     ylab = paste0("Area (", attr(results$area, "units"),
                                   ")"),
                     ylim = c(0, 1.1*max(areas$mean + 2*areas$sd)))
      graphics::lines(0:time_steps, areas$mean + 2*areas$sd, lty = "dashed")
      graphics::lines(0:time_steps, pmax(0, areas$mean - 2*areas$sd),
                      lty = "dashed")
      invisible(grDevices::dev.off())

      # Plot occupancy
      if (include_occupancy) {
        grDevices::png(filename = "total_occupancy.png")
        graphics::plot(0:time_steps, occup$mean, type = "l",
                       main = "Total occupancy (mean +/- 2 SD)",
                       xlab = plot_x_label,
                       ylab = "Occupancy",
                       ylim = c(0, 1.1*max(occup$mean + 2*occup$sd)))
        graphics::lines(0:time_steps, occup$mean + 2*occup$sd, lty = "dashed")
        graphics::lines(0:time_steps, pmax(0, occup$mean - 2*occup$sd),
                        lty = "dashed")
        invisible(grDevices::dev.off())
      }

    } else { # plot values

      # Collect totals and areas
      totals <- array(sapply(results$total, function(tot) tot),
                      c(result_stages, time_steps + 1))
      areas <- sapply(results$area, function(area) area)

      # Collect total occupancy
      if (include_occupancy) {
        occup <- sapply(results$total_occup, function(occup) occup)
      }

      # Plot totals (per result stage)
      for (s in 1:result_stages) {
        filename <- paste0(tot_label$file, pop_label, stage_file[s], ".png")
        main_title <- paste0(tot_label$text, stage_label[s], pop_label)
        main_title <- paste0(toupper(substr(main_title, 1, 1)),
                             substr(main_title, 2, nchar(main_title)))
        y_label <- paste0(toupper(substr(pop_label, 1, 1)),
                          substr(pop_label, 2, nchar(pop_label)))
        grDevices::png(filename = filename)
        graphics::plot(0:time_steps, totals[s,], type = "l",
                       main = main_title,
                       xlab = plot_x_label, ylab = y_label,
                       ylim = c(0, 1.1*max(totals[s,])))
        invisible(grDevices::dev.off())
      }

      # Plot areas
      if (region$spatially_implicit()) {
        filename <- "area_occupied.png"
        main_title <- "Area occupied"
      } else {
        filename <- "total_area_occupied.png"
        main_title <- "Total area occupied"
      }
      grDevices::png(filename = filename)
      graphics::plot(0:time_steps, areas, type = "l",
                     main = main_title,
                     xlab = plot_x_label,
                     ylab = paste0("Area (", attr(results$area, "units"),
                                   ")"),
                     ylim = c(0, 1.1*max(areas)))
      invisible(grDevices::dev.off())

      # Plot total occupancy
      if (include_occupancy) {
        grDevices::png(filename = "total_occupancy.png")
        graphics::plot(0:time_steps, occup, type = "l",
                       main = "Total occupancy",
                       xlab = plot_x_label,
                       ylab = "Occupancy",
                       ylim = c(0, 1.1*max(occup)))
        invisible(grDevices::dev.off())
      }
    }
  }

  return(self)
}
