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
#'   }
#' @include Region.R
#' @export
Results <- function(region, population_model,
                    time_steps = 1,
                    step_duration = 1,
                    step_units = "years",
                    collation_steps = 1,
                    replicates = 1,
                    combine_stages = NULL, ...) {
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
                           combine_stages = NULL, ...) {

  # Validate population model
  if (is.null(population_model) || !inherits(population_model, "Population")) {
    stop("Population model must be a 'Population' or inherited class object.",
         call. = FALSE)
  }

  # Population stages (NULL or number of stages)
  stages <- population_model$get_stages()

  # Initialize result lists
  results <- list(collated = list(), total = list(), area = list())
  if (region$get_type() == "grid") {
    cell_areas <-
      terra::cellSize(region$get_template())[region$get_indices()][,1]
    attr(results$area, "units") <- "square meters"
  } else if (region$spatially_implicit()) {
    attr(results$area, "units") <- "square meters"
  } else {
    attr(results$area, "units") <- "patches"
  }
  zeros <- list(collated = as.integer(population_model$make(initial = 0L)),
                total = 0L, area = 0L)
  if (is.numeric(stages) && is.null(combine_stages)) {
    zeros$collated <- matrix(zeros$collated, ncol = stages)
    zeros$total <- array(0L, c(1, stages))
  }
  if (replicates > 1) { # summaries
    zeros$collated <- list(mean = zeros$collated, sd = zeros$collated)
    zeros$total <- list(mean = zeros$total, sd = zeros$total)
    zeros$area <- list(mean = zeros$area, sd = zeros$area)
  }
  for (tm in as.character(c(0, seq(collation_steps, time_steps,
                                   by = collation_steps)))) {
    results$collated[[tm]] <- zeros$collated
  }
  for (tm in as.character(0:time_steps)) {
    results$total[[tm]] <- zeros$total
    results$area[[tm]] <- zeros$area
  }
  rm(zeros)

  # Create a class structure
  self <- structure(list(), class = "Results")

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
      n <- rowSums(n[,combine_stages, drop = FALSE])
    }

    # Shape total when population is staged
    if (is.numeric(stages) && is.null(combine_stages)) {
      total_n <- array(colSums(n), c(1, stages))
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

    } else {

      # Population at each location
      if (tm %% collation_steps == 0) {
        results$collated[[tmc]] <<- n
      }

      # Total population size
      results$total[[tmc]] <<- total_n

      # Total area occupied
      results$area[[tmc]] <<- total_area
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

      # Save rasters for each time step
      for (tmc in names(results$collated)) {
        for (s in summaries) {

          # Create a raster with a layer per stage
          if (is.null(stages) || is.numeric(combine_stages)) {
            stages <- 1
          }
          output_rast <- list()
          for (i in 1:stages) {
            output_rast[[i]] <- region$get_template()
          }
          output_rast <- terra::rast(output_rast)

          # Copy results into the raster
          if (replicates > 1) {
            output_rast[region$get_indices()] <- results$collated[[tmc]][[s]]
            s <- paste0("_", s)
          } else {
            output_rast[region$get_indices()] <- results$collated[[tmc]]
          }

          # Write raster to file
          filename <- sprintf(paste0("result_t%0",
                                     nchar(as.character(time_steps)),
                                     "d%s.tif"), as.integer(tmc), s)
          terra::writeRaster(output_rast, filename, ...)
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
      summaries <- ""
    }

    # Resolve stages
    if (is.null(stages) || is.numeric(combine_stages)) {
      stages <- 1
    }

    # Collated results for patch only
    if (region$get_type() == "patch") {

      # Location coordinates and labels
      coords <- region$get_coords(extra_cols = TRUE)
      coords <- coords[, c("lon", "lat",
                           names(which(sapply(coords, is.character))))]

      # Save CSV for each time step
      for (tmc in names(results$collated)) {
        for (s in summaries) {

          # Combine coordinates, labels, and results
          if (replicates > 1) {
            collated_results <- cbind(coords, results$collated[[tmc]][[s]])
            s <- paste0("_", s)
          } else {
            collated_results <- cbind(coords, results$collated[[tmc]])
          }

          # Write to CSV file
          filename <- sprintf(paste0("result_t%0",
                                     nchar(as.character(time_steps)),
                                     "d%s.csv"), as.integer(tmc), s)
          utils::write.csv(collated_results, filename, row.names = FALSE)
        }
      }
    }

    # Population totals and area occupied
    for (s in summaries) {

      # Collect totals and area occupied
      if (replicates > 1) {
        totals <- array(sapply(results$total, function(tot) tot[[s]]),
                        c(stages, time_steps + 1))
        areas <- array(sapply(results$area, function(area) area[[s]]),
                       c(1, time_steps + 1))
        s <- paste0("_", s)
      } else {
        totals <- array(sapply(results$total, function(tot) tot),
                        c(stages, time_steps + 1))
        areas <- array(sapply(results$area, function(area) area),
                       c(1, time_steps + 1))
      }
      time_steps_labels <- sprintf(paste0("t%0",
                                          nchar(as.character(time_steps)),
                                          "d"), as.integer(0:time_steps))
      colnames(totals) <- time_steps_labels
      colnames(areas) <- time_steps_labels

      # Write to CSV files
      filename <- sprintf("result_totals%s.csv", s)
      utils::write.csv(totals, filename, row.names = FALSE)
      filename <- sprintf("result_areas%s.csv", s)
      utils::write.csv(areas, filename, row.names = FALSE)
    }
  }

  return(self)
}
