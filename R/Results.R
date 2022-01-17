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
#' @param ... Additional parameters.
#' @return A \code{Results} class object (list) containing functions for
#'   calculating and collating results, as well as accessing lists of results
#'   and the simulation parameters used to produce them:
#'   \describe{
#'     \item{\code{collate(r, tm, n)}}{Collate the results at simulation
#'       replicate \code{r} and time step \code{tm} using the current
#'       vector/array \code{n} representing the population at each location.}
#'     \item{\code{finalize()}}{Finalize the results collation
#'       (summary calculations).}
#'     \item{\code{as_list()}}{Return the results as a list (collated, total,
#'       area).}
#'     \item{\code{get_params()}}{Get the simulation parameters used.}
#'   }
#' @include Region.R
#' @export
Results <- function(region, population_model,
                    time_steps = 1,
                    step_duration = 1,
                    step_units = "years",
                    collation_steps = 1,
                    replicates = 1, ...) {
  UseMethod("Results")
}

#' @name Results
#' @export
Results.Region <- function(region, population_model,
                           time_steps = 1,
                           step_duration = 1,
                           step_units = "years",
                           collation_steps = 1,
                           replicates = 1, ...) {

  # Validate population model
  if (is.null(population_model) || !inherits(population_model, "Population")) {
    stop("Population model must be a 'Population' or inherited class object.",
         call. = FALSE)
  }

  # Population stages (NULL or number of stages)
  stages <- population_model$get_stages()

  # Initialize result lists
  results <- list(collated = list(), total = list(), area = list())
  zeros <- list(collated = as.integer(population_model$make(initial = 0L)),
                total = 0L, area = 0L)
  if (is.numeric(stages)) {
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

    # Clean n (remove attributes etc.)
    if (population_model$get_type() == "presence_only") {
      n <- as.integer(n)
    }

    # Shape total when population is staged
    if (is.numeric(stages)) {
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

    } else {

      # Population at each location
      if (tm %% collation_steps == 0) {
        results$collated[[tmc]] <<- n
      }

      # Total population size
      results$total[[tmc]] <<- total_n
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
                replicates = replicates))
  }

  return(self)
}
