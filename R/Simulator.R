#' Simulator class builder
#'
#' Builds a class to configure and run replicate discrete-time spread
#' simulations over a given spatial region using (sub)models for population
#' growth and dispersal.
#'
#' @param region A \code{raster::RasterLayer}, \code{terra::SpatRaster}, or
#'   \code{Region} or inherited class object representing the spatial region
#'   (template) for spread simulations.
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
#' @param initializer A \code{Initializer} or inherited class object for
#'   generating the initial invasive species population distribution or
#'   incursion locations, as well as optionally generating subsequent
#'   incursions during the spread simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation and growth functionality for the
#'   spread simulations.
#' @param dispersal_models A list of \code{Dispersal} or inherited class
#'   objects defining the dispersal functionality for the different spread
#'   vectors to be simulated.
#' @param ... Additional parameters.
#' @return A \code{Simulator} class object (list) containing functions for
#'   setting objects (in the function environment) and running the simulations:
#'   \describe{
#'     \item{\code{set_initializer(object)}}{Set the initializer object.}
#'     \item{\code{set_population_model(model)}}{Set the population model
#'       object.}
#'     \item{\code{set_dispersal_models(models)}}{Set the list of dispersal
#'       model objects.}
#'     \item{\code{run()}}{Run the simulations and return the results.}
#'   }
#' @include Region.R
#' @export
Simulator <- function(region,
                      time_steps = 1,
                      step_duration = 1,
                      step_units = "years",
                      collation_steps = 1,
                      replicates = 1,
                      initializer = NULL,
                      population_model = NULL,
                      dispersal_models = list(), ...) {
  UseMethod("Simulator")
}

#' @name Simulator
#' @export
Simulator.Raster <- function(region, ...) {
  # Call Region class version
  Simulator(Region(region), ...)
}

#' @name Simulator
#' @export
Simulator.SpatRaster <- function(region, ...) {
  # Call Region class version
  Simulator(Region(region), ...)
}

#' @name Simulator
#' @export
Simulator.Region <- function(region,
                             time_steps = 1,
                             step_duration = 1,
                             step_units = "years",
                             collation_steps = 1,
                             replicates = 1,
                             initializer = NULL,
                             population_model = NULL,
                             dispersal_models = list(), ...) {

  # Validation
  time_params <- c(time_steps, step_duration, collation_steps, replicates)
  if (!is.numeric(time_params) || !all(time_params > 0)) {
    stop("Time step and replicate parameters should be numeric values > 0.",
         call. = FALSE)
  }
  validate_objects <- function() {
    if (!is.null(initializer) &&
        !inherits(initializer, "Initializer")) {
      initializer <<- NULL
      stop("Initializer must be a 'Initializer' or inherited class object.",
           call. = FALSE)
    }
    if (!is.null(population_model) &&
        !inherits(population_model, "Population")) {
      population_model <<- NULL
      stop("Population model must be a 'Population' or inherited class object.",
           call. = FALSE)
    }
    if (length(dispersal_models) &&
        !all(unlist(lapply(dispersal_models, inherits, "Dispersal")))) {
      dispersal_models <<- list()
      stop("Dispersal models must be 'Dispersal' or inherited class objects.",
           call. = FALSE)
    }
  }
  validate_objects()

  # Create a class structure
  self <- structure(list(), class = "Simulator")

  # Set initializer
  self$set_initializer <- function(object) {
    initializer <<- object
    validate_objects()
  }

  # Set population model
  self$set_population_model <- function(model) {
    population_model <<- model
    validate_objects()
  }

  # Set dispersal models
  self$set_dispersal_models <- function(models) {
    dispersal_models <<- models
    validate_objects()
  }

  # Run simulator
  self$run <- function() {

    # Should at least have an initializer and a population model
    object_absence <- c(is.null(initializer), is.null(population_model))
    if (any(object_absence)) {
      stop(sprintf("The simulator requires the %s be set.",
                 paste(c("initializer", "population model")[object_absence],
                       collapse = " and ")),
           call. = FALSE)
    }

  }

  return(self)
}