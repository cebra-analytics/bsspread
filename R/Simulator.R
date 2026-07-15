#' Simulator class builder
#'
#' Builds a class to configure and run replicate discrete-time incursion
#' management simulations over a given spatial region using (sub)models for
#' simulating population growth and dispersal, calculating impacts, and
#' applying management actions, such as detection, control and removal.
#'
#' @param region A \code{raster::RasterLayer}, \code{terra::SpatRaster}, or
#'   \code{Region} or inherited class object representing the spatial region
#'   (template) for the spread simulations.
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
#' @param result_stages Optionally combine (sum) specified stages (a vector of
#'   stage indices) of stage-based population results. The default
#'   (\code{NULL}) maintains results for each stage.
#' @param parallel_cores Number of cores available for parallel processing.
#'   The default NULL implies no parallel processing.
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
#' @param impacts A list of \code{Impacts} class objects specifying how to
#'   calculate various impacts of the simulated population at each time step.
#' @param actions A list of \code{Actions} or inherited class objects for
#'   applying invasive species management actions, such as detection, control,
#'   and removal.
#' @param user_function An optional user-defined \code{function(n, r, tm)} that
#'   is applied to the population vector or matrix \code{n} (returning a
#'   transformed \code{n}) prior to collating the results at simulation
#'   replicate \code{r} and time step \code{tm}.
#' @param class Character class name for inherited classes. Default is
#'   \code{NULL}.
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
#' @include Results.R
#' @export
Simulator <- function(region,
                      time_steps = 1,
                      step_duration = 1,
                      step_units = "years",
                      collation_steps = 1,
                      replicates = 1,
                      result_stages = NULL,
                      parallel_cores = NULL,
                      initializer = NULL,
                      population_model = NULL,
                      dispersal_models = list(),
                      impacts = list(),
                      actions = list(),
                      user_function = NULL,
                      class = character(), ...) {
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
                             result_stages = NULL,
                             parallel_cores = NULL,
                             initializer = NULL,
                             population_model = NULL,
                             dispersal_models = list(),
                             impacts = list(),
                             actions = list(),
                             user_function = NULL,
                             class = character(), ...) {

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
    if (!is.null(user_function) && !is.function(user_function)) {
      stop("User-defined function must define a function for transforming ",
           "the population at each simulation time step.", call. = FALSE)
    }
  }
  validate_objects()
  if (!is.null(result_stages) && !is.null(population_model) &&
      population_model$get_type() == "stage_structured") {
    if (!is.numeric(result_stages) ||
        !all(result_stages %in% 1:population_model$get_stages())) {
      stop("Result stages must be a numeric vector of stage indices ",
           "consistent with that defined in the population model.",
           call. = FALSE)
    }
  }
  # Check parallel cores
  if (!is.null(parallel_cores) && (!is.numeric(parallel_cores) ||
                                   parallel_cores <= 0)) {
    stop("The number of parallel cores should be a numeric value > 0.",
         call. = FALSE)
  }

  # Check impact objects and set impact ids
  if (length(impacts)) {
    if (!all(unlist(lapply(impacts, inherits, "Impacts")))) {
      stop("Impacts must be a list of 'Impacts' objects.", call. = FALSE)
    }
    for (i in 1:length(impacts)) {
      impacts[[i]]$set_id(i)
    }
  }

  # Check actions objects and set action ids
  if (length(actions)) {
    if (!all(unlist(lapply(actions, inherits, "Actions")))) {
      stop("Actions must be a list of 'Actions' objects.", call. = FALSE)
    }
    for (i in 1:length(actions)) {
      actions[[i]]$set_id(i)
    }
  }

  # Set parallel cores in region and dispersal objects
  set_cores <- function(cores = NULL) {
    region$set_cores(cores)
    if (length(dispersal_models)) {
      for (i in 1:length(dispersal_models)) {
        dispersal_models[[i]]$set_cores(cores)
      }
    }
  }
  set_cores(cores = parallel_cores)

  # Create a class structure
  self <- structure(list(), class = c(class, "Simulator"))

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
    set_cores(cores = parallel_cores)
  }

  # Run simulator
  results <- NULL # DEBUG ####
  self$run <- function() {

    # Should at least have an initializer and a population model
    object_absence <- c(is.null(initializer), is.null(population_model))
    if (any(object_absence)) {
      stop(sprintf("The simulator requires the %s to be set.",
                   paste(c("initializer", "population model")[object_absence],
                         collapse = " and ")),
           call. = FALSE)
    }

    # Continued incursions function
    continued_incursions <- initializer$continued_incursions()

    # Results setup
    results <<- Results(region, population_model, # DEBUG ####
                        impacts = impacts,
                        actions = actions,
                        time_steps = time_steps,
                        step_duration = step_duration,
                        step_units = step_units,
                        collation_steps = collation_steps,
                        replicates = replicates,
                        combine_stages = result_stages)

    # Replicates
    for (r in 1:replicates) {

      # Initialize population array
      n <- initializer$initialize()

      # Set diffusion attributes when spatially implicit (single patch)
      if (region$spatially_implicit()) {

        # Diffusion model
        if (any(sapply(dispersal_models,
                       function(dm) inherits(dm, "Diffusion")))) {
          idx <- which(sapply(dispersal_models,
                              function(dm) inherits(dm, "Diffusion")))[1]
          attr(n, "initial_n") <- n
          attr(n, "diffusion_rate") <-
            dispersal_models[[idx]]$get_diffusion_rate()
          attr(n, "diffusion_radius") <- 0
        }

        # Area spread model
        if (any(sapply(dispersal_models,
                       function(dm) inherits(dm, "AreaSpread")))) {
          capacity <- population_model$get_capacity()
          capacity_area <- attr(capacity, "area")
          if (population_model$get_type() == "stage_structured") {
            stages <- population_model$get_capacity_stages()
            attr(n, "spread_area") <-
              sum(n[,stages])/as.numeric(capacity)*capacity_area
          } else { # unstructured
            attr(n, "spread_area") <- n/as.numeric(capacity)*capacity_area
          }
        }
      }

      # Calculate impacts
      if (length(impacts)) {

        # Calculate each impact and update recovery delays where applicable
        for (i in 1:length(impacts)) {
          n <- impacts[[i]]$calculate(n, 0)
          n <- impacts[[i]]$update_recovery_delay(n, 0)
        }

        # Apply any dynamically linked impacts to capacity
        population_model$set_capacity_mult(n)
      }

      # Apply actions
      if (length(actions)) {
        for (i in 1:length(actions)) {
          n <- actions[[i]]$apply(n, 0)
        }
      }

      # Initial results (t = 0)
      results$collate(r, 0, n)

      # Time steps
      for (tm in 1:time_steps) {

        # Population growth
        n <- population_model$grow(n, tm)

        # Dispersal for each spread vector
        if (length(dispersal_models)) {

          # Pack into list of original, remaining and relocated populations
          n <- dispersal_models[[1]]$pack(n)

          # Perform dispersal for each spread vector
          for (i in 1:length(dispersal_models)) {
            n <- dispersal_models[[i]]$disperse(n, tm)
          }

          # Unpack population array from separated list
          n <- dispersal_models[[1]]$unpack(n)
        }

        # Calculate impacts
        if (length(impacts)) {

          # Calculate each impact and update recovery delays where applicable
          for (i in 1:length(impacts)) {
            n <- impacts[[i]]$calculate(n, tm)
            n <- impacts[[i]]$update_recovery_delay(n, tm)
          }

          # Apply any dynamically linked impacts to capacity
          population_model$set_capacity_mult(n)
        }

        # Apply actions
        if (length(actions)) {

          # Clear attributes
          for (i in 1:length(actions)) {
            n <- actions[[i]]$clear_attributes(n)
          }

          # Apply sequentially
          for (i in 1:length(actions)) {
            n <- actions[[i]]$apply(n, tm)
          }
        }

        # User-defined function
        if (is.function(user_function)) {
          if (length(formals(user_function)) == 3) {
            n <- user_function(n, r, tm)
          } else { # previously just n
            n <- user_function(n)
          }
        }

        # Collate results
        results$collate(r, tm, n)

        # Continued incursions
        if (is.function(continued_incursions)) {
          n <- continued_incursions(tm, n)
        }

      } # time steps

    } # replicates

    # Finalize results
    results$finalize()

    # Return results
    return(results)
  }

  return(self)
}
