#' Initializer class builder
#'
#' Builds a class to initialize a spread simulation, given a specified spatial
#' region and population representation, via either a predefined initial
#' population distribution, or via stochastic generation of invasive species
#' incursions.
#'
#' @param x A \code{raster::Raster*}, \code{terra::SpatRaster}, or numeric
#'   vector or array defining the initial population distribution (where layers
#'   or columns optionally represent stages), or an \code{Incursions} or
#'   inherited class object for generating initial, and optionally continued,
#'   invasive species incursions. The initial population distribution may be
#'   defined with actual populations sizes at appropriate locations (default
#'   with no \code{size} attribute), or via the combination of a binary mask
#'   indicating initial locations and an optional attached \code{size}
#'   attribute indicating the size of the total initial population.
#' @param region A \code{Region} or inherited class object defining the spatial
#'   locations included in the spread simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation for the spread simulations.
#' @param class Character class name for inherited classes. Default is empty.
#' @param ... Additional parameters.
#' @return An \code{Initializer} class object (list) containing a function to
#'   appropriately initialize the population of each location:
#'   \describe{
#'     \item{\code{initialize()}}{Generates an appropriately represented
#'       initial population for each simulated location via either predefined
#'       values, or stochastic incursions.}
#'     \item{\code{continued_incursions()}}{Produces a function for generating
#'       continued incursions when configured to do so. The function produced
#'       has the form: \code{function(tm, n)}, where \code{tm} is the
#'       simulation time step and \code{n} is the array representing the
#'       population at each location at that time. The function returns a
#'       population array with merged incursions.}
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
#' @export
Initializer <- function(x,
                        region = NULL,
                        population_model = NULL,
                        class = character(), ...) {
  UseMethod("Initializer")
}

#' @name Initializer
#' @export
Initializer.Raster <- function(x, ...) {

  # Call the terra version of the function
  x_rast <- terra::rast(x)
  attr(x_rast, "size") <- attr(x, "size")
  Initializer(x_rast, ...)
}

#' @name Initializer
#' @export
Initializer.SpatRaster <- function(x,
                                   region = NULL, ...) {
  if (!is.null(region)) {

    # Check region
    if (!inherits(region, "Region")) {
      stop("Region model must be a 'Region' or inherited class object.",
           call. = FALSE)
    }
    if (!region$is_compatible(x)) {
      stop("The spatial object x should be compatible with that defining the ",
           "region.", call. = FALSE)
    }

    # Extract values from locations defined by region
    x_matrix <- as.matrix(x[region$get_indices()])

  } else { # Use all values
    x_matrix <- as.matrix(x[])
  }
  attr(x_matrix, "size") <- attr(x, "size")
  Initializer(x_matrix, ...)
}

#' @name Initializer
#' @export
Initializer.default <- function(x,
                                region = NULL,
                                population_model = NULL,
                                class = character(), ...) {

  # Check region and x
  if (!is.null(region)) {
    if (!inherits(region, "Region")) {
      stop("Region model must be a 'Region' or inherited class object.",
           call. = FALSE)
    }
    if (nrow(as.matrix(x)) != region$get_locations()) {
      stop("Vector (or array) x length (or number of rows) must be equal to ",
           "the number of region locations.", call. = FALSE)
    }
  }

  # Check population model and x
  if (!is.null(population_model)) {
    if (!inherits(population_model, "Population")) {
      stop("Population model must be a 'Population' or inherited class object.",
           call. = FALSE)
    }
    if (!ncol(as.matrix(x)) %in% c(1, population_model$get_stages())) {
      stop("The number of columns in array x must be consistent with the ",
           "number of stages defined in the population model.", call. = FALSE)
    }
  }

  # Check for population size (thus distribution mask)
  pop_size <- attr(x, "size")

  # Collapse when single column
  if (is.matrix(x) && ncol(x) == 1) {
    x <- x[,1]
  }

  # Convert distribution mask to weightings via capacity/establishment prob.
  if (is.numeric(pop_size) &&
      population_model$get_type() %in% c("unstructured", "stage_structured")) {
    if (is.numeric(population_model$get_capacity())) {
      x <- (x > 0)*population_model$get_capacity()
    } else if (is.numeric(population_model$get_establish_pr())) {
      x <- (x > 0)*population_model$get_establish_pr()
    } else {
      x <- +(x > 0)
    }
  }

  # Create a class structure
  self <- structure(list(), class = c(class, "Initializer"))

  # Initialize using the population model make function
  self$initialize <- function() {
    if (!is.null(population_model)) {
      if (is.numeric(pop_size) &&
          population_model$get_type() %in% c("unstructured",
                                             "stage_structured")) {
        return(population_model$make(
          initial = stats::rmultinom(1, size = pop_size, prob = x)))
      } else {
        return(population_model$make(initial = x))
      }
    } else {
      return(x)
    }
  }

  # Empty continued incursion generator
  self$continued_incursions <- function() {
    return(NULL)
  }

  return(self)
}

#' @name Initializer
#' @export
Initializer.Incursions <- function(x,
                                   region = NULL,
                                   population_model = NULL,
                                   class = character(), ...) {

  # Check region
  if (!inherits(region, "Region")) {
    stop("Region model must be a 'Region' or inherited class object.",
         call. = FALSE)
  }

  # Check population model
  if (!is.null(population_model) &&
      !inherits(population_model, "Population")) {
    stop("Population model must be a 'Population' or inherited class object.",
         call. = FALSE)
  }

  # Set incursion population mean and stages where applicable
  if (!is.null(x$get_incursion_mean())) {
    population_model$set_incursion_mean(x$get_incursion_mean())
  }
  if (!is.null(x$get_incursion_stages()) &&
      is.function(population_model$set_incursion_stages)) {
    population_model$set_incursion_stages(x$get_incursion_stages())
  }

  # Create a class structure
  self <- structure(list(), class = c(class, "Initializer"))

  # Initialize with generated incursions transformed to the population model
  self$initialize <- function() {
    x_i <- x$generate()
    if (!is.null(population_model)) {
      return(population_model$make(incursion = x_i))
    } else {
      return(x_i)
    }
  }

  # Continued incursion generator
  self$continued_incursions <- function() {
    if (x$get_continued()) {
      return(
        function(tm, n) {
          if (tm > 1 && tm %% x$get_time_steps() == 0) {
            x_i <- x$generate()
            if (!is.null(population_model)) {
              return(population_model$make(current = n, incursion = x_i))
            } else {
              return(as.logical(n) | x_i)
            }
          } else {
            return(n)
          }
        }
      )
    } else {
      return(NULL)
    }
  }

  return(self)
}
