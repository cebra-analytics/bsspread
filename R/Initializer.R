#' Initializer class builder
#'
#' Builds a class to initialize a spread simulation, given a specified spatial
#' region and population representation, via either a predefined initial
#' population distribution, or via stochastic generation of invasive species
#' incursions.
#'
#' @param x A \code{raster::RasterLayer}, \code{terra::SpatRaster}, or numeric
#'   vector defining the initial population distribution, or an
#'   \code{Incursions} or inherited class object for generating initial, and
#'   optionally continued, invasive species incursions.
#' @param region A \code{Region} or inherited class object defining the spatial
#'   locations included in the spread simulations.
#' @param population_model A \code{Population} or inherited class object
#'   defining the population representation and growth functionality for the
#'   spread simulations.
#' @param class Character class name for inherited classes. Default is empty.
#' @param ... Additional parameters.
#' @return An \code{Initializer} class object (list) containing functions for
#'   accessing attributes (of the function environment) TODO:
#'   \describe{
#'     \item{\code{TODO()}}{TODO.}
#'   }
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
  Initializer(terra::rast(x), ...)
}

#' @name Initializer
#' @export
Initializer.SpatRaster <- function(x,
                                   region = NULL, ...) {
  if (!is.null(region)) {

    # Check region
    if (!inherits(population_model, "Region")) {
      stop("Region model must be a 'Region' or inherited class object.",
           call. = FALSE)
    }
    if (!region$is_compatible(x)) {
      stop("The spatial object x should be compatible with that defining the ",
           "region.", call. = FALSE)
    }

    # Extract values from locations defined by region
    Incursions(x[region$get_indices()], ...)

  } else { # Use all values
    Incursions(x[], ...)
  }
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
    if (length(x) != region$get_locations()) {
      stop("Vector x length must be equal to the number of region locations.",
           call. = FALSE)
    }
  }

  # Check population model
  if (!is.null(population_model) &&
      !inherits(population_model, "Population")) {
    stop("Population model must be a 'Population' or inherited class object.",
         call. = FALSE)
  }

  # Create a class structure
  self <- structure(list(), class = c(class, "Initializer"))

  # Initialize with x conformed to the population model
  self$initialize <- function() {
    if (!is.null(population_model)) {
      return(population_model$conform(x))
    } else {
      return(x)
    }
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
  if (!inherits(population_model, "Region")) {
    stop("Region model must be a 'Region' or inherited class object.",
         call. = FALSE)
  }

  # Check population model
  if (!is.null(population_model) &&
      !inherits(population_model, "Population")) {
    stop("Population model must be a 'Population' or inherited class object.",
         call. = FALSE)
  }

  # Create a class structure
  self <- structure(list(), class = c(class, "Initializer"))

  # Initialize with generated incursions transformed to the population model
  self$initialize <- function() {
    incursions <- x$generate()
    if (!is.null(population_model)) {
      return(population_model$transform(incursions))
    } else {
      return(incursions)
    }
  }

  return(self)
}
