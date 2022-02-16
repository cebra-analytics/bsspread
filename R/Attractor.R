#' Attractor class builder
#'
#' Builds a class to represent spatially distributed weights or scaling factors
#' assigned to cells/locations to represent an "attractor" that modify the
#' likelihood of dispersal from source or to destination locations. Attractor
#' values may be used to scale the proportion of each population dispersing
#' from each source location, or to scale the number of dispersal events
#' originating at each source location. Alternatively, attractor weight values
#' may be used to modify the relative likelihood of dispersal to each
#' destination location. The \code{type} attribute of an attractor may thus be
#' defined as "source", "destination", or "both".
#'
#' @param x A \code{raster::RasterLayer} or \code{terra::SpatRaster}
#'   object, or numeric vector, representing spatial attractor values across
#'   the simulation region.
#' @param region A \code{Region} or inherited class object defining the spatial
#'   locations included in the spread simulations.
#' @param type One of \code{"source"}, \code{"destination"}, or \code{"both"},
#'   indicating whether attractor values represent scale factors for dispersal
#'   from source locations, or dispersal likelihood weightings for destination
#'   locations, or both.
#' @param multiplier Optional numeric multiplier for scaling or weighting
#'   attractor values appropriately (for convenience).
#' @param ... Additional parameters.
#' @return A \code{Attractor} class object (list) containing functions for
#'   accessing attributes:
#'   \describe{
#'     \item{\code{get_rast()}}{Get the attractor \code{terra::SpatRaster}
#'       object.}
#'     \item{\code{get_values(cells = NULL, na.incl = FALSE)}}{Get the
#'       attractor values (multiplied by the specified multiplier when present)
#'       for optionally specified region cell indices, which should only refer
#'       to all raster indices (including NA-value cells) when specified (via
#'       \code{na.incl = TRUE}), otherwise the index sequence should skip
#'       NA-value cells (as per \code{Region get_indices}). Default is all
#'       non-NA cell values.}
#'     \item{\code{get_aggr_values(cells = NULL, na.incl = FALSE)}}{Get the
#'       raster aggregated attractor values (multiplied by the specified
#'       multiplier when present) for optionally specified aggregated region
#'       cell indices, which should only refer to all raster indices (including
#'       NA-value cells) when specified (via \code{na.incl = TRUE}), otherwise
#'       the index sequence should skip NA-value cells. Default is all
#'       non-NA cell values. Only available for an attractor built with a
#'       spatial raster layer, and functional when region is two-tiered.}
#'     \item{\code{get_type())}}{Get the attractor type ("source",
#'       "destination", or "both").}
#'   }
#' @include Region.R
#' @export
Attractor <- function(x, region,
                      type = c("source",
                               "destination",
                               "both"),
                      multiplier = NULL, ...) {
  UseMethod("Attractor")
}

#' @name Attractor
#' @export
Attractor.Raster <- function(x, ...) {
  # Call the terra version of the function
  Attractor(terra::rast(x), ...)
}

#' @name Attractor
#' @export
Attractor.SpatRaster <- function(x, region,
                                 type = c("source",
                                          "destination",
                                          "both"),
                                 multiplier = NULL, ...) {

  # Check region
  if (!inherits(region, "Region")) {
    stop("Region model must be a 'Region' or inherited class object.",
         call. = FALSE)
  }
  if (!region$is_compatible(x)) {
    stop("The spatial object x should be compatible with that defining the ",
         "region.", call. = FALSE)
  }

  # Type is source, destination, or both
  type <- match.arg(type)

  # Multiplier
  if (!is.numeric(multiplier)) {
    multiplier <- 1
  }

  # Aggregate
  aggr_rast <- NULL

  # Create a class structure
  self <- structure(list(), class = "Attractor")

  # Get the attractor raster
  self$get_rast <- function() {
    return(x*multiplier)
  }

  # Get the values for specified region (non-NA) cell indices
  self$get_values <- function(cells = NULL, na.incl = FALSE) {
    if (is.numeric(cells)) {
      if (na.incl) {
        return(x[cells][,1]*multiplier)
      } else {
        return(x[region$get_indices()[cells]][,1]*multiplier)
      }
    } else {
      if (na.incl) {
        return(x[][,1]*multiplier)
      } else {
        return(x[region$get_indices()][,1]*multiplier)
      }
    }
  }

  # Get aggregate values for specified aggregate region cell indices
  self$get_aggr_values <- function(cells = NULL, na.incl = FALSE) {
    if (region$two_tier()) {
      aggr <- region$get_aggr()
      if (is.null(aggr_rast) && region$two_tier()) {
        aggr_rast <<- terra::aggregate(x, fact = aggr$factor, fun = "mean",
                                       na.rm = TRUE)
      }
      if (is.numeric(cells)) {
        if (na.incl) {
          return(aggr_rast[cells][,1]*multiplier)
        } else {
          return(aggr_rast[aggr$indices[cells]][,1]*multiplier)
        }
      } else {
        if (na.incl) {
          return(aggr_rast[][,1]*multiplier)
        } else {
          return(aggr_rast[aggr$indices][,1]*multiplier)
        }
      }
    }
  }

  # Get the attractor type
  self$get_type <- function() {
    return(type)
  }

  return(self)
}

#' @name Attractor
#' @export
Attractor.numeric <- function(x, region,
                              type = c("source",
                                       "destination",
                                       "both"),
                              multiplier = NULL, ...) {

  # Check region
  if (!inherits(region, "Region")) {
    stop("Region model must be a 'Region' or inherited class object.",
         call. = FALSE)
  }
  if (length(x) != region$get_locations()) {
    stop("The vector x should have values for each location defined in the ",
         "region.", call. = FALSE)
  }

  # Type is source, destination, or both
  type <- match.arg(type)

  # Multiplier
  if (!is.numeric(multiplier)) {
    multiplier <- 1
  }

  # Create a class structure
  self <- structure(list(), class = "Attractor")

  # Get the attractor raster
  self$get_rast <- function() {
    x_rast <- region$get_template(empty = TRUE)
    x_rast[region$get_indices()] <- x*multiplier
    return(x_rast)
  }

  # Get the values for specified region cell indices
  self$get_values <- function(cells = NULL) {
    if (is.numeric(cells)) {
      return(x[cells]*multiplier)
    } else {
      return(x*multiplier)
    }
  }

  # Get the attractor type
  self$get_type <- function() {
    return(type)
  }

  return(self)
}
