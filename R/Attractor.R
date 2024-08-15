#' Attractor class builder
#'
#' Builds a class to represent spatially distributed weights or scaling factors
#' assigned to cells/locations to represent an "attractor" that modify the
#' relative likelihood of dispersal to each destination location.
#'
#' @param x A \code{raster::RasterLayer} or \code{terra::SpatRaster}
#'   object, or numeric vector, representing spatial attractor values across
#'   the simulation region.
#' @param region A \code{Region} or inherited class object defining the spatial
#'   locations included in the spread simulations.
#' @param ... Additional parameters.
#' @return A \code{Attractor} class object (list) containing functions for
#'   accessing attributes:
#'   \describe{
#'     \item{\code{get_rast()}}{Get the attractor \code{terra::SpatRaster}
#'       object.}
#'     \item{\code{get_values(cells = NULL, na.incl = FALSE)}}{Get the
#'       attractor values for optionally specified region cell indices, which
#'       should only refer to all raster indices (including NA-value cells)
#'       when specified (via \code{na.incl = TRUE}), otherwise the index
#'       sequence should skip NA-value cells (as per
#'       \code{Region get_indices}). Default is all non-NA cell values.}
#'     \item{\code{get_aggr_values(cells = NULL, na.incl = FALSE)}}{Get the
#'       raster aggregated attractor values for optionally specified aggregated
#'       region cell indices, which should only refer to all raster indices
#'       (including NA-value cells) when specified (via \code{na.incl = TRUE}),
#'       otherwise the index sequence should skip NA-value cells. Default is
#'       all non-NA cell values. Only available for an attractor built with a
#'       spatial raster layer, and functional when region is two-tiered.}
#'   }
#' @references
#'   Bossenbroek, J. M., Kraft, C. E., & Nekola, J. C. (2001). Prediction of
#'   Long-Distance Dispersal Using Gravity Models: Zebra Mussel Invasion of
#'   Inland Lakes. \emph{Ecological Applications}, 11(6), 1778–1788.
#'   \doi{10.2307/3061095}
#'
#'   Carrasco, L. R., Mumford, J. D., MacLeod, A., Harwood, T.,
#'   Grabenweger, G., Leach, A. W., Knight, J. D., & Baker, R. H. A. (2010).
#'   Unveiling human-assisted dispersal mechanisms in invasive alien insects:
#'   Integration of spatial stochastic simulation and phenology models.
#'   \emph{Ecological Modelling}, 221(17), 2068–2075.
#'   \doi{10.1016/j.ecolmodel.2010.05.012}
#'
#'   Crespo-Pérez, V., Rebaudo, F., Silvain, J.-F., & Dangles, O. (2011).
#'   Modeling invasive species spread in complex landscapes: the case of potato
#'   moth in Ecuador. \emph{Landscape Ecology}, 26(10), 1447.
#'   \doi{10.1007/s10980-011-9649-4}
#'
#'   Muirhead, J. R., Leung, B., Overdijk, C., Kelly, D. W., Nandakumar, K.,
#'   Marchant, K. R., & MacIsaac, H. J. (2006). Modelling local and
#'   long-distance dispersal of invasive emerald ash borer Agrilus planipennis
#'   (Coleoptera) in North America. \emph{Diversity & Distributions}, 12(1),
#'   71–79. \doi{10.1111/j.1366-9516.2006.00218.x}
#' @include Region.R
#' @export
Attractor <- function(x, region, ...) {
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
Attractor.SpatRaster <- function(x, region, ...) {

  # Check region
  if (!inherits(region, "Region")) {
    stop("Region model must be a 'Region' or inherited class object.",
         call. = FALSE)
  }
  if (!region$is_compatible(x)) {
    stop("The spatial object x should be compatible with that defining the ",
         "region.", call. = FALSE)
  }

  # Aggregate
  aggr_rast <- NULL

  # Create a class structure
  self <- structure(list(), class = "Attractor")

  # Get the attractor raster
  self$get_rast <- function() {
    return(x)
  }

  # Get the values for specified region (non-NA) cell indices
  self$get_values <- function(cells = NULL, na.incl = FALSE) {
    if (is.numeric(cells)) {
      if (na.incl) {
        return(x[cells][,1])
      } else {
        return(x[region$get_indices()[cells]][,1])
      }
    } else {
      if (na.incl) {
        return(x[][,1])
      } else {
        return(x[region$get_indices()][,1])
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
          return(aggr_rast[cells][,1])
        } else {
          return(aggr_rast[aggr$indices[cells]][,1])
        }
      } else {
        if (na.incl) {
          return(aggr_rast[][,1])
        } else {
          return(aggr_rast[aggr$indices][,1])
        }
      }
    }
  }

  return(self)
}

#' @name Attractor
#' @export
Attractor.numeric <- function(x, region, ...) {

  # Check region
  if (!inherits(region, "Region")) {
    stop("Region model must be a 'Region' or inherited class object.",
         call. = FALSE)
  }
  if (length(x) != region$get_locations()) {
    stop("The vector x should have values for each location defined in the ",
         "region.", call. = FALSE)
  }

  # Create a class structure
  self <- structure(list(), class = "Attractor")

  # Get the attractor raster
  self$get_rast <- function() {
    if (region$get_type() == "grid") {
      x_rast <- region$get_template(empty = TRUE)
      x_rast[region$get_indices()] <- x
      return(x_rast)
    } else {
      return(NULL)
    }
  }

  # Get the values for specified region cell indices
  self$get_values <- function(cells = NULL) {
    if (is.numeric(cells)) {
      return(x[cells])
    } else {
      return(x)
    }
  }

  return(self)
}
