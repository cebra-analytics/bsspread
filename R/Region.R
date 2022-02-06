#' Region class builder
#'
#' Builds a class to represent a spatial region for a spread simulation,
#' defined by a raster layer with active (non-NA) cells, or by a network of
#' locations or patches. When defined by a raster layer, the region can be
#' configured for two-tier dispersal, whereby long-distance dispersal is
#' calculated using a courser (aggregate) spatial resolution for destinations
#' outside of an inner radius, within which local dispersal is calculated using
#' the original spatial resolution of the region. Functions for calculating
#' distances and directions between (local and aggregate) cells are also
#' provided, as well as functionality for calculating network paths (graphs)
#' between cells (local and aggregate).
#'
#' @param x A \code{raster::RasterLayer} or \code{terra::SpatRaster}
#'   object representing the spatial region (template) for spread simulations.
#' @param ... Additional parameters.
#' @return A \code{Region} class object (list) containing functions for
#'   accessing attributes and checking compatibility of objects with the
#'   region:
#'   \describe{
#'     \item{\code{get_template(empty = FALSE)}}{Get the spatial template with
#'      either zeros in non-NA locations (default), or with no values when
#'      \code{empty = TRUE}.}
#'     \item{\code{get_indices()}}{Get cell indices of grid locations that are
#'       included in the simulation.}
#'     \item{\code{get_locations()}}{Get the number of locations (cells or
#'       patches) that are included in the simulation.}
#'     \item{\code{get_type()}}{Get the type of spatial representation: grid
#'       (raster cells) or patch (network).}
#'     \item{\code{is_compatible(y)}}{Check the compatibility of object
#'       \code{y} with the region defined by \code{x}.}
#'     \item{\code{two_tier()}}{Check if the region is configured for two-tier
#'       dispersal (local plus aggregate cells) when type is grid.}
#'     \item{\code{get_aggr()}}{Get a list of two-tier dispersal aggregation
#'       components \code{factor}, \code{inner_radius}, \code{rast}
#'       (\code{terra::SpatRaster}), \code{indices} (non-NA), and \code{pts}
#'       (\code{terra::SpatVector}).}
#'     \item{\code{set_aggr(aggr_factor, inner_radius)}}{Configure a two-tier
#'       dispersal aggregation of the region via an aggregation factor and
#'       an inner radius (in m) to define the boundary between local dispersal
#'       at the original resolution and long-distance dispersal at an aggregate
#'       resolution.}
#'     \item{\code{get_reachable(cells, max_distance = NULL)}}{Get the indices
#'       of cells that are reachable for each cell index in \code{cells},
#'       including local (\code{$cell}) and aggregate (\code{$aggr}) cells when
#'       configured, and optionally limited via a maximum distance (in m).}
#'     \item{\code{get_distances(cells, max_distance = NULL)}}{Get the distance
#'       (in m) of cells that are reachable from each cell index in
#'       \code{cells}, including local (\code{$cell}) and aggregate
#'       (\code{$aggr}) cells when configured, and optionally limited via a
#'       maximum distance (in m).}
#'     \item{\code{get_directions(cells, max_distance = NULL)}}{Get the
#'       direction (0-360 degrees) of cells that are reachable from each cell
#'       index in \code{cells}, including local (\code{$cell}) and aggregate
#'       (\code{$aggr}) cells when configured, and optionally limited via a
#'       maximum distance (in m).}
#'   }
#' @export
Region <- function(x, ...) {
  UseMethod("Region")
}

#' @name Region
#' @export
Region.Raster <- function(x, ...) {
  # Call the terra version of the function
  Region(terra::rast(x), ...)
}

#' @name Region
#' @export
Region.SpatRaster <- function(x, ...) {

  # Non-NA cell indices and points (terra::SpatVector)
  indices <- which(!is.na(x[]))
  region_pts <- terra::as.points(x, values = FALSE)

  # Create a class structure
  self <- structure(list(), class = "Region")

  # Get spatial template with zero/NA or empty values
  self$get_template <- function(empty = FALSE) {
    if (empty) {
      template <- terra::rast(x)
    } else {
      template <- x*0
      names(template) <- "value"
    }
    return(template)
  }

  # Get cell indices
  self$get_indices <- function() {
    return(indices)
  }

  # Get the number of active cell locations
  self$get_locations <- function() {
    return(length(indices))
  }

  # Get the spatial region type
  self$get_type <- function() {
    return("grid")
  }

  # Check compatibility of a spatial raster y with the region defined by x
  self$is_compatible <- function(y) {
    y <- terra::rast(y)
    return(terra::crs(y) == terra::crs(x) &&
             terra::ext(y) == terra::ext(x) &&
             all(terra::res(y) == terra::res(x)))
  }

  # Check, get, and set two-tier aggregation (list) for dispersal
  aggr <- NULL
  self$two_tier <- function() {
    return(is.list(aggr))
  }
  self$get_aggr <- function() {
    return(aggr)
  }
  self$set_aggr <- function(aggr_factor, inner_radius) {
    aggr <<- list(factor = aggr_factor, inner_radius = inner_radius)
    aggr$rast <<- terra::aggregate(x, fact = aggr_factor, fun = "mean",
                                   na.rm = TRUE)
    aggr$indices <<- which(!is.na(aggr$rast[]))
    aggr$pts <<- terra::as.points(aggr$rast, values = FALSE)
  }

  # Path list of distances, directions, and graphs for dispersal calculations
  paths <- list(idx = list(), distances = list(), directions = list(),
                graphs = list())

  # Calculate reachable cells for region cells (indices) [internal function]
  calculate_reachable <- function(cells, max_distance) {

    for (cell in cells) {

      # Calculate when not present
      if (!as.character(cell) %in% names(paths$idx)) {

        # Map path lists use cells as characters
        cell_char <- as.character(cell)

        if (is.list(aggr)) { # two-tier approach

          # Select aggregate cells in/on inner circle
          inner_vect <- terra::buffer(region_pts[cell,],
                                      width = aggr$inner_radius,
                                      quadsegs = 180)
          aggr_idx <- terra::cells(aggr$rast, inner_vect, touches = TRUE)[,2]

          # Create a polygon from aggregate cells in/on inner circle
          aggr_poly <- terra::rast(aggr$rast)
          aggr_poly[aggr_idx] <- 0
          aggr_poly <- terra::as.polygons(aggr_poly)

          # Local cell indices inside aggregate polygon
          paths$idx[[cell_char]] <<- list(
            cell = which(indices %in% terra::cells(x, aggr_poly)[,2]))

          # Aggregate cell indices outside polygon
          if (is.numeric(max_distance)) {
            outer_vect <- terra::buffer(region_pts[cell,],
                                        width = max_distance, quadsegs = 180)
            outer_idx <- terra::cells(aggr$rast, outer_vect,
                                      touches = TRUE)[,2]
            outer_idx <- outer_idx[which(!outer_idx %in% aggr_idx)]
            paths$idx[[cell_char]]$aggr <<- which(aggr$indices %in% outer_idx)
          } else {
            paths$idx[[cell_char]]$aggr <<- which(!aggr$indices %in% aggr_idx)
          }

        } else { # cell approach

          # Select cells within range when applicable
          if (is.numeric(max_distance)) {
            range_vect <- terra::buffer(region_pts[cell,],
                                        width = max_distance, quadsegs = 180)
            paths$idx[[cell_char]] <<- list(
              cell = which(indices %in% terra::cells(x, range_vect,
                                                     touches = TRUE)[,2]))
          } else {
            paths$idx[[cell_char]] <<- list(cell = indices)
          }
        }
      }
    }
  }

  # Get reachable cells for region cells (indices)
  self$get_reachable <- function(cells, max_distance = NULL) {

    # Calculate reachable cell indices (when not present)
    calculate_reachable(cells, max_distance)

    # Return from paths list
    return(paths$idx[as.character(cells)])
  }

  # Get distances for region cells (indices)
  self$get_distances <- function(cells, max_distance = NULL) {

    # Calculate reachable cell indices (when not present)
    calculate_reachable(cells, max_distance)

    # Get distances to reachable cells for each cell
    for (cell in cells) {

      # Map path lists use cells as characters
      cell_char <- as.character(cell)

      # Calculate distances when not present
      if (!cell_char %in% names(paths$distances)) {

        # Calculate (local) cell distances
        paths$distances[[cell_char]] <<- list(
          cell = as.integer(round(as.numeric(
            terra::distance(region_pts[cell],
                            region_pts[paths$idx[[cell_char]]$cell])))))

        # Calculate aggregate cell distances when applicable
        if (is.list(aggr)) {
          paths$distances[[cell_char]]$aggr <<-
            as.integer(round(as.numeric(
              terra::distance(region_pts[cell],
                              aggr$pts[paths$idx[[cell_char]]$aggr]))))
        }
      }
    }

    # Return from paths list
    return(paths$distances[as.character(cells)])
  }

  # Get directions for region cells (indices)
  self$get_directions <- function(cells, max_distance = NULL) {

    # Calculate reachable cell indices (when not present)
    calculate_reachable(cells, max_distance)

    # Get directions to reachable cells for each cell
    for (cell in cells) {

      # Map path lists use cells as characters
      cell_char <- as.character(cell)

      # Calculate directions when not present
      if (!cell_char %in% names(paths$directions)) {

        # Calculate (local) cell directions
        xy_diff <- (terra::crds(region_pts[cell])[
          rep(1, length(paths$idx[[cell_char]]$cell)),] -
            terra::crds(region_pts[paths$idx[[cell_char]]$cell]))
        paths$directions[[cell_char]] <<- list(cell = as.integer(round(
          atan2(xy_diff[,"y"], xy_diff[,"x"])*180/pi + 180)))

        # Calculate aggregate cell directions when applicable
        if (is.list(aggr)) {
          xy_diff <- (terra::crds(region_pts[cell])[
            rep(1, length(paths$idx[[cell_char]]$aggr)),] -
              terra::crds(aggr$pts[paths$idx[[cell_char]]$aggr]))
          paths$directions[[cell_char]]$aggr <<- as.integer(round(
            atan2(xy_diff[,"y"], xy_diff[,"x"])*180/pi + 180))
        }
      }
    }

    # Return from paths list
    return(paths$directions[as.character(cells)])
  }

  # Get graph (shortest path) weights
  self$get_graph_weights <- function() {
    # TODO
  }

  # IN PROGRESS

  return(self)
}
