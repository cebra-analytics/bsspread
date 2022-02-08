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
#'   accessing attributes, checking compatibility of objects with the
#'   region, and to maintain common spatial data for dispersal:
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
#'     \item{\code{configure_paths(directions = FALSE, graphs = FALSE,
#'       max_distance = NULL)}}{Configure the inclusion of path directions from
#'       occupied to reachable cells, the inclusion of \code{graphs} for
#'       representing weighted (permeability) paths, and an optional maximum
#'       distance (in m) for the paths, which are used in dispersal
#'       calculations. Defaults infer non-inclusion.}
#'     \item{\code{calculate_paths(cells)}}{Calculate the cells that are
#'       reachable for each cell index in \code{cells}, including local and
#'       aggregate cells when configured, and a graph connecting these paths
#'       when configured. These are both limited via a maximum distance (in m)
#'       when configured.}
#'     \item{\code{get_paths(cells, directions = FALSE, max_distance = NULL)}}{
#'       Get a list of indices, distances (in m), and directions (0-360
#'       degrees) when configured, of/to reachable cells for each cell index in
#'       \code{cells}, including local (\code{$cell}) and aggregate
#'       (\code{$aggr}) cells when configured, and optionally (further) limited
#'       via a maximum distance (in m).}
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

  # Path list of reachable indices and distances, and optionally configured
  # directions, graphs and maximum distance for dispersal calculations
  paths <- list(idx = list(), distances = list())
  self$configure_paths <- function(directions = FALSE, graphs = FALSE,
                                   max_distance = NULL) {
    if (directions) {
      paths$directions <<- list()
    }
    if (graphs) {
      paths$graphs <<- list()
    }
    if (is.numeric(max_distance)) { # set if empty, or replace if larger
      if (is.null(paths$max_distance) ||
          (is.numeric(paths$max_distance) &&
           max_distance > paths$max_distance)) {
        paths$max_distance <<- max_distance
      }
    } else { # use full extent of region
      paths$max_distance <<- Inf
    }
  }

  # Calculate reachable cells and graphs for region cells (indices)
  self$calculate_paths <- function(cells) {

    # Reachable cells
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
          if (is.numeric(paths$max_distance) &&
              is.finite(paths$max_distance)) {
            outer_vect <- terra::buffer(region_pts[cell,],
                                        width = paths$max_distance,
                                        quadsegs = 180)
            outer_idx <- terra::cells(aggr$rast, outer_vect,
                                      touches = TRUE)[,2]
            outer_idx <- outer_idx[which(!outer_idx %in% aggr_idx)]
            paths$idx[[cell_char]]$aggr <<- which(aggr$indices %in% outer_idx)
          } else {
            paths$idx[[cell_char]]$aggr <<- which(!aggr$indices %in% aggr_idx)
          }

        } else { # cell approach

          # Select cells within range when applicable
          if (is.numeric(paths$max_distance) &&
              is.finite(paths$max_distance)) {
            range_vect <- terra::buffer(region_pts[cell,],
                                        width = paths$max_distance,
                                        quadsegs = 180)
            paths$idx[[cell_char]] <<- list(
              cell = which(indices %in% terra::cells(x, range_vect,
                                                     touches = TRUE)[,2]))
          } else {
            paths$idx[[cell_char]] <<- list(cell = indices)
          }
        }
      }
    }

    # Graphs
    if (is.list(paths$graphs)) {

      # Maintain a graph to all region cells within reach
      if (is.list(aggr) ||
          (is.numeric(paths$max_distance) && is.finite(paths$max_distance))) {

        if (is.list(aggr)) { # two-tier approach

          # Select aggregate cells in/on intersected inner circles
          inner_vect <- terra::buffer(region_pts[cells,],
                                      width = aggr$inner_radius,
                                      quadsegs = 180)
          aggr_idx <- terra::cells(aggr$rast, inner_vect, touches = TRUE)[,2]

          # Create a polygon from aggregate cells in/on the inner circles
          new_poly <- terra::rast(aggr$rast)
          new_poly[aggr_idx] <- 0
          new_poly <- terra::as.polygons(new_poly)

          # Build a graph for the full aggregation (once)
          if (is.null(paths$graphs$aggr)) {
            # TODO ####
          }

        } else { # cell approach

          # Create a polygon from intersected range circles
          new_poly <- terra::buffer(region_pts[cells,],
                                    width = paths$max_distance,
                                    quadsegs = 180)
        }

        # Determine new graph coverage and update total via polygons
        if (!is.null(paths$graphs$poly)) {
          new_poly <- terra::erase(new_poly, paths$graphs$poly)
          paths$graphs$poly <- terra::union(paths$graphs$poly, new_poly)
        } else {
          paths$graphs$poly <- new_poly
        }

        # Calculate region cell indices inside new polygon
        cell_idx <- terra::cells(x, new_poly, touches = TRUE)[,2]
        if (length(cell_idx)) {

          # Find adjacency of region cells
          cell_adj_df <- rbind(
            cbind(terra::adjacent(x, cells = cell_idx, directions = "rook",
                                  pairs = TRUE), weight = 1),
            cbind(terra::adjacent(x, cells = cell_idx, directions = "bishop",
                                  pairs = TRUE), weight = sqrt(2)))

          # Create and set or update graph
          new_graph <- igraph::graph_from_data_frame(cell_adj_df,
                                                     directed = FALSE)
          if (is.null(paths$graphs$cell)) {
            paths$graphs$cell <- new_graph
          } else { # union
            paths$graphs$cell <- igraph::union(paths$graphs$cell,
                                               new_graph, byname = TRUE)
          }
        }

      } else if (is.null(paths$graphs)) { # static graph for all cells
        # TODO ####
      }
    }
  }

  # Get a list of indices, distances and directions of/to reachable cells for
  # region cells (indices), optionally (further) limited via maximum distance
  self$get_paths <- function(cells, directions = FALSE, max_distance = NULL) {

    # Calculate reachable cell indices (when not present)
    self$calculate_paths(cells)

    # Prepare selected paths as a nested list
    selected <- list(idx = list(), distances = list())

    # Get distances and directions to reachable cells for each cell
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

      # Limit to paths within a maximum distance when required
      if (is.numeric(max_distance) && max_distance < paths$max_distance) {
        cell_ids <- which(paths$distances[[cell_char]]$cell <= max_distance)
        selected$idx[[cell_char]] <- list(
          cell = paths$idx[[cell_char]]$cell[cell_ids])
        selected$distances[[cell_char]] <- list(
          cell = paths$distances[[cell_char]]$cell[cell_ids])
        if (is.list(aggr)) {
          aggr_ids <- which(paths$distances[[cell_char]]$aggr <= max_distance)
          selected$idx[[cell_char]]$aggr <-
            paths$idx[[cell_char]]$aggr[aggr_ids]
          selected$distances[[cell_char]]$aggr <-
            paths$distances[[cell_char]]$aggr[aggr_ids]
        }
      } else { # use full range
        selected$idx[[cell_char]] <- paths$idx[[cell_char]]
        selected$distances[[cell_char]] <- paths$distances[[cell_char]]
      }

      # Get directions when required
      if (is.list(paths$directions) && directions) {

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

        # Limit to paths within a maximum distance when required (again)
        if (is.numeric(max_distance) && max_distance < paths$max_distance) {
          selected$directions[[cell_char]] <- list(
            cell = paths$directions[[cell_char]]$cell[cell_ids])
          if (is.list(aggr)) {
            selected$directions[[cell_char]]$aggr <-
              paths$directions[[cell_char]]$aggr[aggr_ids]
          }
        } else { # use full range
          selected$directions[[cell_char]] <- paths$directions[[cell_char]]
        }
      }
    }

    return(selected)
  }

  # Get graph (shortest path) weights
  self$get_graph_weights <- function() {
    # TODO
  }

  # IN PROGRESS

  return(self)
}
