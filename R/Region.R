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
#' provided, as well as functionality for calculating weighted paths (graphs)
#' between cells (local and aggregate) via permeability layers. When used,
#' the inverse of the permeability (0-1) of cells is used to (linearly) scale
#' the actual distance between adjacent cells. The effective distance
#' between any two cells in the region is thus calculated via the shortest
#' weighted path between the cells.
#'
#' @param x A \code{raster::RasterLayer} or \code{terra::SpatRaster}
#'   object representing the spatial region (template) for spread simulations.
#' @param ... Additional parameters.
#' @return A \code{Region} class object (list) containing functions for
#'   accessing attributes, checking compatibility of objects with the
#'   region, and to maintain and calculate spatial data for dispersal:
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
#'     \item{\code{is_included(indices)}}{Check if cell \code{indices} of grid
#'       locations are included (non-NA cells) in the simulation (when type is
#'       grid), and return a logical vector indicating the inclusion of each
#'       index.}
#'     \item{\code{two_tier()}}{Check if the region is configured for two-tier
#'       dispersal (local plus aggregate cells) when type is grid.}
#'     \item{\code{get_aggr()}}{Get a list of two-tier dispersal aggregation
#'       components \code{factor}, \code{inner_radius}, \code{rast}
#'       (\code{terra::SpatRaster}), \code{indices} (non-NA), \code{pts}
#'       (\code{terra::SpatVector}), \code{cells} (list of region cells), and a
#'       function \code{get_cells(indices)} that returns the region cells
#'       (indices) within the specified aggregate cell (\code{indices}).}
#'     \item{\code{set_aggr(aggr_factor, inner_radius)}}{Configure a two-tier
#'       dispersal aggregation of the region via an aggregation factor and
#'       an inner radius (in m) to define the boundary between local dispersal
#'       at the original resolution and long-distance dispersal at an aggregate
#'       resolution.}
#'     \item{\code{configure_paths(directions = FALSE, max_distance = NULL,
#'       permeability = NULL)}}{Configure the inclusion of path directions from
#'       occupied to reachable cells, an optional maximum distance (in m) for
#'       the paths, and an optional \code{Permeability} class object containing
#'       a spatial layer for calculating weighted paths, each of which are used
#'       in dispersal calculations. Defaults infer non-inclusion.}
#'     \item{\code{calculate_paths(cells)}}{Calculate the indices, distances
#'       (in m), and directions (0-360 degrees) when required, of cells that
#'       are reachable for each cell index in \code{cells}, including local and
#'       aggregate cells when using two-tier configuration, and a graph
#'       connecting these paths when configured with permeability layers.
#'       The reachable cells are limited via a maximum distance (in m) and/or
#'       permeability layers when either are configured.}
#'     \item{\code{get_paths(cells, directions = FALSE, max_distance = NULL,
#'       perm_id = NULL)}}{Get a list of indices of, and distances (in m) and
#'       (optional) directions (0-360 degrees) to, reachable cells for each
#'       cell index in \code{cells}, including local (\code{$cell}) and
#'       aggregate (\code{$aggr}) cells when configured, and optionally
#'       (further) limited via a maximum distance (in m). Distances, and thus
#'       reachable cells, are optionally modified via weighted paths specified
#'       by a configured permeability object id.}
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

  # Check if cell indices are included (non-NA cells) in the simulation
  self$is_included <- function(indices) {
      return(!is.na(x[indices][,1]))
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
    idx_rast <- terra::rast(x)
    idx_rast[indices] <- 1:length(indices)
    aggr$cells <<- list()
    aggr$rast <<- terra::aggregate(idx_rast, fact = aggr_factor,
                                   fun = function(values) {
      values <- values[which(is.finite(values))]
      if (length(values)) {
        aggr$cells <<- c(aggr$cells, list(values))
        return(length(values))
      } else {
        return(NA)
      }
    })
    names(aggr$rast) <<- "cells"
    aggr$indices <<- which(!is.na(aggr$rast[]))
    aggr$get_cells <<- function(indices) {
      return(unlist(aggr$cells[indices]))
    }
    aggr$pts <<- terra::as.points(aggr$rast, values = FALSE)
    rm(idx_rast)
  }

  # Path list of reachable indices and distances, and optionally configured
  # directions, maximum distance and permeability for dispersal calculations
  paths <- list(idx = list(), distances = list())
  self$configure_paths <- function(directions = FALSE, max_distance = NULL,
                                   permeability = NULL) {
    if (directions && is.null(paths$directions)) {
      paths$directions <<- list()
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
    if (inherits(permeability, "Permeability")) {
      if (is.null(paths$graphs)) {
        paths$graphs <<- list()
      }
      if (is.null(paths$weights)) {
        paths$weights <<- list()
      }
      if (is.null(paths$perms)) {
        paths$perms <<- list(permeability)
        permeability$set_id(1)
      } else {
        paths$perms[[length(paths$perms) + 1]] <<- permeability
        permeability$set_id(length(paths$perms))
      }
    }
  }

  # Calculate reachable region cell indices, distances, and directions (when
  # required), as well as permeability graphs (when required)
  self$calculate_paths <- function(cells) {

    # Reachable cells
    for (cell in cells) {

      # Map path lists use cells as characters
      cell_char <- as.character(cell)

      # Calculate indices when not present
      if (!cell_char %in% names(paths$idx)) {

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
            cell = which(indices %in% terra::cells(x, aggr_poly)[,2] &
                           indices != indices[cell]))

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
                                                     touches = TRUE)[,2] &
                             indices != indices[cell]))
          } else {
            paths$idx[[cell_char]] <<- list(cell = indices[-cell])
          }
        }
      }

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

      # Calculate directions when required and not present
      if (is.list(paths$directions) &&
          !cell_char %in% names(paths$directions)) {

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

    # Graphs for permeability
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

            # Find adjacency of aggregate cells
            cell_adj_df <- rbind(
              cbind(terra::adjacent(aggr$rast,
                                    cells = 1:terra::ncell(aggr$rast),
                                    directions = "rook", pairs = TRUE),
                    weight = 1),
              cbind(terra::adjacent(aggr$rast,
                                    cells = 1:terra::ncell(aggr$rast),
                                    directions = "bishop", pairs = TRUE),
                    weight = sqrt(2)))

            # Create graph for all aggregate cells (including NAs)
            paths$graphs$aggr <<- igraph::graph_from_data_frame(
              cell_adj_df, directed = FALSE)

            # Extract adjacency from graph (different order than above)
            cell_adj_df <- igraph::as_data_frame(paths$graphs$aggr)

            # Calculate path weights for aggregate base/permeability layers
            paths$weights$aggr <<- list(base = cell_adj_df$weight) # faster
            for (perm in paths$perms) {

              # Aggregate the permeability raster layer
              aggr_rast <- terra::aggregate(perm$get_rast(),
                                            fact = aggr$factor, fun = "mean",
                                            na.rm = TRUE)

              # Calculate permeability weights
              perm_weight <- rowMeans(
                1/as.matrix(cbind(aggr_rast[as.integer(cell_adj_df$from)],
                                  aggr_rast[as.integer(cell_adj_df$to)]))
                )*cell_adj_df$weight
              perm_weight[which(is.na(perm_weight))] <- Inf

              # Add to list
              if (is.null(paths$weights$aggr$perms)) {
                paths$weights$aggr$perms <<- list(perm_weight)
              } else {
                index <- length(paths$weights$aggr$perms) + 1
                paths$weights$aggr$perms[[index]] <<- perm_weight
              }
            }
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
          if (length(new_poly)) {
            paths$graphs$poly <<- terra::aggregate(
              terra::union(paths$graphs$poly, new_poly))
          }
        } else {
          paths$graphs$poly <<- new_poly
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

          # Create or update graph
          if (is.null(paths$graphs$cell)) {

            # Create initial graph
            paths$graphs$cell <<-
              igraph::graph_from_data_frame(cell_adj_df, directed = FALSE)

            # Extract adjacency from graph (different order than above)
            cell_adj_df <- igraph::as_data_frame(paths$graphs$cell)

          } else { # union

            # Create graph for new cells
            colnames(cell_adj_df)[3] <- "new_weight"
            new_graph <- igraph::graph_from_data_frame(cell_adj_df,
                                                       directed = FALSE)

            # Union with existing graph (keeps both weight columns)
            new_graph <- igraph::union(paths$graphs$cell, new_graph,
                                       byname = TRUE)

            # Extract adjacency from graph (different order than above)
            cell_adj_df <- igraph::as_data_frame(new_graph)

            # Merge new with existing weights
            new_idx <- which(is.na(cell_adj_df$weight))
            cell_adj_df$weight[new_idx] <- cell_adj_df$new_weight[new_idx]
            igraph::edge_attr(new_graph, "weight") <- cell_adj_df$weight
            paths$graphs$cell <<- igraph::delete_edge_attr(new_graph,
                                                          "new_weight")
          }

          # Create or update path weights for base and permeability layers
          paths$weights$cell <<- list(base = cell_adj_df$weight) # faster
          for (perm in paths$perms) {

            # Calculate permeability weights
            perm_weight <- rowMeans(
              1/as.matrix(cbind(perm$get_rast()[as.integer(cell_adj_df$from)],
                                perm$get_rast()[as.integer(cell_adj_df$to)]))
            )*cell_adj_df$weight
            perm_weight[which(is.na(perm_weight))] <- Inf

            # Add to list
            if (is.null(paths$weights$cell$perms)) {
              paths$weights$cell$perms <<- list(perm_weight)
            } else {
              index <- length(paths$weights$cell$perms) + 1
              paths$weights$cell$perms[[index]] <<- perm_weight
            }
          }
        }

      } else if (is.null(paths$graphs$cell)) { # static graph for all cells

        # Find adjacency of region cells
        cell_adj_df <- rbind(
          cbind(terra::adjacent(x, cells = 1:terra::ncell(x),
                                directions = "rook", pairs = TRUE),
                weight = 1),
          cbind(terra::adjacent(x, cells = 1:terra::ncell(x),
                                directions = "bishop", pairs = TRUE),
                weight = sqrt(2)))

        # Create graph for all region cells (including NAs)
        paths$graphs$cell <<- igraph::graph_from_data_frame(
          cell_adj_df, directed = FALSE)

        # Extract adjacency from graph (different order than above)
        cell_adj_df <- igraph::as_data_frame(paths$graphs$cell)

        # Create path weights for base and permeability layers
        paths$weights$cell <<- list(base = cell_adj_df$weight) # faster
        for (perm in paths$perms) {

          # Calculate permeability weights
          perm_weight <- rowMeans(
            1/as.matrix(cbind(perm$get_rast()[as.integer(cell_adj_df$from)],
                              perm$get_rast()[as.integer(cell_adj_df$to)]))
          )*cell_adj_df$weight
          perm_weight[which(is.na(perm_weight))] <- Inf

          # Add to list
          if (is.null(paths$weights$cell$perms)) {
            paths$weights$cell$perms <<- list(perm_weight)
          } else {
            index <- length(paths$weights$cell$perms) + 1
            paths$weights$cell$perms[[index]] <<- perm_weight
          }
        }
      }
    }
  }

  # Get a list of indices, distances and directions of/to reachable cells for
  # region cells (indices), optionally (further) limited via maximum distance
  # Distances (and reachable cells) are optionally modified via permeability
  self$get_paths <- function(cells, directions = FALSE, max_distance = NULL,
                             perm_id = NULL) {

    # Calculate reachable cell indices (when not present)
    self$calculate_paths(cells)

    # Prepare selected paths as a nested list
    selected <- list(idx = list(), distances = list())

    # Get distances and directions to reachable cells for each cell
    for (cell in cells) {

      # Map path lists use cells as characters
      cell_char <- as.character(cell)

      # Select indices and distances
      selected$idx[[cell_char]] <- paths$idx[[cell_char]]
      selected$distances[[cell_char]] <- paths$distances[[cell_char]]

      # Modify distances via permeability graph weights when required
      limit_paths <- FALSE
      if (is.numeric(perm_id) && is.list(paths$perms) &&
          length(paths$perms) >= perm_id) {

        # Get base (no permeability) weight distance to reachable inner cells
        base_dist <- as.vector(igraph::distances(
          paths$graphs$cell,
          v = as.character(indices[cell]),
          to = as.character(indices[paths$idx[[cell_char]]$cell]),
          weights = paths$weights$cell$base))

        # Get permeability weight distance to reachable inner cells
        perm_dist <- as.vector(igraph::distances(
          paths$graphs$cell,
          v = as.character(indices[cell]),
          to = as.character(indices[paths$idx[[cell_char]]$cell]),
          weights = paths$weights$cell$perms[[perm_id]]))

        # Scale the distance to (otherwise) reachable inner cells
        selected$distances[[cell_char]]$cell <- round(
          selected$distances[[cell_char]]$cell*as.vector(perm_dist/base_dist))

        # Modify aggregate cell distances when applicable
        if (is.list(aggr)) {

          # Find the aggregate cell index that contains the region cell
          aggr_i <- terra::cells(aggr$rast, region_pts[cell],
                                 touches = TRUE)[, "cell"]

          # Get base weight distance to reachable aggregate cells
          base_dist <- as.vector(igraph::distances(
            paths$graphs$aggr,
            v = as.character(aggr_i),
            to = as.character(aggr$indices[paths$idx[[cell_char]]$aggr]),
            weights = paths$weights$aggr$base))

          # Get permeability weight distance to reachable aggregate cells
          perm_dist <- as.vector(igraph::distances(
            paths$graphs$aggr,
            v = as.character(aggr_i),
            to = as.character(aggr$indices[paths$idx[[cell_char]]$aggr]),
            weights = paths$weights$aggr$perms[[perm_id]]))

          # Scale the distance to (otherwise) reachable inner cells
          selected$distances[[cell_char]]$aggr <- round(
            selected$distances[[cell_char]]$aggr*as.vector(perm_dist/base_dist))
        }

        limit_paths <- TRUE
      }

      # Limit to paths within a maximum distance or reachable cells if required
      if (limit_paths ||
          is.numeric(max_distance) && max_distance < paths$max_distance) {

        # Resolve maximum distance when present
        if (is.numeric(max_distance) && is.numeric(paths$max_distance)) {
          max_distance <- min(max_distance, paths$max_distance)
        } else if (is.numeric(paths$max_distance)) {
          max_distance <- paths$max_distance
        }

        # Find cells within range
        if (is.numeric(max_distance)) {
          cell_ids <-
            which(selected$distances[[cell_char]]$cell <= max_distance)
        } else {
          cell_ids <- which(is.finite(selected$distances[[cell_char]]$cell))
        }

        # Update selected indices and distances to in-range cells
        selected$idx[[cell_char]]$cell <-
          selected$idx[[cell_char]]$cell[cell_ids]
        selected$distances[[cell_char]]$cell <-
          as.integer(selected$distances[[cell_char]]$cell[cell_ids])

        # Update selected aggregate indices and distances
        if (is.list(aggr)) {

          # Find cells within range
          if (is.numeric(max_distance)) {
            aggr_ids <-
              which(selected$distances[[cell_char]]$aggr <= max_distance)
          } else {
            aggr_ids <- which(is.finite(selected$distances[[cell_char]]$aggr))
          }

          # Update selected indices and distances to in-range aggregate cells
          selected$idx[[cell_char]]$aggr <-
            selected$idx[[cell_char]]$aggr[aggr_ids]
          selected$distances[[cell_char]]$aggr <-
            as.integer(selected$distances[[cell_char]]$aggr[aggr_ids])
        }

        limit_paths <- TRUE
      }

      # Get directions when required
      if (is.list(paths$directions) && directions) {

        # Select directions
        selected$directions[[cell_char]] <- paths$directions[[cell_char]]

        # Limit within range paths when required
        if (limit_paths) {
          selected$directions[[cell_char]]$cell <-
            selected$directions[[cell_char]]$cell[cell_ids]
          if (is.list(aggr)) {
            selected$directions[[cell_char]]$aggr <-
              selected$directions[[cell_char]]$aggr[aggr_ids]
          }
        }
      }
    }

    return(selected)
  }

  return(self)
}
