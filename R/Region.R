#' Region class builder
#'
#' Builds a class to represent a spatial region for a spread simulation,
#' defined by a raster layer with active (non-NA) cells, or by a network of
#' locations or patches. When defined by a grid-based raster layer, the region
#' can be configured for two-tier dispersal, whereby long-distance dispersal is
#' calculated using a courser (aggregate) spatial resolution for destinations
#' outside of an inner radius, within which local dispersal is calculated using
#' the original spatial resolution of the region. Functions for calculating
#' distances and directions between (local and aggregate) cells or patches are
#' also provided, as well as functionality for calculating weighted paths
#' (graphs) between cells (local and aggregate) or patches via permeability
#' layers/networks. When used, the inverse of the permeability (0-1) of
#' cells/patches is used to (linearly) scale the actual distance between
#' adjacent/connected cells/patched. The effective distance between any two
#' cells/patches in the region is thus calculated via the shortest weighted
#' path between the cells/patches.
#'
#' @param x A \code{raster::RasterLayer} or \code{terra::SpatRaster} object
#'   representing a grid-based spatial region (template) for spread simulations.
#'   Alternatively, a network of locations or patches may be defined via a data
#'   frame (or matrix) of location coordinates in longitude and latitude
#'   (WGS84) with explicitly named columns "lon" and "lat".
#' @param ... Additional parameters.
#' @return A \code{Region} class object (list) containing functions for
#'   accessing attributes, checking compatibility of objects with the
#'   region, and to maintain and calculate spatial data for dispersal:
#'   \describe{
#'     \item{\code{get_type()}}{Get the type of spatial representation: "grid"
#'       (raster cells) or "patch" (network).}
#'     \item{\code{get_locations()}}{Get the number of locations (cells or
#'       patches) that are included in the simulation.}
#'     \item{\code{is_compatible(y)}}{Check the compatibility of object
#'       \code{y} with the region defined by \code{x}.}
#'     \item{\code{get_template(empty = FALSE)}}{Get the spatial template when
#'      the \code{type} is "grid", with either zeros in non-NA locations
#'      (default), or with no values when \code{empty = TRUE}.}
#'     \item{\code{get_indices()}}{Get cell indices of grid or patch locations
#'       that are included in the simulation when the \code{type} is "grid".}
#'     \item{\code{get_res()}}{Get the spatial cell resolution (in m) of the
#'       region when the \code{type} is "grid".}
#'     \item{\code{is_included(indices)}}{Check if cell \code{indices} of grid
#'       locations are included (non-NA cells) in the simulation when the
#'       \code{type} is "grid", and return a logical vector indicating the
#'       inclusion of each index.}
#'     \item{\code{two_tier()}}{Check if the region is configured for two-tier
#'       dispersal (local plus aggregate cells) when the \code{type} is
#'       "grid".}
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
#'       occupied to reachable locations (cells or patches), an optional
#'       maximum distance (in m) for the paths, and an optional
#'       \code{Permeability} class object containing a spatial layer (for
#'       grids), or weighted connectivity data (for patches), for calculating
#'       weighted paths, each of which are used in dispersal calculations.
#'       Defaults infer non-inclusion.}
#'     \item{\code{calculate_paths(cells|patches)}}{Calculate the indices,
#'       distances (in m), and directions (0-360 degrees) when required, of
#'       locations (cells or patches) that are reachable for each location
#'       index in \code{cells} or \code{patches}, including local and aggregate
#'       cells when using a two-tier grid configuration. Paths are optionally
#'       limited via a maximum distance (in m). When configured with
#'       permeability objects, a graph of path connections is derived and used
#'       to calculate modified distances and accessibility to reachable
#'       locations for each permeability object.}
#'     \item{\code{get_paths(cells|patches, directions = FALSE,
#'       max_distance = NULL, perm_id = NULL)}}{Get a list of indices of, and
#'       distances (in m) and (optional) directions (0-360 degrees) to,
#'       reachable locations (cells or patches) for each location index in
#'       \code{cells} or \code{patches}, including local (\code{$cell}) and
#'       aggregate (\code{$aggr}) cells when using a two-tier grid
#'       configuration, and optionally (further) limited via a maximum distance
#'       (in m) and permeability. When a permeability object (id) is specified,
#'       modified distances (and accessibility) to reachable locations is
#'       included in the returned list (with name (\code{perm_dist}).}
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

  # Get the spatial region type
  self$get_type <- function() {
    return("grid")
  }

  # Get the number of active cell locations
  self$get_locations <- function() {
    return(length(indices))
  }

  # Check compatibility of a spatial raster y with the region defined by x
  self$is_compatible <- function(y) {
    y <- terra::rast(y)
    return(terra::crs(y) == terra::crs(x) &&
             terra::ext(y) == terra::ext(x) &&
             all(terra::res(y) == terra::res(x)))
  }

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

  # Get the spatial cell resolution
  self$get_res <- function() {
    if (terra::is.lonlat(x)) { # EPSG:4326
      corners <- array(terra::ext(x2), c(2, 2))
      diagonal <- terra::distance(corners[1,,drop = FALSE],
                                  corners[2,,drop = FALSE], lonlat = TRUE)
      return(diagonal/sqrt(terra::nrow(x2)^2 + terra::ncol(x2)^2))
    } else {
      return(mean(terra::res(x)[1]))
    }
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
    aggr$rast <<- terra::aggregate(
      idx_rast, fact = aggr_factor, fun = function(values) {
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
        if (is.list(aggr)) {
          idx_rast <- terra::rast(x)
          idx_rast[] <- 1:terra::ncell(x)
          paths$graphs$agg_idx <<- list()
          idx_rast <-
            terra::aggregate(idx_rast, fact = aggr$factor,
                             fun = function(values) {
                               paths$graphs$agg_idx <<-
                                 c(paths$graphs$agg_idx, list(values))
                               return(NA)
                             })
          paths$graphs$get_cells <<- function(indices) {
            return(unlist(paths$graphs$agg_idx[indices]))
          }
          rm(idx_rast)
        }
      }
      if (is.null(paths$weights)) {
        paths$weights <<- list()
      }
      if (is.null(paths$perms)) {
        paths$perms <<- list(permeability)
        permeability$set_id(1)
      } else if (is.null(permeability$get_id())) {
        paths$perms[[length(paths$perms) + 1]] <<- permeability
        permeability$set_id(length(paths$perms))
      }
      if (is.null(paths$perm_dist)) {
        paths$perm_dist <<- list()
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

          # Local cell indices within inner area
          paths$idx[[cell_char]] <<- list(
            cell = aggr$get_cells(which(aggr$indices %in% aggr_idx)))
          cell_i <- which(paths$idx[[cell_char]]$cell == cell)
          paths$idx[[cell_char]]$cell <<- paths$idx[[cell_char]]$cell[-cell_i]

          # Aggregate cells outside inner area
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
          inner_vect <- terra::aggregate(
            terra::buffer(region_pts[cells,], width = aggr$inner_radius,
                          quadsegs = 180))
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
            paths$weights$aggr <<- list(base = cell_adj_df$weight,
                                        perms = list()) # faster
            for (perm_id in 1:length(paths$perms)) {

              # Aggregate the permeability raster layer
              aggr_rast <- terra::aggregate(paths$perms[[perm_id]]$get_rast(),
                                            fact = aggr$factor, fun = "mean",
                                            na.rm = TRUE)

              # Calculate permeability weights
              perm_weight <- rowMeans(
                1/as.matrix(cbind(aggr_rast[as.integer(cell_adj_df$from)],
                                  aggr_rast[as.integer(cell_adj_df$to)]))
              )*cell_adj_df$weight
              perm_weight[which(is.na(perm_weight))] <- Inf

              # Add to list
              paths$weights$aggr$perms[[perm_id]] <<- perm_weight
            }
          }

        } else { # cell approach

          # Create a polygon from intersected range circles
          new_poly <- terra::aggregate(
            terra::buffer(region_pts[cells,], width = paths$max_distance,
                          quadsegs = 180))
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
        if (is.list(aggr)) { # two-tier approach
          aggr_idx <- terra::cells(aggr$rast, new_poly, touches = FALSE)[,2]
          cell_idx <- paths$graphs$get_cells(aggr_idx)
        } else {
          cell_idx <- terra::cells(x, new_poly, touches = TRUE)[,2]
        }

        # Build graph via cell adjacency
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

          } else { # extend via union

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
          paths$weights$cell <<- list(base = cell_adj_df$weight,
                                      perms = list()) # faster
          for (perm_id in 1:length(paths$perms)) {

            # Calculate permeability weights
            perm <- paths$perms[[perm_id]]
            perm_weight <- rowMeans(
              1/as.matrix(cbind(perm$get_rast()[as.integer(cell_adj_df$from)],
                                perm$get_rast()[as.integer(cell_adj_df$to)]))
            )*cell_adj_df$weight
            perm_weight[which(is.na(perm_weight))] <- Inf

            # Add to list
            paths$weights$cell$perms[[perm_id]] <<- perm_weight
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
        paths$weights$cell <<- list(base = cell_adj_df$weight,
                                    perms = list()) # faster
        for (perm_id in 1:length(paths$perms)) {

          # Calculate permeability weights
          perm <- paths$perms[[perm_id]]
          perm_weight <- rowMeans(
            1/as.matrix(cbind(perm$get_rast()[as.integer(cell_adj_df$from)],
                              perm$get_rast()[as.integer(cell_adj_df$to)]))
          )*cell_adj_df$weight
          perm_weight[which(is.na(perm_weight))] <- Inf

          # Add to list
          paths$weights$cell$perms[[perm_id]] <<- perm_weight
        }
      }

      # Calculate permeability-modified distances for reachable cells
      for (cell in cells) {

        # Map path lists use cells as characters
        cell_char <- as.character(cell)

        # Calculate modified distances when not present
        if (!cell_char %in% names(paths$perm_dist)) {

          # List for modified distances
          paths$perm_dist[[cell_char]] <<- list()

          # Get base (no permeability) weight distance to reachable inner cells
          base_dist <- as.vector(igraph::distances(
            paths$graphs$cell,
            v = as.character(indices[cell]),
            to = as.character(indices[paths$idx[[cell_char]]$cell]),
            weights = paths$weights$cell$base))

          # Calculate modified distances for each permeability layer
          paths$perm_dist[[cell_char]]$cell <<- list()
          for (perm_id in 1:length(paths$perms)) {

            # Get permeability weight distance to reachable inner cells
            perm_dist <- as.vector(igraph::distances(
              paths$graphs$cell,
              v = as.character(indices[cell]),
              to = as.character(indices[paths$idx[[cell_char]]$cell]),
              weights = paths$weights$cell$perms[[perm_id]]))

            # Calculate the distance modifiers
            perm_dist <- perm_dist/base_dist
            perm_dist[which(!is.finite(perm_dist))] <- NA

            # Scale the distance to (otherwise) reachable inner cells
            paths$perm_dist[[cell_char]]$cell[[perm_id]] <<- as.integer(
              round(paths$distances[[cell_char]]$cell*perm_dist))
          }

          # Aggregate distance multipliers when applicable
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

            # Calculate modified distances for each permeability layer
            paths$perm_dist[[cell_char]]$aggr <<- list()
            for (perm_id in 1:length(paths$perms)) {

              # Get permeability weight distance to reachable aggregate cells
              perm_dist <- as.vector(igraph::distances(
                paths$graphs$aggr,
                v = as.character(aggr_i),
                to = as.character(aggr$indices[paths$idx[[cell_char]]$aggr]),
                weights = paths$weights$aggr$perms[[perm_id]]))

              # Calculate the distance modifiers
              perm_dist <- perm_dist/base_dist
              perm_dist[which(!is.finite(perm_dist))] <- NA

              # Scale the distance to (otherwise) reachable inner cells
              paths$perm_dist[[cell_char]]$aggr[[perm_id]] <<- as.integer(
                round(paths$distances[[cell_char]]$aggr*perm_dist))
            }
          }
        }
      }
    }
  }

  # Get a list of indices, distances and directions of/to reachable cells for
  # region cells (indices), optionally (further) limited via maximum distance
  # and/or permeability modified distances.
  self$get_paths <- function(cells, directions = FALSE, max_distance = NULL,
                             perm_id = NULL) {

    # Using permeability?
    use_perm <- (is.numeric(perm_id) && is.list(paths$perms) &&
                   length(paths$perms) >= perm_id)

    # Prepare selected paths as a nested list
    selected <- list(idx = list(), distances = list())
    if (use_perm) {
      selected$perm_dist <- list()
    }

    # Get distances and directions to reachable cells for each cell
    for (cell in cells) {

      # Map path lists use cells as characters
      cell_char <- as.character(cell)

      # Select indices and distances
      selected$idx[[cell_char]] <- paths$idx[[cell_char]]
      selected$distances[[cell_char]] <- paths$distances[[cell_char]]

      # Get modified distances for permeability (id)
      limit_paths <- FALSE
      if (use_perm) {
        selected$perm_dist[[cell_char]]$cell <-
          paths$perm_dist[[cell_char]]$cell[[perm_id]]
        if (is.list(aggr)) {
          selected$perm_dist[[cell_char]]$aggr <-
            paths$perm_dist[[cell_char]]$aggr[[perm_id]]
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

        # Find cells within range (including adjacent cells)
        if (use_perm) {
          if (is.numeric(max_distance)) {
            cell_ids <-
              which(selected$perm_dist[[cell_char]]$cell <= max_distance |
                    (selected$distances[[cell_char]]$cell <=
                       1.6*region$get_res()))
            cell_ids <- cell_ids[
              which(is.finite(selected$perm_dist[[cell_char]]$cell[cell_ids]))]
          } else {
            cell_ids <- which(is.finite(selected$perm_dist[[cell_char]]$cell))
          }
        } else {
          if (is.numeric(max_distance)) {
            cell_ids <-
              which(selected$distances[[cell_char]]$cell <= max_distance)
          } else {
            cell_ids <- which(is.finite(selected$distances[[cell_char]]$cell))
          }
        }

        # Update selected indices and distances to in-range cells
        selected$idx[[cell_char]]$cell <-
          selected$idx[[cell_char]]$cell[cell_ids]
        selected$distances[[cell_char]]$cell <-
          selected$distances[[cell_char]]$cell[cell_ids]
        if (use_perm) {
          selected$perm_dist[[cell_char]]$cell <-
            selected$perm_dist[[cell_char]]$cell[cell_ids]
        }

        # Update selected aggregate indices and distances
        if (is.list(aggr)) {

          # Find aggregate cells within range
          if (use_perm) {
            if (is.numeric(max_distance)) {
              aggr_ids <-
                which(selected$perm_dist[[cell_char]]$aggr <= max_distance)
            } else {
              aggr_ids <- which(is.finite(selected$perm_dist[[cell_char]]$aggr))
            }
          } else {
            if (is.numeric(max_distance)) {
              aggr_ids <-
                which(selected$distances[[cell_char]]$aggr <= max_distance)
            } else {
              aggr_ids <- which(is.finite(selected$distances[[cell_char]]$aggr))
            }
          }

          # Update selected indices and distances to in-range aggregate cells
          selected$idx[[cell_char]]$aggr <-
            selected$idx[[cell_char]]$aggr[aggr_ids]
          selected$distances[[cell_char]]$aggr <-
            selected$distances[[cell_char]]$aggr[aggr_ids]
          if (use_perm) {
            selected$perm_dist[[cell_char]]$aggr <-
              selected$perm_dist[[cell_char]]$aggr[aggr_ids]
          }
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

#' @name Region
#' @export
Region.matrix <- function(x, ...) {

  # Call the data frame version of the function
  Region(as.data.frame(x), ...)
}

#' @name Region
#' @export
Region.data.frame <- function(x, ...) {

  # Check data frame
  if (!all(c("lon", "lat") %in% names(x))) {
    stop("Coordinate data frame must contain columns named 'lon' and 'lat'.",
         call. = FALSE)
  }

  # Region points (terra::SpatVector)
  region_pts <- terra::vect(x[, c("lon", "lat")], crs = "EPSG:4326")

  # Create a class structure
  self <- structure(list(), class = "Region")

  # Get the spatial region type
  self$get_type <- function() {
    return("patch")
  }

  # Get the number of patch locations
  self$get_locations <- function() {
    return(nrow(x))
  }

  # Check compatibility of vector, matrix, or adjacency data frame y
  # with the region defined by x
  self$is_compatible <- function(y) {
    if (is.data.frame(y)) {
      return(ncol(y) == 3 && all(unique(unlist(y[,1:2])) %in% 1:nrow(x)))
    } else {
      y <- as.matrix(y)
      return(nrow(y) == nrow(x) && ncol(y) %in% c(1, nrow(x)))
    }
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
      if (is.null(paths$perms)) {
        paths$perms <<- list(permeability)
        permeability$set_id(1)
      } else if (is.null(permeability$get_id())) {
        paths$perms[[length(paths$perms) + 1]] <<- permeability
        permeability$set_id(length(paths$perms))
      }
    }
  }

  # Calculate reachable region patch indices, distances, and directions (when
  # required), as well as permeability graphs (when required)
  self$calculate_paths <- function(patches) {

    # Reachable patches
    for (patch in patches) {

      # Map path lists use patches as characters
      patch_char <- as.character(patch)

      # Calculate indices when not present
      if (!patch_char %in% names(paths$idx)) {

          # Select patches within range when applicable
          if (is.numeric(paths$max_distance) &&
              is.finite(paths$max_distance)) {
            range_vect <- terra::buffer(region_pts[patch,],
                                        width = paths$max_distance,
                                        quadsegs = 180)
            paths$idx[[patch_char]] <<- which(
              terra::relate(region_pts, range_vect, "within")[,1] &
                1:nrow(x) != patch)
          } else {
            paths$idx[[patch_char]] <<- (1:nrow(x))[-patch]
          }
      }

      # Calculate distances when not present
      if (!patch_char %in% names(paths$distances)) {
        paths$distances[[patch_char]] <<- as.integer(round(as.numeric(
            terra::distance(region_pts[patch],
                            region_pts[paths$idx[[patch_char]]]))))
      }

      # Calculate directions when required and not present
      if (is.list(paths$directions) &&
          !patch_char %in% names(paths$directions)) {
        xy_diff <- (terra::crds(region_pts[patch])[
          rep(1, length(paths$idx[[patch_char]])),] -
            terra::crds(region_pts[paths$idx[[patch_char]]]))
        paths$directions[[patch_char]] <<- as.integer(round(
          atan2(xy_diff[,"y"], xy_diff[,"x"])*180/pi + 180))
      }
    }

    # Graphs for permeability
    if (is.list(paths$perms)) {

      # Build a graph for all region patches (once)
      if (is.null(paths$graphs)) {

        # Create a merged graph
        for (perm_id in 1:length(paths$perms)) {

          # Get permeability
          perm <- paths$perms[[perm_id]]

          # Create or update graph with indexed weight attribute
          new_graph <- igraph::graph_from_data_frame(perm$get_data()[,1:2],
                                                     directed = FALSE)
          igraph::edge_attr(new_graph,
                            paste0("w", perm_id)) <- perm$get_data()[,3]
          if (is.null(paths$graphs)) {
            paths$graphs <<- new_graph
          } else {
            paths$graphs <<- igraph::union(paths$graphs, new_graph)
          }
        }

        # Get merged permeability data (in graph order)
        graph_df <- igraph::as_data_frame(paths$graphs)

        # Calculate distances between connected patches
        graph_df$distance <- as.integer(round(
          terra::distance(region_pts[as.integer(graph_df$from)],
                          region_pts[as.integer(graph_df$to)],
                          pairwise = TRUE)))

        # Calculated weighted equivalent distances (used as graph weights)
        paths$weights <<- list()
        for (perm_id in 1:length(paths$perms)) {
          paths$weights[[perm_id]] <<- as.integer(round(
            graph_df$distance/graph_df[,paste0("w", perm_id)]))
        }
      }

      # IN PROGRESS ####

      # Calculate permeability-modified distances for reachable cells
      for (cell in cells) {

        # Map path lists use cells as characters
        cell_char <- as.character(cell)

        # Calculate modified distances when not present
        if (!cell_char %in% names(paths$perm_dist)) {

          # List for modified distances
          paths$perm_dist[[cell_char]] <<- list()

          # Get base (no permeability) weight distance to reachable inner cells
          base_dist <- as.vector(igraph::distances(
            paths$graphs$cell,
            v = as.character(indices[cell]),
            to = as.character(indices[paths$idx[[cell_char]]$cell]),
            weights = paths$weights$cell$base))

          # Calculate modified distances for each permeability layer
          paths$perm_dist[[cell_char]]$cell <<- list()
          for (perm_id in 1:length(paths$perms)) {

            # Get permeability weight distance to reachable inner cells
            perm_dist <- as.vector(igraph::distances(
              paths$graphs$cell,
              v = as.character(indices[cell]),
              to = as.character(indices[paths$idx[[cell_char]]$cell]),
              weights = paths$weights$cell$perms[[perm_id]]))

            # Calculate the distance modifiers
            perm_dist <- perm_dist/base_dist
            perm_dist[which(!is.finite(perm_dist))] <- NA

            # Scale the distance to (otherwise) reachable inner cells
            paths$perm_dist[[cell_char]]$cell[[perm_id]] <<- as.integer(
              round(paths$distances[[cell_char]]$cell*perm_dist))
          }

          # Aggregate distance multipliers when applicable
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

            # Calculate modified distances for each permeability layer
            paths$perm_dist[[cell_char]]$aggr <<- list()
            for (perm_id in 1:length(paths$perms)) {

              # Get permeability weight distance to reachable aggregate cells
              perm_dist <- as.vector(igraph::distances(
                paths$graphs$aggr,
                v = as.character(aggr_i),
                to = as.character(aggr$indices[paths$idx[[cell_char]]$aggr]),
                weights = paths$weights$aggr$perms[[perm_id]]))

              # Calculate the distance modifiers
              perm_dist <- perm_dist/base_dist
              perm_dist[which(!is.finite(perm_dist))] <- NA

              # Scale the distance to (otherwise) reachable inner cells
              paths$perm_dist[[cell_char]]$aggr[[perm_id]] <<- as.integer(
                round(paths$distances[[cell_char]]$aggr*perm_dist))
            }
          }
        }
      }
    }
  }

  return(self)
}
