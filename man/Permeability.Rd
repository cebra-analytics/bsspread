% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Permeability.R
\name{Permeability}
\alias{Permeability}
\alias{Permeability.Raster}
\alias{Permeability.SpatRaster}
\alias{Permeability.matrix}
\alias{Permeability.data.frame}
\title{Permeability class builder}
\usage{
Permeability(x, region, ...)

\method{Permeability}{Raster}(x, ...)

\method{Permeability}{SpatRaster}(x, region, ...)

\method{Permeability}{matrix}(x, ...)

\method{Permeability}{data.frame}(x, region, ...)
}
\arguments{
\item{x}{A \code{raster::RasterLayer} or \code{terra::SpatRaster}
object representing the grid-based spatial permeability, or a 3-column
data frame (or matrix) of patch index pairs and permeability values
(weights) for connected (adjacent) paths between patches, consistent with
the spread simulation region.}

\item{region}{A \code{Region} or inherited class object defining the spatial
locations included in the spread simulations.}

\item{...}{Additional parameters.}
}
\value{
A \code{Permeability} class object (list) containing functions for
accessing attributes:
\describe{
\item{\code{get_id()}}{Get the object numeric identifier.}
\item{\code{set_id(id)}}{Set the object numeric identifier.}
\item{\code{get_rast()}}{Get the permeability \code{terra::SpatRaster}
object when the simulation region is grid-based.}
\item{\code{get_data()}}{Get the path data specifying the weights
between patches when the simulation region is patch-based.}
}
}
\description{
Builds a class to represent the permeability of cells within a grid-based
spatial region, or paths between patches in a patch-based region, of a
spread simulation. The inverse of the permeability (0-1) of cells/paths is
used to (linearly) scale the actual distance between adjacent cells/patches.
For example, a permeability value of \code{0.5} results in an effective
distance twice that of the actual distance, a value of \code{0} prevents
spread to or through a cell, or spread along a network path to a patch, and
a value of \code{1} does not modify the effective distance.
}
