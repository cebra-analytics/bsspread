% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulator.R
\name{Simulator}
\alias{Simulator}
\alias{Simulator.Raster}
\alias{Simulator.SpatRaster}
\alias{Simulator.Region}
\title{Simulator class builder}
\usage{
Simulator(
  region,
  time_steps = 1,
  step_duration = 1,
  step_units = "years",
  collation_steps = 1,
  replicates = 1,
  initializer = NULL,
  population_model = NULL,
  dispersal_models = list(),
  ...
)

\method{Simulator}{Raster}(region, ...)

\method{Simulator}{SpatRaster}(region, ...)

\method{Simulator}{Region}(
  region,
  time_steps = 1,
  step_duration = 1,
  step_units = "years",
  collation_steps = 1,
  replicates = 1,
  initializer = NULL,
  population_model = NULL,
  dispersal_models = list(),
  ...
)
}
\arguments{
\item{region}{A \code{raster::RasterLayer}, \code{terra::SpatRaster}, or
\code{Region} or inherited class object representing the spatial region
(template) for spread simulations.}

\item{time_steps}{The number of discrete time steps to simulate. Default is
1.}

\item{step_duration}{The duration of the simulation time steps in units
specified by \code{step_units}. Default is 1.}

\item{step_units}{The units for the simulation step duration
(\code{step_duration}) as a character string. Default is "years".}

\item{collation_steps}{The interval in time steps for collating results.
Default is 1, that is, results are collated at every time step.}

\item{replicates}{The number of replicate or repeated simulations to be run.
Default is 1. Note that replicate simulations results are collated as
summary statistics across simulations.}

\item{initializer}{A \code{Initializer} or inherited class object for
generating the initial invasive species population distribution or
incursion locations, as well as optionally generating subsequent
incursions during the spread simulations.}

\item{population_model}{A \code{Population} or inherited class object
defining the population representation and growth functionality for the
spread simulations.}

\item{dispersal_models}{A list of \code{Dispersal} or inherited class
objects defining the dispersal functionality for the different spread
vectors to be simulated.}

\item{...}{Additional parameters.}
}
\value{
A \code{Simulator} class object (list) containing functions for
setting objects (in the function environment) and running the simulations:
\describe{
\item{\code{set_initializer(object)}}{Set the initializer object.}
\item{\code{set_population_model(model)}}{Set the population model
object.}
\item{\code{set_dispersal_models(models)}}{Set the list of dispersal
model objects.}
\item{\code{run()}}{Run the simulations and return the results.}
}
}
\description{
Builds a class to configure and run replicate discrete-time spread
simulations over a given spatial region using (sub)models for population
growth and dispersal.
}