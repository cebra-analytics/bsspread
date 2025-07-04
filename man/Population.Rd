% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Population.R
\name{Population}
\alias{Population}
\alias{Population.Region}
\title{Population model class builder}
\usage{
Population(
  region,
  type = c("presence_only", "unstructured", "stage_structured"),
  growth = NULL,
  capacity = NULL,
  capacity_area = NULL,
  establish_pr = NULL,
  incursion_mean = NULL,
  class = character(),
  ...
)

\method{Population}{Region}(
  region,
  type = c("presence_only", "unstructured", "stage_structured"),
  growth = NULL,
  capacity = NULL,
  capacity_area = NULL,
  establish_pr = NULL,
  incursion_mean = NULL,
  class = character(),
  ...
)
}
\arguments{
\item{region}{A \code{Region} or inherited class object representing the
spatial region (template) for spread simulations.}

\item{type}{One of \code{"presence_only"} (default), \code{"unstructured"},
or \code{"stage_structured"} to indicate how populations are represented.}

\item{growth}{Numeric growth rate or lambda (e.g. 1.2 for 20\% growth per
time step) when \code{type = "unstructured"}, or age/stage-based
transition rates (matrix) for each time step when
\code{type = "stage_structured"}. Default is \code{NULL} for when
\code{type = "presence_only"}.}

\item{capacity}{A (static) vector or matrix (containing temporal columns) of
carrying capacity values of the invasive species at each location (row)
specified by the \code{region}, or per unit area defined by
\code{capacity_area} when the \code{region} is spatially implicit. Default
is \code{NULL} for when \code{type = "presence_only"} or growth is not
capacity-limited. The number of columns for temporal capacity should
coincide with the number of simulation time steps, or a cyclic pattern
(e.g. 12 columns for seasonal variation with monthly time steps).}

\item{capacity_area}{Carrying capacity area (in m^2) specified for
capacity-limited growth in a spatially implicit (single patch)
\code{region}. For example, use a value of \code{1e+06} when
\code{capacity} is specified in km^2. Default is \code{NULL}, although a
value is required when \code{capacity} is specified for a spatially
implicit region.}

\item{establish_pr}{An optional (static) vector or matrix (containing
temporal columns) of probability values (0-1) to represent the likelihood
of establishment at each location (row) specified by the \code{region}.
This may be used to avoid transient/unsuccessful incursion or migration
arrivals from being presented in the simulation results, and/or from
subsequently contributing to spread in presence-only models. Default is
\code{NULL}. The number of columns for temporal capacity should
coincide with the number of simulation time steps, or a cyclic pattern
(e.g. 12 columns for seasonal variation with monthly time steps).}

\item{incursion_mean}{Numeric mean population size for unstructured and
stage structured populations applied at incursion locations. The
population size is sampled from the Poisson distribution for each
incursion location.}

\item{class}{Character class name for inherited classes. Default is
\code{NULL}.}

\item{...}{Additional parameters.}
}
\value{
A \code{Population} class object (list) containing functions for
accessing attributes, conforming data to the appropriate population
representation, and simulating growth:
\describe{
\item{\code{get_region()}}{Get the region object.}
\item{\code{get_type()}}{Get the population representation type.}
\item{\code{get_growth()}}{Get the growth rate or age/stage transition
matrix.}
\item{\code{get_stages()}}{Get the number of stages in a age/stage-based
representation.}
\item{\code{get_capacity(cells = NULL, tm = NULL)}}{Get the carrying
capacity as a vector of values for each region location or optionally
specified region locations \code{cells} (indices) at (optional)
simulation time step \code{tm} (for temporally defined capacity).}
\item{\code{get_establish_pr(cells = NULL, tm = NULL)}}{Get the
establishment probability as a vector of values for each region
location or optionally specified region locations \code{cells}
(indices) at (optional) simulation time step \code{tm} (for
temporally defined establishment probability).}
\item{\code{set_incursion_mean(m)}}{Set the incursion mean.}
\item{\code{make(initial, current, incursion, tm)}}{Make a population
vector or array (rows:stages x columns:locations) via the defined
population representation using vectors or arrays of the
\code{initial} or \code{current} and \code{incursion} population at
each region location at simulation time step \code{tm}.}
\item{\code{grow(x, tm)}}{Performs growth or age/stage-based transitions
on population \code{x} vector/array (stages by locations) at
simulation time step \code{tm}, and returns the transformed
vector/array.}
}
}
\description{
Builds a class for population representation and growth functionality for
spread simulations.
}
\references{
Beverton, R. J. H., & Holt, S. J. (1957). On the dynamics of exploited
fish populations. \emph{Fisheries Investigations}, 19, 1-533.

Bradhurst, R., Spring, D., Stanaway, M., Milner, J., & Kompas, T. (2021).
A generalised and scalable framework for modelling incursions,
surveillance and control of plant and environmental pests.
\emph{Environmental Modelling & Software}, 139, N.PAG.
\doi{10.1016/j.envsoft.2021.105004}

García Adeva, J. J., Botha, J. H., & Reynolds, M. (2012). A simulation
modelling approach to forecast establishment and spread of Bactrocera
fruit flies. \emph{Ecological Modelling}, 227, 93–108.
\doi{10.1016/j.ecolmodel.2011.11.026}

Jongejans, E., Skarpaas, O., & Shea, K. (2008). Dispersal, demography and
spatial population models for conservation and control management.
\emph{Perspectives In Plant Ecology Evolution And Systematics}, 9(3–4),
153–170. \doi{10.1016/j.ppees.2007.09.005}

Ricker, W. E. (1958). Handbook of computations for biological statistics
of fish populations. \emph{Bulletin (Fisheries Research Board of Canada)},
119. \url{https://waves-vagues.dfo-mpo.gc.ca/Library/1485.pdf}
}
