% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UnstructPopulation.R
\name{UnstructPopulation}
\alias{UnstructPopulation}
\title{Unstructured population model class builder}
\usage{
UnstructPopulation(
  region,
  growth = 1,
  capacity = NULL,
  capacity_area = NULL,
  establish_pr = NULL,
  incursion_mean = NULL,
  ...
)
}
\arguments{
\item{region}{A \code{Region} or inherited class object representing the
spatial region (template) for the spread simulations.}

\item{growth}{Numeric growth rate or lambda (e.g. 1.2 for 20\% growth per
time step). Default is \code{1.0} for no increase.}

\item{capacity}{A vector of carrying capacity values of the invasive species
at each location specified by the \code{region}, or per unit area defined
by \code{capacity_area} when the \code{region} is spatially implicit.
Default is \code{NULL} for when growth is not capacity-limited.}

\item{capacity_area}{Carrying capacity area (in m^2) specified for
capacity-limited growth in a spatially implicit (single patch)
\code{region}. For example, use a value of \code{1e+06} when
\code{capacity} is specified in km^2. Default is \code{NULL}, although a
value is required when \code{capacity} is specified for a spatially
implicit region.}

\item{establish_pr}{An optional vector of probability values (0-1) to
represent the likelihood of establishment at each location specified by
the \code{region}. This may be used to avoid transient/unsuccessful
incursions or migrations from being presented in the simulation results.
Default is \code{NULL}.}

\item{incursion_mean}{Numeric mean population size for incursion locations.
The population size is sampled from the Poisson distribution for each
incursion location.}

\item{...}{Additional parameters.}
}
\value{
A \code{UnstructPopulation} class object (list) containing functions
for accessing attributes and simulating growth:
\describe{
\item{\code{get_region()}}{Get the region object.}
\item{\code{get_type()}}{Get the population representation type.}
\item{\code{get_growth()}}{Get the growth rate.}
\item{\code{get_capacity()}}{Get the carrying capacity as a vector of
values for each location.}
\item{\code{get_establish_pr()}}{Get the establishment probability as a
vector of values for each location.}
\item{\code{set_incursion_mean(m)}}{Set the incursion mean.}
\item{\code{make(initial, current, incursion)}}{Make a population vector
via using vectors of the \code{initial} or \code{current} and
\code{incursion} population at each region location.}
\item{\code{grow(x)}}{Performs logistic (capacity-limited) growth on the
population \code{x} vector, and returns the transformed vector.}
}
}
\description{
Builds a class for representing unstructured populations in spread
simulations. Simulates logistic (capacity-limited) growth.
}
\references{
Beverton, R. J. H., & Holt, S. J. (1957). On the dynamics of exploited
fish populations. \emph{Fisheries Investigations}, 19, 1-533.

Bradhurst, R., Spring, D., Stanaway, M., Milner, J., & Kompas, T. (2021).
A generalised and scalable framework for modelling incursions,
surveillance and control of plant and environmental pests.
\emph{Environmental Modelling & Software}, 139, N.PAG.
\doi{10.1016/j.envsoft.2021.105004}

Jongejans, E., Skarpaas, O., & Shea, K. (2008). Dispersal, demography and
spatial population models for conservation and control management.
\emph{Perspectives In Plant Ecology Evolution And Systematics}, 9(3–4),
153–170. \doi{10.1016/j.ppees.2007.09.005}

Ricker, W. E. (1958). Handbook of computations for biological statistics
of fish populations. \emph{Bulletin (Fisheries Research Board of Canada)},
119. \url{https://waves-vagues.dfo-mpo.gc.ca/Library/1485.pdf}
}
