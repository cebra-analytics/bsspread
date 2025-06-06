#' Unstructured population model class builder
#'
#' Builds a class for representing unstructured populations in spread
#' simulations. Simulates logistic (capacity-limited) growth.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the spread simulations.
#' @param growth Numeric growth rate or lambda (e.g. 1.2 for 20% growth per
#'   time step). Default is \code{1.0} for no increase.
#' @param capacity A (static) vector or matrix (containing temporal columns) of
#'   carrying capacity values of the invasive species at each location (row)
#'   specified by the \code{region}, or per unit area defined by
#'   \code{capacity_area} when the \code{region} is spatially implicit. Default
#'   is \code{NULL} for when growth is not capacity-limited. The number of
#'   columns for temporal capacity should coincide with the number of
#'   simulation time steps, or a cyclic pattern (e.g. 12 columns for seasonal
#'   variation with monthly time steps).
#' @param capacity_area Carrying capacity area (in m^2) specified for
#'   capacity-limited growth in a spatially implicit (single patch)
#'   \code{region}. For example, use a value of \code{1e+06} when
#'   \code{capacity} is specified in km^2. Default is \code{NULL}, although a
#'   value is required when \code{capacity} is specified for a spatially
#'   implicit region.
#' @param establish_pr An optional (static) vector or matrix (containing
#'   temporal columns) of probability values (0-1) to represent the likelihood
#'   of establishment at each location (row) specified by the \code{region}.
#'   This may be used to avoid transient/unsuccessful incursion or migration
#'   arrivals from being presented in the simulation results. Default is
#'   \code{NULL}. The number of columns for temporal capacity should
#'   coincide with the number of simulation time steps, or a cyclic pattern
#'   (e.g. 12 columns for seasonal variation with monthly time steps).
#' @param incursion_mean Numeric mean population size for incursion locations.
#'   The population size is sampled from the Poisson distribution for each
#'   incursion location.
#' @param ... Additional parameters.
#' @return A \code{UnstructPopulation} class object (list) containing functions
#'   for accessing attributes and simulating growth:
#'   \describe{
#'     \item{\code{get_region()}}{Get the region object.}
#'     \item{\code{get_type()}}{Get the population representation type.}
#'     \item{\code{get_growth()}}{Get the growth rate.}
#'     \item{\code{get_capacity(cells = NULL, tm = NULL)}}{Get the carrying
#'       capacity as a vector of values for each region location or optionally
#'       specified region locations \code{cells} (indices) at (optional)
#'       simulation time step \code{tm} (for temporally defined capacity).}
#'     \item{\code{get_establish_pr(cells = NULL, tm = NULL)}}{Get the
#'       establishment probability as a vector of values for each region
#'       location or optionally specified region locations \code{cells}
#'       (indices) at (optional) simulation time step \code{tm} (for
#'       temporally defined establishment probability).}
#'     \item{\code{set_incursion_mean(m)}}{Set the incursion mean.}
#'     \item{\code{make(initial, current, incursion, tm)}}{Make a population
#'       vector via using vectors of the \code{initial} or \code{current} and
#'       \code{incursion} population at each region location at simulation time
#'       step \code{tm}.}
#'     \item{\code{grow(x, tm)}}{Performs logistic (capacity-limited) growth on
#'       the population \code{x} vector at simulation time step \code{tm}, and
#'       returns the transformed vector.}
#'   }
#' @references
#'   Beverton, R. J. H., & Holt, S. J. (1957). On the dynamics of exploited
#'   fish populations. \emph{Fisheries Investigations}, 19, 1-533.
#'
#'   Bradhurst, R., Spring, D., Stanaway, M., Milner, J., & Kompas, T. (2021).
#'   A generalised and scalable framework for modelling incursions,
#'   surveillance and control of plant and environmental pests.
#'   \emph{Environmental Modelling & Software}, 139, N.PAG.
#'   \doi{10.1016/j.envsoft.2021.105004}
#'
#'   Jongejans, E., Skarpaas, O., & Shea, K. (2008). Dispersal, demography and
#'   spatial population models for conservation and control management.
#'   \emph{Perspectives In Plant Ecology Evolution And Systematics}, 9(3–4),
#'   153–170. \doi{10.1016/j.ppees.2007.09.005}
#'
#'   Ricker, W. E. (1958). Handbook of computations for biological statistics
#'   of fish populations. \emph{Bulletin (Fisheries Research Board of Canada)},
#'   119. \url{https://waves-vagues.dfo-mpo.gc.ca/Library/1485.pdf}
#' @include Population.R
#' @export
UnstructPopulation <- function(region,
                               growth = 1.0,
                               capacity = NULL,
                               capacity_area = NULL,
                               establish_pr = NULL,
                               incursion_mean = NULL, ...) {

  # Build via base class
  self <- Population(region,
                     type = "unstructured",
                     growth = growth,
                     capacity = capacity,
                     capacity_area = capacity_area,
                     establish_pr = establish_pr,
                     incursion_mean = incursion_mean,
                     class = "UnstructPopulation")

  # Grow method - override for logistic (capacity-limited) growth
  self$grow <- function(x, tm) {

    # Indices of occupied locations
    indices <- which(x > 0)

    # Get capacity at time step
    capacity_tm <- self$get_capacity(tm = tm)

    # Calculate logistic growth rates
    if (is.numeric(capacity_tm)) { # capacity-limited

      # Remove populations at locations having zero capacity
      if (any(capacity_tm[indices] <= 0)) {
        zero_idx <- indices[which(capacity_tm[indices] <= 0)]
        x[zero_idx] <- 0
        indices <- indices[!indices %in% zero_idx]
      }

      # Calculate capacity for spatially implicit diffusion or area spread
      if (region$spatially_implicit()) {

        # Diffusion
        if (is.numeric(attr(x, "diffusion_rate")) &&
            is.numeric(attr(x, "diffusion_radius"))) {

          # Calculate capacity of diffusion area
          capacity_radius <-
            attr(x, "diffusion_radius") + attr(x, "diffusion_rate")
          area_capacity <- capacity_tm*pi*capacity_radius^2/capacity_area

          # Calculate capacity-limited growth rate
          r <- exp(log(growth)*(1 - x/area_capacity))

        } else if (is.numeric(capacity_area) &&
                   is.numeric(region$get_max_implicit_area())) {

          # Calculate capacity of maximum area
          area_capacity <- (capacity_tm*region$get_max_implicit_area()/
                              capacity_area)

          # Calculate capacity-limited growth rate
          r <- exp(log(growth)*(1 - x/area_capacity))

        } else {
          r <- growth # unlimited
        }

      } else {

        # Calculate capacity-limited growth rates for each occupied location
        r <- exp(log(growth)*(1 - x[indices]/capacity_tm[indices]))
      }

    } else {
      r <- growth
    }

    # Sample the new population values via the Poisson distribution
    if (length(indices)) {
      x[indices] <- stats::rpois(length(indices), r*x[indices])
    }

    return(x)
  }

  return(self)
}
