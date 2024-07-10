#' Area spread dispersal model class builder
#'
#' Builds a class for representing spatially implicit area spread, whereby the
#' simulated occupied area grows in proportion to the size of an unstructured
#' or stage-based population given its capacity per unit area. When a maximum
#' total area occupied is specified (via the \code{Region} object), then both
#' the population and its area of occupancy will follow a typical
#' capacity-limited sigmoid growth curve.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatially implicit (single location) region (template) for spread
#'   simulations. The region object contains functionality for setting the
#'   maximum total area occupied.
#' @param population_model A \code{UnstructPopulation}, \code{StagePopulation},
#'   or inherited class object defining the population representation for the
#'   spread simulations.
#' @param ... Additional parameters.
#' @return A \code{AreaSpread} class object (list), containing inherited and
#'   extended functions from the generic \code{Dispersal} class for accessing
#'   attributes (of the function environment) and performing area spread
#'   dispersal:
#'   \describe{
#'     \item{\code{pack(n)}}{Packs a population vector or matrix \code{n} into
#'       a list containing a vector of occupied \code{cells} (indices), the
#'       \code{original} population values at the occupied locations only, the
#'       \code{remaining} occupied values (initially a duplicate
#'       \code{original}), and a vector or matrix for the \code{relocated}
#'       population values at all locations (initially all zero).}
#'     \item{\code{unpack(n)}}{Unpacks a population list by combining the
#'       \code{remaining} and \code{relocated} population values to form a
#'       new post-dispersal population vector or matrix.}
#'     \item{\code{disperse(n)}}{Perform area spread dispersal on a list
#'       \code{n} of vectors or matrices, representing the single occupied
#'       location (\code{cells}) (indices), the \code{original} population, the
#'       \code{remaining} population, and the \code{relocated} population, and
#'       return the transformed list of vectors or matrices. Spatially implicit
#'       spread attaches a \code{diffusion_radius} attribute to \code{n} at
#'       each time step indicative of the area of occupancy.}
#'   }
#' @references
#'   Andow, D. A., Kareiva, P. M., Levin, S. A., & Okubo, A. (1990). Spread of
#'   invading organisms. \emph{Landscape Ecology}, 4(2–3), 177.
#'   \doi{10.1007/bf00132860}
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
#'   Robinet, C., Kehlenbeck, H., Kriticos, D. J., Baker, R. H. A.,
#'   Battisti, A., Brunel, S., Dupin, M., Eyre, D., Faccoli, M., Ilieva, Z.,
#'   Kenis, M., Knight, J., Reynaud, P., Yart, A., & van der Werf, W. (2012).
#'   A Suite of Models to Support the Quantitative Assessment of Spread in Pest
#'   Risk Analysis. \emph{PLoS ONE}, 7(10), 1–18.
#'   \doi{10.1371/journal.pone.0043366}
#' @include Region.R
#' @include Dispersal.R
#' @export
AreaSpread <- function(region, population_model, ...) {

  # Check that region is spatially implicit
  if (!region$spatially_implicit()) {
    stop("The region must be spatially implicit for area spread.",
         call. = FALSE)
  }

  # Check that the population is unstructured or staged
  if (!population_model$get_type() %in% c("unstructured",
                                          "stage_structured")) {
    stop(paste("The population model must be unstructured or stage structured",
               "for area spread."), call. = FALSE)
  }

  # Check that the population capacity is specified
  if (!is.numeric(population_model$get_capacity())) {
    stop("The population capacity must be specified for area spread.",
         call. = FALSE)
  }

  # Build via base class
  self <- Dispersal(region, population_model,
                    dispersal_stages = NULL,
                    proportion = NULL,
                    events = NULL,
                    distance_function = NULL,
                    direction_function = NULL,
                    combined_function = NULL,
                    distance_adjust = NULL,
                    attractors = list(),
                    permeability = NULL,
                    max_distance = NULL,
                    class = "AreaSpread", ...)

  # Override disperse function
  self$disperse <- function(n) {

    # Calculate spread area from population size
    capacity <- population_model$get_capacity()
    capacity_area <- attr(capacity, "area")
    if (population_model$get_type() == "stage_structured") {
      stages <- population_model$get_capacity_stages()
      spread_area <- sum(n$remaining[,stages])/capacity*capacity_area
    } else { # unstructured
      spread_area <- n$remaining/capacity*capacity_area
    }

    # Limit via maximum area when applicable
    if (is.numeric(region$get_max_implicit_area())) {
      spread_area <- min(spread_area, region$get_max_implicit_area())
    }

    # Attach attribute for spread area
    attr(n$relocated, "spread_area") <- spread_area

    return(n)
  }

  return(self)
}
