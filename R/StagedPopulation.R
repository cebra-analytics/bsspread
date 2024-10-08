#' Stage or age-based population model class builder
#'
#' Builds a class for representing stage or age-based populations in spread
#' simulations. Simulates logistic (capacity-limited) growth.
#'
#' @param region A \code{Region} or inherited class object representing the
#'   spatial region (template) for the spread simulations.
#' @param growth A stage or age-based matrix containing reproduction and
#'   survival rates for each time step. An optional attribute \code{survivals}
#'   may be used to differentiate the survival rates, otherwise survival rates
#'   are assumed to be any non-zero values on the lower (left) triangle of the
#'   matrix, including those on the diagonal. The stages/ages may also be
#'   labelled via an optional attribute \code{labels}.
#' @param capacity A vector of carrying capacity values of the invasive species
#'   at each location specified by the \code{region}, or per unit area defined
#'   by \code{capacity_area} when the \code{region} is spatially implicit.
#'   Default is \code{NULL} for when growth is not capacity-limited.
#' @param capacity_area Carrying capacity area (in m^2) specified for
#'   capacity-limited growth in a spatially implicit (single patch)
#'   \code{region}. For example, use a value of \code{1e+06} when
#'   \code{capacity} is specified in km^2. Default is \code{NULL}, although a
#'   value is required when \code{capacity} is specified for a spatially
#'   implicit region.
#' @param capacity_stages A vector of stage or age indices (as per the
#'   rows/columns of \code{growth} matrix) indicating which life-stages/ages
#'   are applicable for capacity-limited growth (and survival). Default is
#'   \code{NULL} for when growth is not capacity-limited, or when all
#'   stages/ages are applicable for capacity-limited growth.
#' @param establish_pr An optional vector of probability values (0-1) to
#'   represent the likelihood of establishment at each location specified by
#'   the \code{region}. This may be used to avoid transient/unsuccessful
#'   incursions or migrations from being presented in the simulation results.
#'   Default is \code{NULL}.
#' @param incursion_mean Numeric mean population size for incursion locations.
#'   The population size is sampled from the Poisson distribution for each
#'   incursion location.
#' @param incursion_stages A vector of stage or age indices (as per the
#'   rows/columns of \code{growth} matrix) indicating which life-stages/ages
#'   are applicable for incursions. Default is \code{NULL} for when all
#'   stages/ages are applicable for incursions.
#' @param ... Additional parameters.
#' @return A \code{StagedPopulation} class object (list) containing functions
#'   for accessing attributes and simulating growth:
#'   \describe{
#'     \item{\code{get_region()}}{Get the region object.}
#'     \item{\code{get_type()}}{Get the population representation type.}
#'     \item{\code{get_growth()}}{Get the growth stage/age matrix.}
#'     \item{\code{get_growth_r()}}{Get the equivalent (single value) growth
#'       rate for the \code{growth} stage/age matrix.}
#'     \item{\code{get_capacity()}}{Get the carrying capacity as a vector of
#'       values for each location.}
#'     \item{\code{get_capacity_stages()}}{Get the stage/age (indices) that are
#'       applicable for capacity-limited growth.}
#'     \item{\code{get_establish_pr()}}{Get the establishment probability as a
#'       vector of values for each location.}
#'     \item{\code{set_incursion_mean(m)}}{Set the incursion mean.}
#'     \item{\code{set_incursion_stages(s)}}{Set the incursion stages.}
#'     \item{\code{make(initial, current, incursion)}}{Make a population matrix
#'       (with a row per location and a column per stage/age) via using vectors
#'       or matrices of the \code{initial} or \code{current} and
#'       \code{incursion} population at each region location.}
#'     \item{\code{grow(x)}}{Performs logistic (capacity-limited) growth on the
#'       population \code{x} matrix (having a row per location and a column per
#'       stage/age), and returns the transformed matrix.}
#'   }
#' @references
#'   Beverton, R. J. H., & Holt, S. J. (1957). On the dynamics of exploited
#'   fish populations. \emph{Fisheries Investigations}, 19, 1-533.
#'
#'   García Adeva, J. J., Botha, J. H., & Reynolds, M. (2012). A simulation
#'   modelling approach to forecast establishment and spread of Bactrocera
#'   fruit flies. \emph{Ecological Modelling}, 227, 93–108.
#'   \doi{10.1016/j.ecolmodel.2011.11.026}
#'
#'   Jongejans, E., Skarpaas, O., & Shea, K. (2008). Dispersal, demography and
#'   spatial population models for conservation and control management.
#'   \emph{Perspectives In Plant Ecology Evolution And Systematics}, 9(3–4),
#'   153–170. \doi{10.1016/j.ppees.2007.09.005}
#'
#'   Lefkovitch, L. P. (1965). The Study of Population Growth in Organisms
#'   Grouped by Stages. \emph{Biometrics}, 21(1), 1–18.
#'   \doi{doi.org/10.2307/2528348}
#'
#'   Ricker, W. E. (1958). Handbook of computations for biological statistics
#'   of fish populations. \emph{Bulletin (Fisheries Research Board of Canada)},
#'   119. \url{https://waves-vagues.dfo-mpo.gc.ca/Library/1485.pdf}
#' @include Population.R
#' @export
StagedPopulation <- function(region, growth,
                             capacity = NULL,
                             capacity_area = NULL,
                             capacity_stages = NULL,
                             establish_pr = NULL,
                             incursion_mean = NULL,
                             incursion_stages = NULL, ...) {

  # Build via base class
  self <- Population(region,
                     type = "stage_structured",
                     growth = growth,
                     capacity = capacity,
                     capacity_area = capacity_area,
                     establish_pr = establish_pr,
                     incursion_mean = incursion_mean,
                     class = "StagedPopulation")

  # Number of stages or age classes
  stages = self$get_stages()

  # Get growth (overridden)
  self$get_growth <- function() {
    return(growth)
  }

  # Equivalent growth rate for stage/age matrix (growth)
  growth_r <- Re((sort(eigen(growth, only.values = TRUE)$values,
                       decreasing = TRUE))[1])
  self$get_growth_r <- function() {
    return(growth_r)
  }

  # Growth survivals from attribute when present or lower triangle
  if ("survivals" %in% names(attributes(growth))) {
    survivals <- attr(growth, "survivals")
    if (!is.matrix(survivals) || !is.numeric(survivals) ||
        !all(dim(survivals) == dim(growth))) {
      stop("Growth survivals attribute should be a numeric matrix with the ",
           "same dimensions as the growth (stage/age) matrix.", call. = FALSE)
    }
  } else {
    survivals <- lower.tri(growth, diag = TRUE)*growth
  }

  # Stage/age labels from attribute when present
  if ("labels" %in% names(attributes(growth))) {
    labels <- attr(growth, "labels")
    if (!is.character(labels) || !all(length(labels) == dim(growth))) {
      stop("Stage/age labels attribute should be a character vector ",
           "compatible with the dimensions of the growth (stage/age) matrix.",
           call. = FALSE)
    }
  } else {
    labels <- paste("stage", 1:stages)
    attr(growth, "labels") <- labels
  }

  # Growth reproductions
  reproductions <- growth - survivals
  attr(reproductions, "survivals") <- NULL

  # Validate capacity stages
  if (is.null(capacity_stages)) {
    if (!is.null(capacity)) { # default values
      capacity_stages <- 1:stages
    }
  } else {
    capacity_stages <- unique(capacity_stages)
    if (!is.numeric(capacity_stages) ||
        !all(capacity_stages %in% 1:stages)) {
      stop(paste0("Capacity stages should specify index values between 1 and ",
                  stages, "."), call. = FALSE)
    }
  }

  # Validate incursion stages
  if (!is.null(incursion_stages)) {
    incursion_stages <- unique(incursion_stages)
    if (!is.numeric(incursion_stages) ||
        !all(incursion_stages %in% 1:stages)) {
      stop(paste0("Incursion stages should specify index values between 1 and ",
                  stages, "."), call. = FALSE)
    }
  }

  # Get capacity stages
  self$get_capacity_stages <- function() {
    return(capacity_stages)
  }

  # Set the incursion stages
  self$set_incursion_stages <- function(s) {
    incursion_stages <<- s
  }

  # Make method (extends inherited function from Population class)
  inherited_make <- self$make
  self$make <- function(initial = NULL, current = NULL, incursion = NULL) {

    # Run inherited function for initial or incursion only
    values <- inherited_make(initial = initial, current = NULL,
                             incursion = incursion)

    # Distribute single location values across stages if required
    if (ncol(as.matrix(values)) == 1) {

      # Calculate ratio for stable distribution of stages
      ratio <- abs(Re((eigen(growth)$vectors)[,1]))

      # Select initial stages based on initial age
      if (!is.null(initial) && !is.null(attr(initial, "stages"))) {
        if (!is.null(attr(initial, "age"))) {
          initial_age <- attr(initial, "age") + as.numeric(initial)*0
        } else {
          initial_age <- as.numeric(initial)*0
        }
        unique_ages <- sort(unique(c(0, initial_age)))
        ratio_list <- lapply(unique_ages,
                             function(age) array(0, c(stages, 1)))
        names(ratio_list) <- as.character(unique_ages)
        ratio_list[["0"]][attr(initial, "stages"),] <-
          ratio[attr(initial, "stages")]
        for (age in unique_ages[-1]) {
          ratio_list[[as.character(age)]] <-
            (growth %*% ratio_list[[as.character(age - 1)]] > 0)*ratio
        }
      }

      # Select incursion stages only
      if (!is.null(incursion) && !is.null(incursion_stages)) {
        ratio[which(!1:stages %in% incursion_stages)] <- 0
      }

      # Distribute across stages for occupied indices
      values <- array(values, c(length(values), stages))
      for (i in which(values[,1] > 0)) {
        if (!is.null(initial) && !is.null(attr(initial, "stages"))) {
          ratio <- ratio_list[[as.character(initial_age[i])]]
        }
        values[i,] <- stats::rmultinom(1, size = values[i,], prob = ratio)
      }
    }

    # Combine incursion and current values
    if (!is.null(incursion) && !is.null(current)) {
      if (ncol(as.matrix(current)) == ncol(as.matrix(values))) {
        values <- values + current
      } else {
        stop("Cannot combine incursion with current population array as ",
             "the columns are inconsistent.", call. = FALSE)
      }
    }

    # Label columns
    colnames(values) <- labels

    return(values)
  }

  # Look-up table of matrix multipliers for equivalent growth rates
  if (is.numeric(capacity)) { # capacity-limited

    # Maximum multiplier when 100% survival
    max_mult <- max(1/max(colSums(survivals), 0.001), 1)

    # Construct look-up table
    r_mult <- data.frame(r = NA, mult = (1:trunc(max_mult*1000))/1000)
    r_mult$r <- sapply(r_mult$mult, function(m) {
      Re((sort(eigen(growth*m, only.values = TRUE)$values,
               decreasing = TRUE))[1])
    })
  }

  # Grow method - override for logistic (capacity-limited) growth
  self$grow <- function(x) {

    # Indices of occupied locations
    indices <- which(rowSums(x) > 0)

    if (length(indices)) { # grow occupied populations

      # Calculate logistic growth rates and their equivalent matrix multipliers
      if (is.numeric(capacity)) { # capacity-limited

        # Remove populations at locations having zero capacity
        if (any(capacity[indices] <= 0)) {
          zero_idx <- indices[which(capacity[indices] <= 0)]
          x[zero_idx,] <- 0
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
            area_capacity <- capacity*pi*capacity_radius^2/capacity_area

            # Calculate capacity-limited growth rate
            r <- exp(log(growth_r)*(1 - sum(x[capacity_stages])/area_capacity))

          } else if (is.numeric(capacity_area) &&
                     is.numeric(region$get_max_implicit_area())) {

            # Calculate capacity of maximum area
            area_capacity <- (capacity*region$get_max_implicit_area()/
                                capacity_area)

            # Calculate capacity-limited growth rate
            r <- exp(log(growth_r)*(1 - sum(x[capacity_stages])/area_capacity))

          } else {
            r <- growth_r # unlimited
          }

        } else {

          # Calculate capacity-limited growth rates for each occupied location
          r <- exp(log(growth_r)*(1 - (rowSums(x[indices, capacity_stages,
                                                 drop = FALSE])/
                                         capacity[indices])))
        }

        # Look-up stage/age matrix multiplier for each occupied location
        mult <- sapply(r, function(r) {
          r_mult$mult[which.min(abs(r_mult$r - r))]})

      } else {
        mult <- rep(1, length(indices))
      }

      # Sample the new population values (reproduction and survival)
      if (length(indices)) {
        new_x <- array(0L, c(length(indices), stages))
        for (stage in 1:stages) {

          # Sample stage reproduction via a Poisson distribution
          new_x <- new_x + stats::rpois(
            length(indices)*stages,
            (x[indices, stage]*
               t(reproductions[, rep(stage, length(indices))])*mult))

          # Sample stage survival via a binomial distribution
          stage_surv <- stats::rbinom(length(indices), x[indices, stage],
                                      rep(sum(survivals[, stage]),
                                          length(indices))*mult)

          # Distribute survivals across stages via multinomial sampling
          new_x <- new_x + t(sapply(1:length(indices), function(i) {
            if (any(survivals[, stage]*mult[i] > 0)) {
              stats::rmultinom(1, size = stage_surv[i],
                               prob = survivals[, stage]*mult[i])
            } else {
              matrix(0, nrow = nrow(survivals), ncol = 1)
            }
          }))
        }

        # Update populations at occupied locations
        x[indices,] <- new_x
      }
    }

    return(x)
  }

  return(self)
}
