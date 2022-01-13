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
#'   matrix, including those on the diagonal.
#' @param capacity A vector of carrying capacity values of the invasive species
#'   at each location specified by the \code{region}. Default is \code{NULL}
#'   for when growth is not capacity-limited.
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
#' @param incursion_values Defines how incursion locations are populated.
#'   A matrix of fixed values with a row for each region location, or a single
#'   row applicable for all locations, and with columns for each stage/age.
#'   Alternatively, a list of matrices (as before), labelled \code{min} and
#'   \code{max}, be used to specifying ranges of values for uniform random
#'   sampling of values.
#' @param ... Additional parameters.
#' @return A \code{StagedPopulation} class object (list) containing functions
#'   for accessing attributes and simulating growth:
#'   \describe{
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
#'     \item{\code{make(initial, current, incursion)}}{Make a population matrix
#'       (with a row per location and a column per stage/age) via using vectors
#'       or matrices of the \code{initial} or \code{current} and
#'       \code{incursion} population at each region location.}
#'     \item{\code{grow(x)}}{Performs logistic (capacity-limited) growth on the
#'       population \code{x} matrix (having a row per location and a column per
#'       stage/age), and returns the transformed matrix.}
#'   }
#' @include Population.R
#' @export
StagedPopulation <- function(region, growth,
                             capacity = NULL,
                             capacity_stages = NULL,
                             establish_pr = NULL,
                             incursion_values = NULL, ...) {

  # Build via base class
  self <- Population(region,
                     type = "stage_structured",
                     growth = growth,
                     capacity = capacity,
                     establish_pr = establish_pr,
                     incursion_values = incursion_values,
                     class = "StagedPopulation")

  # Number of stages or age classes
  stages = self$get_stages()

  # Equivalent growth rate for stage/age matrix (growth)
  growth_r <- Re((eigen(growth, only.values = TRUE)$values)[1])
  get_growth_r <- function() {
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

  # Get capacity stages
  self$get_capacity_stages <- function() {
    return(capacity_stages)
  }

  # Make method (extends inherited function from Population class)
  inherited_make <- self$make
  self$make <- function(initial = NULL, current = NULL, incursion = NULL) {

    # Run inherited function
    values <- inherited_make(initial = initial, current = current,
                             incursion = incursion)

    # Distribute single location values across stages if required
    if (ncol(as.matrix(values)) == 1) {

      # Calculate ratio for stable distribution of stages
      ratio <- Re((eigen(growth)$vectors)[,1])

      # Distribute across stages for occupied indices
      values <- array(values, c(length(values), stages))
      for (i in which(values[,1] > 0)) {
        values[i,] <- stats::rmultinom(1, size = values[i,], prob = ratio)
      }
    }

    return(values)
  }

  # Look-up table of matrix multipliers for equivalent growth rates
  if (is.numeric(capacity)) { # capacity-limited

    # Maximum multiplier when 100% survival
    max_mult <- max(1/max(colSums(survivals), 0.001), 1)

    # Construct look-up table
    r_mult <- data.frame(r = NA, mult = (1:trunc(max_mult*1000))/1000)
    r_mult$r <- sapply(r_mult$mult, function(m) {
      Re((eigen(growth*m, only.values = TRUE)$values)[1])
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

        # Calculate capacity-limited growth rates for each occupied location
        r <- exp(log(growth_r)*(1 - (rowSums(x[indices, capacity_stages])/
                                       capacity[indices])))

        # Look-up stage/age matrix multiplier for each occupied location
        mult <- sapply(r, function(r) {
          r_mult$mult[which.min(abs(r_mult$r - r))]})

      } else {
        mult <- rep(1, length(indices))
      }

      # Sample the new population values (reproduction and survival)
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
          stats::rmultinom(1, size = stage_surv[i],
                           prob = survivals[, stage]*mult[i])
        }))
      }

      # Update populations at occupied locations
      x[indices,] <- new_x
    }

    return(x)
  }

  return(self)
}
