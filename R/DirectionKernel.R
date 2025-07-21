#' Direction Kernel functions generation class builder
#'
#' Builds a class for generating direction functions (or kernels) for
#' calculating the (relative) probability of dispersal, used by
#' \code{Dispersal} models and inherited class objects. The kernel functions
#' include the Beta distribution function and a function for defining the
#' probability of directions from a table of values. Directions may be defined
#' in degrees (0-360) using trigonometric angles or navigational directions.
#' The orientation of directions may also be specified to describe
#' whether direction probabilities refer to the likelihood of dispersing to or
#' from specified directions.
#'
#' @param direction_type Specifies whether directions (0-360 degrees) are
#'   defined using \code{"trigonometric"} angles (i.e. anticlockwise from East
#'   at 0 degrees) or \code{"navigational"} directions (i.e. clockwise from
#'   North at 0 degrees). Default is \code{"trigonometric"}.
#' @param orientation Specifies whether direction probabilities refer to the
#'   likelihood of dispersing \code{"to"} (e.g. bird migration) or
#'   \code{"from"} (e.g. wind) specified directions. Default is \code{"to"}.
#' @param multiplier Optional numeric multiplier for scaling the (relative)
#'   probabilities returned by the kernel functions.
#' @param ... Additional parameters.
#' @return A \code{DirectionKernel} class object (list) containing functions for
#'   generating direction kernel functions for calculating the (relative)
#'   probability of dispersal. The generated (returned) functions of form
#'   \code{function(x)}, where \code{x} is a vector of directions (0-360
#'   degrees), which return the (relative) probability for each value in
#'   \code{x}:
#'   \describe{
#'     \item{\code{get_beta_function(alpha, beta, shift = 0)}}{Get a Beta
#'       direction kernel function specified with distribution parameters
#'       \code{alpha} and \code{beta} scaled for directions \code{x} 0-360
#'       degrees, as well as an optional parameter to \code{shift} (default is
#'       0) the Beta distribution along the \code{x} axis, wrapping the
#'       distribution from 360 degrees back to zero.}
#'     \item{\code{get_lookup_function(table)}}{Get a table look-up kernel
#'       function specified with parameter \code{table}, a two-column data
#'       frame (or matrix) of direction (0-360 degrees) \code{x} values (first
#'       column) and corresponding (relative) probability values (second
#'       column). Look-up functionality performs linear interpolation between
#'       discrete values.}
#'   }
#' @references
#'   Carrasco, L. R., Harwood, T. D., Toepfer, S., MacLeod, A., Levay, N.,
#'   Kiss, J., Baker, R. H. A., Mumford, J. D., & Knight, J. D. (2010).
#'   Dispersal kernels of the invasive alien western corn rootworm and the
#'   effectiveness of buffer zones in eradication programmes in Europe.
#'   \emph{Annals of Applied Biology}, 156(1), 63–77.
#'   \doi{10.1111/j.1744-7348.2009.00363.x}
#'
#'   Jongejans, E., Skarpaas, O., & Shea, K. (2008). Dispersal, demography and
#'   spatial population models for conservation and control management.
#'   \emph{Perspectives In Plant Ecology Evolution And Systematics}, 9(3–4),
#'   153–170. \doi{10.1016/j.ppees.2007.09.005}
#'
#'   Shaw, M. W. (1995). Simulation of Population Expansion and Spatial Pattern
#'   when Individual Dispersal Distributions do not Decline Exponentially with
#'   Distance. \emph{Proceedings B: Biological Sciences}, 259(1356), 243–248.
#'   \doi{10.1098/rspb.1995.0036}
#' @export
DirectionKernel <- function(direction_type = c("trigonometric",
                                               "navigational"),
                            orientation = c("to", "from"),
                            multiplier = NULL, ...) {

  # Match direction type and orientation
  direction_type <- match.arg(direction_type)
  orientation <- match.arg(orientation)

  # Default multiplier is 1
  if (is.null(multiplier)) {
    multiplier <- 1
  }

  # Build via base class
  self <- Kernels(multiplier = multiplier,
                  class = "DirectionKernel")
  super <- list(get_lookup_function = self$get_lookup_function)

  # Get a Beta direction kernel function
  self$get_beta_function <- function(alpha, beta, shift = 0) {
    if (shift > 360) {
      stop("The Beta direction kernel shift parameter must be <= 360 degrees.",
           call. = FALSE)
    }
    return(
      function(x) {
        if (direction_type == "navigational") {
          x <- (x <= 90)*(90 - x) + (x > 90)*(450 - x)
        }
        if (orientation == "from") {
          x <- (x < 180)*(x + 180) + (x >= 180)*(x - 180)
        }
        x <- (x >= shift)*(x - shift) + (x < shift)*(x - shift + 360)
        return(multiplier*(alpha/(alpha + beta))*
                 stats::dbeta(x/360, shape1 = alpha, shape2 = beta))
      }
    )
  }

  # Get a table look-up direction kernel function
  self$get_lookup_function <- function(table) {
    if (direction_type == "navigational") {
      x <- table[,1]
      table[,1] <- (x <= 90)*(90 - x) + (x > 90)*(450 - x)
    }
    if (orientation == "from") {
      x <- table[,1]
      table[,1] <- (x < 180)*(x + 180) + (x >= 180)*(x - 180)
    }
    if (any(c(0, 360) %in% table[,1])) {
      if (0 %in% table[,1] && !(360 %in% table[,1])) {
        table <- rbind(table, c(360, table[which(table[,1] == 0),2]))
      }
      if (360 %in% table[,1] && !(0 %in% table[,1])) {
        table <- rbind(table, c(0, table[which(table[,1] == 360),2]))
      }
    } else {
      i_min <- which.min(table[,1])
      i_max <- which.max(table[,1])
      prob_0 <- (
        table[i_max, 2] + (table[i_min, 2] - table[i_max, 2])*
          (360 - table[i_max, 1])/(360 + table[i_min, 1] - table[i_max, 1]))
      table <- rbind(c(0, prob_0), table, c(360, prob_0))
    }
    return(super$get_lookup_function(table))
  }

  return(self)
}
