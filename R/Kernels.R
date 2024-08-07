#' Kernel functions generation class builder
#'
#' Builds a class for generating distance or direction functions (or kernels)
#' for calculating the (relative) probability of dispersal, used by
#' \code{Dispersal} models and inherited class objects. The kernel functions
#' include Beta, Cauchy, (negative) Exponential, Gaussian (normal), Lognormal,
#' and Weibull distributions. Also included is a function for defining the
#' probability of distances or directions from a table of values.
#'
#' @param multiplier Optional numeric multiplier for scaling the (relative)
#'   probabilities returned by the kernel functions.
#' @param ... Additional parameters.
#' @return A \code{Kernel} class object (list) containing functions for
#'   generating distance or direction kernel functions for calculating the
#'   (relative) probability of dispersal. The generated (returned) functions
#'   of form \code{function(x)}, where \code{x} is a vector of distances (in m)
#'   or directions (0-360 degrees), which return the (relative) probability for
#'   each value in \code{x}:
#'   \describe{
#'     \item{\code{get_beta_function(alpha, beta, upper = 1, shift = 0)}}{Get a
#'       Beta kernel function specified with distribution parameters
#'       \code{alpha} and \code{beta}, as well as optional parameters for the
#'       \code{upper} limit (default is 1) of \code{x}, and to \code{shift}
#'       (default is 0) the Beta distribution along the \code{x} axis, wrapping
#'       the distribution from the \code{upper} limit back to zero, which is
#'       useful for defining direction kernels.}
#'     \item{\code{get_cauchy_function(scale)}}{Get a (half) Cauchy kernel
#'       function specified with distribution parameter \code{scale} (for
#'       scaling \code{x}).}
#'     \item{\code{get_exp_function(mean)}}{Get a (negative) Exponential kernel
#'       function specified with distribution parameter \code{mean}.}
#'     \item{\code{get_gaussian_function(sd)}}{Get a (half) Gaussian (normal)
#'       kernel function specified with distribution parameter \code{sd}
#'       (standard deviation).}
#'     \item{\code{get_lognormal_function(mean, sd)}}{Get a Lognormal kernel
#'       function specified with distribution parameters \code{mean} and
#'       \code{sd} (standard deviation).}
#'     \item{\code{get_weibull_function(shape, scale)}}{Get a Weibull kernel
#'       function specified with distribution parameters \code{shape} and
#'       \code{scale} (for scaling \code{x}).}
#'     \item{\code{get_lookup_function(table)}}{Get a table look-up kernel
#'       function specified with parameter \code{table}, a two-column data
#'       frame (or matrix) of \code{x} values (first column) and corresponding
#'       (relative) probability values (second column). Look-up functionality
#'       performs linear interpolation between discrete values.}
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
Kernels <- function(multiplier = NULL, ...) {

  # Check/set multiplier
  if (!is.null(multiplier) && (!is.numeric(multiplier) || multiplier <= 0)) {
    stop("The multiplier parameter must be numeric and > 0.", call. = FALSE)
  }
  if (is.null(multiplier)) {
    multiplier <- 1
  }

  # Create a class structure
  self <- structure(list(), class = "Kernels")

  # Get a Beta kernel function
  self$get_beta_function <- function(alpha, beta, upper = 1, shift = 0) {
    if (shift > upper) {
      stop("The Beta kernel shift parameter must be <= the upper parameter.",
           call. = FALSE)
    }
    return(
      function(x) {
        x <- (x >= shift)*(x - shift) + (x < shift)*(x - shift + upper)
        return(multiplier*(alpha/(alpha + beta))*
                 stats::dbeta(x/upper, shape1 = alpha, shape2 = beta))
      }
    )
  }

  # Get a (half) Cauchy kernel function
  self$get_cauchy_function <- function(scale) {
    return(
      function(x) {
        return(2*multiplier*scale*stats::dcauchy(x, scale = scale))
      }
    )
  }

  # Get a (negative) Exponential kernel function
  self$get_exp_function <- function(mean) {
    return(
      function(x) {
        return(multiplier*mean*stats::dexp(x, rate = 1/mean))
      }
    )
  }

  # Get a (half) Gaussian (normal) kernel function
  self$get_gaussian_function <- function(sd) {
    return(
      function(x) {
        return(2*multiplier*sd*stats::dnorm(x, sd = sd))
      }
    )
  }

  # Get a Lognormal kernel function
  self$get_lognormal_function <- function(mean, sd) {
    meanlog <- log(mean^2/sqrt(mean^2 + sd^2))
    sdlog <- sqrt(log(1 + sd^2/mean^2))
    return(
      function(x) {
        return(multiplier*mean*stats::dlnorm(x, meanlog = meanlog,
                                             sdlog = sdlog))
      }
    )
  }

  # Get a Weibull kernel function
  self$get_weibull_function <- function(shape, scale) {
    return(
      function(x) {
        return(multiplier*scale*stats::dweibull(x, shape = shape,
                                                scale = scale))
      }
    )
  }

  # Get a table look-up kernel function
  self$get_lookup_function <- function(table) {
    table <- table[order(table[,1]),]
    return(
      function(x) {

        # Look-up x values in table
        idx <- pmin(as.numeric(cut(x, breaks = c(-Inf, table[,1], Inf))),
                    nrow(table))

        # Get range of x and probability values
        d <- as.data.frame(cbind(x, table[pmax(idx - 1, 1),], table[idx,]))
        names(d)[2:5] <- c("x1", "p1", "x2", "p2")

        # Interpolate probability values using x range when in-between
        d$p <- d$p2
        idx <- which(x >= d$x1 & x < d$x2)
        d$p[idx] <- (d$p1 + (d$x - d$x1)/(d$x2 - d$x1)*(d$p2 - d$p1))[idx]

        return(multiplier*d$p)
      }
    )
  }

  return(self)
}
