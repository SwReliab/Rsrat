#' Class for software reliability model with s-metrics
#'
#' @docType class
#' @name sGLM
#' @return Object of \code{\link{R6Class}} with methods for software reliability model with s-metrics.
#' @format \code{\link{R6Class}} object.
#' @field srms A list for srms.
#' @field params A numeric vector for the regression coefficients.
#' @field df An integer for the degrees of freedom of the model (total).
#' @field data Data to esimate parameters.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{print()}}{This method prints model parameters.}
#'   \item{\code{omega()}}{This method returns a vector of the number of total faults.}
#'   \item{\code{coefficients()}}{This method returns a vector of the coefficients for s-metrics.}
#'   \item{\code{init_params(data)}}{This method changes the model parameters based on a given data.
#'          This is used to set the initial value for the fitting algorithm.}
#'   \item{\code{set_params(params)}}{This method sets the model parameters.}
#'   \item{\code{set_data(data)}}{This method sets data.}
#'   \item{\code{em(params, data)}}{This method returns a list with an updated parameter vector (param),
#'          absolute difference of parameter vector (pdiff),
#'          log-likelihood function for a given parameter vector (llf),
#'          the number of total faults (total) via EM algorithm for a given data.}
#'   \item{\code{llf(data)}}{This method returns the log-likelihood function for a given data.}
#' }
#' @seealso \code{\link{fit.srm.logit}}
NULL

#' @rdname sGLM
sGLM <- R6::R6Class("sGLM",
  private = list(
    linkfun = NA
  ),
  public = list(
    srms = NA,
    params = NA,
    df = NA,
    data = NA,
    print = function(digits = max(3, getOption("digits") - 3), ...) {
      cat(gettextf("Link function: %s\n", private$linkfun))
      for (m in srms) {
        print(m)
      }
    },
    omega = function() { self$params[1L] },
    coefficients = function() { self$params[2L:length(self$params)] },
    initialize = function(omega = 1, coefficients = c(1)) {
      self$params <- c(omega, coefficients)
      self$df <- length(self$params)
    },
    init_params = function(data) {
      self$params <- numeric(1L + data$nmetrics)
      self$params[1] <- data$total + 1.0
      self$df <- length(self$params)
    },
    set_params = function(params) { self$params <- params },
    set_data = function(data) { self$data <- data },
    em = function(params, data, ...) {
      omega <- params[1]
      coefficients <- params[2L:length(params)]
      family <- binomial(link = private$linkfun)
      eta <- data$metrics %*% coefficients + data$offset
      mu <- family$linkinv(eta)
      residual <- omega * prod(1-mu)
      total <- sum(data$fault) + residual
      rfault <- total - cumsum(data$fault)
      wopt <- getOption("warn")
      options(warn = -1)
      result <- glm.fit(data$metrics, cbind(data$fault, rfault),
        family=binomial(link=private$linkfun), offset=data$offset, ...)
      options(warn = wopt)
      newparams <- c(total, result$coefficients)
      names(newparams) <- c("omega", names(result$coefficients))
      pdiff <- abs(params - newparams)
      llf <- self$llf(data, omega=omega, mu=mu)
      list(param=newparams, pdiff=pdiff, llf=llf, total=total)
    },
    llf = function(data, fault, omega, mu) {
      if (missing(omega)) {
        omega <- self$omega()
      }
      if (missing(mu)) {
        family <- binomial(link = private$linkfun)
        eta <- data$metrics %*% self$coefficients() + data$offset
        mu <- family$linkinv(eta)
      }
      if (missing(fault)) {
        fault <- data$fault
      }
      nonzeron <- fault != 0
      rfault <- sum(fault) - cumsum(fault)
      nonzeror <- rfault != 0
      sum((fault * log(mu))[nonzeron]) + sum((rfault * log(1-mu))[nonzeror]) -
        sum(lgamma(fault+1)) + sum(fault) * log(omega) - omega * (1 - prod(1-mu))
    }
  )
)

#' @rdname dGLM
#' @export
dGLM.logit <- R6::R6Class("dGLM.logit",
  inherit = dGLM,
  private = list(
    linkfun = "logit"
  ),
  public = list(
    name = "dGLM.logit"
  )
)

#' @rdname dGLM
#' @export
dGLM.probit <- R6::R6Class("dGLM.probit",
  inherit = dGLM,
  private = list(
    linkfun = "probit"
  ),
  public = list(
    name = "dGLM.probit"
  )
)

#' @rdname dGLM
#' @export
dGLM.cloglog <- R6::R6Class("dGLM.cloglog",
  inherit = dGLM,
  private = list(
    linkfun = "cloglog"
  ),
  public = list(
    name = "dGLM.cloglog"
  )
)
