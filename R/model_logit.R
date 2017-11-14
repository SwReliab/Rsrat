#' Class for NHPP-based software reliability model with d-metrics
#'
#' @docType class
#' @name dGLM
#' @return Object of \code{\link{R6Class}} with methods for NHPP-based software reliability model with d-metrics.
#' @format \code{\link{R6Class}} object.
#' @field name A character string for the name of model.
#' @field params A numeric vector for the model parameters.
#' @field df An integer for the degrees of freedom of the model.
#' @field data Data to esimate parameters.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{print()}}{This method prints model parameters.}
#'   \item{\code{omega()}}{This method returns the number of total faults.}
#'   \item{\code{coefficients()}}{This method returns a vector for the coefficients.}
#'   \item{\code{mvf(t, data = NULL)}}{This method returns the mean value function at time t.
#'          The d-metrics is given from \code{data}. If \code{data} is NULL, the d-metrics for the estimation is used.}
#'   \item{\code{dmvf(t, data = NULL)}}{This method returns the intensity function at time t.
#'          The d-metrics is given from \code{data}. If \code{data} is NULL, the d-metrics for the estimation is used.}
#'   \item{\code{residual(t, data = NULL)}}{This method returns the expected residual number of faults at time t.
#'          The d-metrics is given from \code{data}. If \code{data} is NULL, the d-metrics for the estimation is used.}
#'   \item{\code{ffp(t, data = NULL)}}{This method returns the fault-free probability at time t.
#'          The d-metrics is given from \code{data}. If \code{data} is NULL, the d-metrics for the estimation is used.}
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

#' @rdname dGLM
dGLM <- R6::R6Class("dGLM",
  private = list(
    linkfun = NA
  ),
  public = list(
    name = NA,
    params = NA,
    df = NA,
    data = NA,
    print = function(digits = max(3, getOption("digits") - 3), ...) {
      cat(gettextf("Link function: %s\n", private$linkfun))
      print.default(format(self$params, digits = digits), print.gap = 2, quote = FALSE)
    },
    omega = function() { self$params[1L] },
    coefficients = function() { self$params[2L:length(self$params)] },
    mvf = function(t, data = NULL) {
      if (is.null(data)) {
        metrics <- self$data$metrics
        offset <- self$data$offset
      } else {
        metrics <- data$metrics
        offset <- data$offset
      }
      mname <- names(self$coefficients())
      if (!all(mname %in% colnames(metrics))) {
        warning("colnames(metrics) do not match to names(coefficients).")
        metrics <- metrics[,1:length(mname)]
        colnames(metrics) <- mname
      }
      else {
        metrics <- metrics[,mname]
      }
      family <- binomial(link = private$linkfun)
      result <- if (length(mname) == 1L) {
        sapply(t, function(t0) {
            eta <- metrics[1L:t0] * self$coefficients() + offset[1L:t0]
            mu <- family$linkinv(eta)
            self$omega() * (1 - prod(1-mu))
          }
        )
      }
      else {
        sapply(t, function(t0) {
            eta <- metrics[1L:t0,] %*% self$coefficients() + offset[1L:t0]
            mu <- family$linkinv(eta)
            self$omega() * (1 - prod(1-mu))
          }
        )
      }
      names(result) <- NULL
      result
    },
    dmvf = function(t, data = NULL) {
      result <- self$mvf(t, data)
      c(result[1], diff(result))
    },
    residual = function(t, data = NULL) {
      if (is.null(data)) {
        metrics <- self$data$metrics
        offset <- self$data$offset
      } else {
        metrics <- data$metrics
        offset <- data$offset
      }
      mname <- names(self$coefficients())
      if (!all(mname %in% colnames(metrics))) {
        warning("colnames(metrics) do not match to names(coefficients).")
        metrics <- metrics[,1:length(mname)]
        colnames(metrics) <- mname
      }
      else {
        metrics <- metrics[,mname]
      }
      family <- binomial(link = private$linkfun)
      result <- if (length(mname) == 1L) {
        sapply(t, function(t0) {
            eta <- metrics[1:t0] * self$coefficients() + offset[1L:t0]
            mu <- family$linkinv(eta)
            self$omega() * prod(1-mu)
          }
        )
      }
      else {
        sapply(t, function(t0) {
            eta <- metrics[1:t0,] %*% self$coefficients() + offset[1L:t0]
            mu <- family$linkinv(eta)
            self$omega() * prod(1-mu)
          }
        )
      }
      names(result) <- NULL
      result
    },
    ffp = function(t, data = NULL) { exp(-self$residual(t, data)) },
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
