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
#' @seealso \code{\link{fit.srm.poireg}}
NULL

#' @rdname sGLM
sGLM <- R6::R6Class("sGLM",
  private = list(
    linkfun = NA,
    x = 1,
    clear.params = function(nm) {
      self$params <- numeric(0)
      private$x <- 1
    },
    add.params = function(nm, params) {
      if (length(params) != 0) {
        self$params <- c(self$params, params)
        self$params.position[[nm]] <- private$x:(private$x+length(params)-1)
        private$x <- private$x + length(params)
      }
      else {
        self$params.position[[nm]] <- c()
      }
    },
    get.params = function(params, nm) {
      params[self$params.position[[nm]]]
    }
  ),
  public = list(
    name = NA,
    names = NULL,
    srms = NULL,
    params = NA,
    params.position = list(),
    df = NULL,
    data = NULL,
    print = function(digits = max(3, getOption("digits") - 3), ...) {
      print(self$data)
      cat(gettextf("\nLink function: %s\n", private$linkfun))
      print(self$coefficients())
      for (nm in self$names) {
        cat(gettextf("\n%s\n", nm))
        print(self$srms[[nm]])
      }
    },
    omega = function() {
      result <- sapply(self$names, function(nm) self$params[self$params.position[[nm]]][1L])
      names(result) <- self$names
      result
    },
    coefficients = function() {
      result <- self$params[self$params.position$coefficients]
      names(result) <- colnames(self$data$metrics)
      result
    },
    initialize = function(srms, names = NULL, coefficients = c()) {
      if (is.null(names)) {
        if (is.null(names(srms))) {
          stop("names or names of srms list should be needed.")
        }
        else {
          self$names <- names(srms)
          self$srms <- srms
        }
      }
      else {
        if (is.null(names(srms))) {
          m <- length(names)
          if (length(srms) < m) {
            stop("length of srms should be greater than or equal to m.")
          }
          else {
            self$names <- names
            self$srms <- lapply(1:m, function(i) srms[[i]])
            names(self$srms) <- self$names
          }
        }
        else {
          stopifnot(all(names %in% names(srms)))
          self$names <- names
          self$srms <- lapply(self$names, function(nm) srms[[nm]])
          names(self$srms) <- self$names
        }
      }

      private$clear.params()
      for (nm in self$names) {
        private$add.params(nm, self$srms[[nm]]$params)
      }
      private$add.params("coefficients", coefficients)

      self$df <- sum(sapply(srms, function(m) m$df-1L))  + length(self$params.position$coefficients)
    },
    init_params = function(data) {
      stopifnot(all(self$names %in% data$names))
      private$clear.params()
      for (nm in self$names) {
        self$srms[[nm]]$init_params(self$srms[[nm]]$data)
        private$add.params(nm, self$srms[[nm]]$params)
      }
      private$add.params("coefficients", numeric(data$nmetrics))
      self$df <- sum(sapply(self$srms, function(m) m$df-1L))  + length(self$params.position$coefficients)

      self$set_params(self$em(self$params, data)$param)
    },
    set_params = function(params) {
      self$params <- params
      for (nm in self$names) {
        self$srms[[nm]]$set_params(private$get.params(params, nm))
      }
    },
    set_data = function(data) { self$data <- data },
    em = function(params, data, ...) {
      newparams <- params
      result <- lapply(self$names, function(nm) self$srms[[nm]]$em(private$get.params(params, nm), self$srms[[nm]]$data, ...))
      names(result) <- self$names

      llf <- sum(sapply(result, function(r) r$llf))
      total <- sapply(result, function(r) r$total)

      for (nm in self$names) {
        newparams[self$params.position[[nm]]] <- result[[nm]]$param
      }

      wopt <- getOption("warn")
      options(warn = -1)
      result <- glm.fit(x=data$metrics, y=total, family=poisson(link=private$linkfun), ...)
      options(warn = wopt)
      for (nm in self$names) {
        newparams[self$params.position[[nm]]][1L] <- result$fitted.values[nm]
      }
      newparams[self$params.position$coefficients] <- result$coefficients

      pdiff <- newparams - params
      list(param=newparams, pdiff=pdiff, llf=llf, total=NULL)
    },
    llf = function(data) {
      sum(sapply(self$names, function(nm) self$srms[[nm]]$llf(self$srms[[nm]]$data)))
    }
  )
)

#' @rdname sGLM
#' @export
sGLM.log <- R6::R6Class("sGLM.log",
  inherit = sGLM,
  private = list(
    linkfun = "log"
  ),
  public = list(
    name = "sGLM.log"
  )
)

#' @rdname sGLM
#' @export
sGLM.identity <- R6::R6Class("sGLM.identity",
  inherit = sGLM,
  private = list(
    linkfun = "identity"
  ),
  public = list(
    name = "sGLM.identity"
  )
)
