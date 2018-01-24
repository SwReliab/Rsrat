
#' @rdname sGLM
sGLM.penalized <- R6::R6Class("sGLM.penalized",
  inherit = sGLM,
  public = list(
    lambda = NULL,
    K = NULL,
    set_penalized = function(lambda, K = NULL) {
      self$lambda <- lambda
      self$K <- K
    },
    em = function(params, data, ...) {
      newparams <- params

      result <- lapply(self$names, function(nm) self$srms[[nm]]$em(private$get.params(params, nm), self$srms[[nm]]$data, ...))
      names(result) <- self$names

      llf <- sum(sapply(result, function(r) r$llf))
      total <- sapply(result, function(r) r$total)
      for (nm in self$names) {
        newparams[self$params.position[[nm]]] <- result[[nm]]$param
      }

      result <- pnglm.fit(x=data$metrics, y=total, family=poisson(link=private$linkfun),
        lambda=self$lambda, K=self$K)
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
sGLM.penalized.log <- R6::R6Class("sGLM.penalized.log",
  inherit = sGLM.penalized,
  private = list(
    linkfun = "log"
  ),
  public = list(
    name = "sGLM.penalized.log"
  )
)

#' @rdname sGLM
#' @export
sGLM.penalized.identity <- R6::R6Class("sGLM.penalized.identity",
  inherit = sGLM.penalized,
  private = list(
    linkfun = "identity"
  ),
  public = list(
    name = "sGLM.penalized.identity"
  )
)
