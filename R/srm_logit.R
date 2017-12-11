#' Regression-based software reliability model with d-metrics
#'
#' Estimate model parameters for d-metrics-based software reliability model.
#'
#' The control argument is a list that can supply any of the following components:
#' \describe{
#'   \item{maxiter}{An integer for the maximum number of iterations in the fitting algorithm.}
#'   \item{reltol}{A numeric value. The algorithm stops if the relative error is
#' less than \emph{reltol} and the absolute error is less than \emph{abstol}.}
#'   \item{abstol}{A numeric value. The algorithm stops if the relative error is
#' less than \emph{reltol} and the absolute error is less than \emph{abstol}.}
#'   \item{stopcond}{A character string. \emph{stopcond} gives the criterion
#' for the stop condition of the algorithm. Either llf or parameter is selected.}
#'   \item{trace}{A logical. If TRUE, the intermediate parameters are printed.}
#'   \item{printsteps}{An integer for print.}
#' }
#' The linkfun argument can take the following strings:
#' \describe{
#'   \item{logit}{A logit function.}
#'   \item{probit}{A probit function.}
#'   \item{cloglog}{A complementary log-log function.}
#' }
#'
#' @param formula An object of class formula. A symbolic description of the
#' model to be fitted. The output variable should be the column for the number
#' of faults.
#' @param data A dataframe for d-metrics and the number of faults.
#' @param linkfun A character string indicating a linkfunction. See Details.
#' @param offset An integer. This can be used to specify an a priori known
#' component to be included in the linear predictor during fitting. This should
#' be NULL or a numeric vector of length equal to the number of cases.
#' @param control A list of control parameters. See Details.
#' @param ... Other parameters.
#' @return A list with components;
#' \item{initial}{A vector for initial parameters.}
#' \item{srm}{A class of NHPP. The SRM with the estiamted parameters.}
#' \item{llf}{A numeric value for the maximum log-likelihood function.}
#' \item{df}{An integer for degrees of freedom.}
#' \item{convergence}{A boolean meaning the alorigthm is converged or not.}
#' \item{iter}{An integer for the number of iterations.}
#' \item{aerror}{A numeric value for absolute error.}
#' \item{rerror}{A numeric value for relative error.}
#' \item{data}{The data used.}
#' \item{linkfun}{The linkfunction used.}
#' \item{formula}{The formula supplied.}
#' \item{ctime}{A numeric value for computation time.}
#' \item{terms}{The terms object used.}
#' \item{call}{The method call.}
#' @export

fit.srm.logit <- function(formula, data, linkfun = "logit", offset = NULL, control = list(), ...) {
  call <- match.call()
  con <- srm.logit.options()
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
  warning("unknown names in control: ", paste(noNms, collapse = ", "))

  ## init parameters
  ldata <- faultdata.dmet(formula, data, offset)
  if (ldata$nmetrics > nrow(ldata$metrics)) {
    warning("too many parameters: the number of parameters might be smaller than the number of data.")
  }

  model <- switch(linkfun,
    "logit"=dGLM.logit$new(),
    "probit"=dGLM.probit$new(),
    "cloglog"=dGLM.cloglog$new(),
    NA
  )

  tres <- system.time(result <- emfit(model, ldata, initialize = TRUE,
    maxiter = con$maxiter, reltol = con$reltol, abstol = con$abstol,
    stopcond = con$stopcond, trace=con$trace, printsteps=con$printsteps))
  result <- c(result,
              list(
                data=data,
                linkfun=linkfun,
                formula=formula,
                ctime=tres[1],
                terms = attr(model.frame(formula, data), "terms"),
                call=call))
  class(result) <- "srm.logit.result"
  result
}

#' Options for fit.srm.logit
#'
#' Generate a list of option values.
#'
#' @return A list of options.
#' @rdname fit.srm.logit
#' @export

srm.logit.options <- function() {
  list(maxiter = 10000,
    reltol = sqrt(.Machine$double.eps),
    abstol = 1.0e+200,
    stopcond = "llf",
    trace = FALSE,
    printsteps = 50)
}

#' Print a fit result of fit.srm.logit (S3 method)
#'
#' Print the fitting result.
#'
#' @param x An object of srm.logit.result.
#' @param ... Other parameters
#' @param digits The minimum number of significant digits.
#' @details This function uses the print for dGLM.
#' @export

print.srm.logit.result <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  if (!is.null(x$srm)) {
    print(x$srm, digits=digits, ...)
  }
  if (!is.null(x$llf))
  cat("Maximum LLF:", x$llf, "\n")
  if (!is.null(x$df))
  cat("AIC:", -2*x$llf+2*x$df, "\n")
  if (!is.null(x$convergence))
  cat("Convergence:", x$convergence, "\n\n")
  invisible(x)
}

#' AIC for a fit result of fit.srm.logit (S3 method)
#'
#' Return AIC for the fitting result.
#'
#' @param fit An object of srm.logit.result.
#' @param ... Other parameters
#' @param scale used in the definition of the AIC statistic for selecting the
#' models. The default value, 0, indicates the scale should be estimated.
#' @param k the multiple of the number of degrees of freedom used for the
#' penalty. Only k = 2 gives the genuine AIC: k = log(n) is
#' sometimes referred to as BIC or SBC.
#' @details This function is used in \code{step}.
#' @export

extractAIC.srm.logit.result <- function(fit, scale, k = 2, ...) {
res <- eval(fit)
aic = -2*res$llf + k * res$df
return(c(res$df, aic))
}

#' Extract the number of observations (S3 method)
#'
#' Extract the number of observations from a model fit. This is principally
#' intended to be used in computing BIC (see AIC).
#'
#' @param object A fitted model object.
#' @param ... Other parameters
#' @details This function is used in \code{step}.
#' @export

nobs.srm.logit.result <- function(object, ...) {
return(nrow(object$data))
}
