#' Estimate SRM parameters.
#'
#' This function provides the maximum likelihood (ML) estiamtes for a given SRM with a given data.
#' The ML estimates are computed with the EM algorithm. The initial parameters for the EM algorithm
#' are automatically decided if the flag \emph{initialize} is TRUE.
#'
#' @param srm A model.
#' @param data A faultdata.
#' @param initialize Either TRUE or FALSE. If TRUE, the model parameters are
#' initilized with a given data before executing the fitting algorithm.
#' @param maxiter An integer for the maximum number of iterations in the fitting algorithm.
#' @param reltol A numeric value. The algorithm stops if the relative error is
#' less than \emph{reltol} and the absolute error is less than \emph{abstol}.
#' @param abstol A numeric value. The algorithm stops if the relative error is
#' less than \emph{reltol} and the absolute error is less than \emph{abstol}.
#' @param stopcond A character string. \emph{stopcond} gives the criterion
#' for the stop condition of the algorithm. Either llf or parameter is selected.
#' @param trace A logical. If TRUE, the intermediate parameters are printed.
#' @param printsteps An integer for print.
#' @param ... A list for other parameters which are sent to the \code{em} method of srm.
#' @return A list with components;
#' \item{initial}{A vector for initial parameters.}
#' \item{srm}{A class of NHPP. The SRM with the estiamted parameters.}
#' \item{llf}{A numeric value for the maximum log-likelihood function.}
#' \item{df}{An integer for degrees of freedom.}
#' \item{convergence}{A boolean meaning the alorigthm is converged or not.}
#' \item{iter}{An integer for the number of iterations.}
#' \item{aerror}{A numeric value for absolute error.}
#' \item{rerror}{A numeric value for relative error.}
#' @examples
#' data(tomcat5)
#' data <- faultdata(time=time, fault=fault, data=tomcat5.catalina)
#' emfit(srm("exp"), data)
#' @export

emfit <- function(srm, data,
  initialize = TRUE,
  maxiter = 2000, reltol = 1.0e-6, abstol = 1.0e-3,
  stopcond = "llf", trace = FALSE, printsteps = 50, ...) {
    ## init
    if (initialize) {
      srm$init_params(data)
    }

    iter <- 1
    conv <- FALSE
    param <- srm$params

    ## stopcond
    if (stopcond == "parameter" || stopcond == "param") {
      term.fn <- function(res0, res1) {
        sdiff <- res1$llf - res0$llf
        para0 <- c(res0$param)
        para1 <- c(res1$param)
        aerror <- max(abs(para1 - para0))
        rerror <- aerror / max(abs(para0))
        return(c(aerror, rerror, sdiff))
      }
    } else if (stopcond == "llf") {
      term.fn <- function(res0, res1) {
        sdiff <- res1$llf - res0$llf
        aerror <- abs(res1$llf - res0$llf)
        rerror <- aerror / abs(res0$llf)
        return(c(aerror, rerror, sdiff))
      }
    } else {
      stop("wrong stopcond condition.")
    }

    ## repeat emstep
    res0 <- list(param=param, llf=-Inf)

    repeat {
      res1 <- srm$em(res0$param, data, ...)
      error <- term.fn(res0, res1)

      if (trace) {
        if (iter %% printsteps == 0) {
          cat("llf=", res1$llf, "(", error[3], ") params=(", res1$param, ")\n")
        }
      }

      if (!is.finite(res1$llf) || !all(is.finite(res1$param))) {
        warning(sprintf("LLF or param becomes +-Inf, NaN or NA: %s %d",
          srm$name, iter))
        res1 <- res0
        break
      }
      if (error[3] < 0) {
        warning(sprintf("LLF decreses: %s %d %e", srm$name,
          iter, error[3]))
      }
      if ((error[1] < abstol) && (error[2] < reltol)) {
        conv <- TRUE
        srm$set_params(res1$param)
        break
      }
      if (iter >= maxiter) {
        warning("Did not converge to MLE by max iteration.")
        srm$set_params(res1$param)
        break
      }
      iter <- iter + 1
      res0 <- res1
    }

    srm$set_data(data)

    result <- list(
      initial = param,
      srm = srm,
      llf = res1$llf,
      df = srm$df,
      convergence = conv,
      iter = iter,
      aerror = error[1],
      rerror = error[2]
    )
    result
}
