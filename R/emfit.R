#' Estimate SRM parameters.
#'
#' This function provides the maximum likelihood (ML) estiamtes for a given SRM with a given data.
#' The ML estimates are computed with the EM algorithm. The initial parameters for the EM algorithm
#' are automatically decided if the flag \emph{initialize} is TRUE.
#'
#' @param srm An NHPP model.
#' @param data A faultdata.
#' @param initialize Either TRUE or FALSE. If TRUE, the model parameters are
#' initilized with a given data before executing the fitting algorithm.
#' @param maxiter An integer for the maximum number of iterations in the fitting algorithm.
#' @param rtol A numeric value. The algorithm stops if the relative error is
#' less than \emph{rtol} and the absolute error is less than \emph{atol}.
#' @param atol A numeric value. The algorithm stops if the relative error is
#' less than \emph{rtol} and the absolute error is less than \emph{atol}.
#' @param termination A character string. \emph{termination} gives the criterion
#' for the stop condition of the algorithm. Either llf or parameter is selected.
#' @return A list with components;
#' \item{initial}{A vector for initial parameters.}
#' \item{srm}{A class of NHPP. The SRM with the estiamted parameters.}
#' \item{llf}{A numeric value for the maximum log-likelihood function.}
#' \item{convergence}{A boolean meaning the alorigthm is converged or not.}
#' \item{iter}{An integer for the number of iterations.}
#' \item{aerror}{A numeric value for absolute error.}
#' \item{rerror}{A numeric value for relative error.}
#' @examples
#' data(musa)
#' emfit(srm("exp"), musa.sys1.group)
#' @export

emfit <- function(srm, data,
  initialize = TRUE,
  maxiter = 2000, rtol = 1.0e-6, atol = 1.0e-3,
  termination = "llf") {
    ## init
    if (initialize) {
      srm$init_params(data)
    }

    iter <- 1
    conv <- FALSE
    param <- srm$params

    ## termination
    if (termination == "parameter" || termination == "param") {
      term.fn <- function(res0, res1) {
        sdiff <- res1$llf - res0$llf
        para0 <- c(res0$param)
        para1 <- c(res1$param)
        aerror <- max(abs(para1 - para0))
        rerror <- aerror / max(abs(para0))
        return(c(aerror, rerror, sdiff))
      }
    } else if (termination == "llf") {
      term.fn <- function(res0, res1) {
        sdiff <- res1$llf - res0$llf
        aerror <- abs(res1$llf - res0$llf)
        rerror <- aerror / abs(res0$llf)
        return(c(aerror, rerror, sdiff))
      }
    } else {
      stop("wrong termination condition.")
    }

    ## repeat emstep
    res0 <- list(param=param, llf=-Inf)

    repeat {
      res1 <- srm$em(res0$param, data)
      error <- term.fn(res0, res1)

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
      if ((error[1] < atol) && (error[2] < rtol)) {
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


    result <- list(
      initial = param,
      srm = srm,
      llf = res1$llf,
      convergence = conv,
      iter = iter,
      aerror = error[1],
      rerror = error[2]
    )
    class(result) <- "Rsrat.result"
    result
}
