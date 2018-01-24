#' Penalized GLM
#'
#' Fit the generalized linear model under a constraint.
#'
#' The control argument is a list that can supply any of the following components:
#' \describe{
#'   \item{maxiter}{An integer for the maximum number of iterations.}
#'   \item{reltol}{A numeric value. The algorithm stops if the relative error is
#' less than \emph{reltol} and the absolute error is less than \emph{abstol}.}
#'   \item{abstol}{A numeric value. The algorithm stops if the relative error is
#' less than \emph{reltol} and the absolute error is less than \emph{abstol}.}
#'   \item{stopcond}{A character string. \emph{stopcond} gives the criterion
#' for the stop condition of the algorithm. Either llf or parameter is selected.}
#'   \item{trace}{A logical. If TRUE, the intermediate parameters are printed.}
#'   \item{printsteps}{An integer for print.}
#' }
#'
#' @param x A model matrix.
#' @param y fault An integer vector for the number of faults detected in time intervals.
#' The fault detected just at the end of time interal is not counted.
#' @param coef A numeric vector indicating starting coefficients.
#' @param offset An integer. This can be used to specify an a priori known
#' component to be included in the linear predictor during fitting. This should
#' be NULL or a numeric vector of length equal to the number of cases.
#' @param family A description of the error distribution and link function to
#' be used in the model.
#' @param lambda A numeric value indicating the penalized parameter. When lambda = 0, coefficients are unconstraint.
#' @param K A matrix to detect the contraint structure. The default is an identity matrix.
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
#' \item{ctime}{A numeric value for computation time.}
#' \item{call}{The method call.}
#' @export

pnglm.fit <- function(x, y, coef = NULL, offset = NULL, family = gaussian(),
    lambda = 0, K = NULL, control = list(), ...) {

    ## TODO: we should make the routine for removing elements with both 0 (see "good" in glm.fit)

    con <- pnglm.options()
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))

    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) {
      rownames(y)
    }
    else names(y)

    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    if (lambda == 0) {
        nobs2 <- nobs
        penalized <- FALSE
    }
    else {
        nobs2 <- nobs + nvars
        penalized <- TRUE
    }

    if (is.null(offset)) {
        offset <- rep(0, nobs)
    }
    weights <- rep(1, nobs)
    if (is.null(coef)) {
        mustart <- NULL
        eval(family$initialize)
        eta <- family$linkfun(mustart)
    }
    else {
        eval(family$initialize)
        eta <- x %*% coef + offset
        if (!family$valideta(eta)) {
            stop(gettextf("eta is invalid"))
        }
    }

    if (penalized && is.null(K)) {
    	K <- diag(rep(1, nvars))
    }

    if (penalized) {
      xdash <- rbind(x, K)
    }
    else {
      xdash <- x
    }

    mu <- family$linkinv(eta)
    if (!family$validmu(mu)) {
        stop(gettextf("mu is invalid"))
    }
    var <- family$variance(mu)
    mu.eta <- family$mu.eta(eta)
    if (penalized) {
      w <- diag(c(c(sqrt(weights * var / (weights * mu.eta)^2)), rep(1.0/(2.0 * lambda), nvars)))
      z <- c(eta - offset + (y - mu) / mu.eta, rep(0, nvars))
    }
    else {
      w <- diag(c(sqrt(weights * var / (weights * mu.eta)^2)))
      z <- eta - offset + (y - mu) / mu.eta
    }
    devold <- sum(family$dev.resids(y, mu, weights))

    res <- dggglm(xdash, w, z)

    coef <- res$x
    coefold <- coef

    for (iter in 1:con$maxiter) {
        eta <- x %*% coef + offset
        if (!family$valideta(eta)) {
            stop(gettextf("eta is invalid at %d", iter))
        }
        mu <- family$linkinv(eta)
        if (!family$validmu(mu)) {
            stop(gettextf("mu is invalid at %d", iter))
        }
        var <- family$variance(mu)
        mu.eta <- family$mu.eta(eta)
        if (penalized) {
          w <- diag(c(c(sqrt(weights * var / (weights * mu.eta)^2)), rep(1.0/(2.0 * lambda), nvars)))
          z <- c(eta - offset + (y - mu) / mu.eta, rep(0, nvars))
        }
        else {
          w <- diag(c(sqrt(weights * var / (weights * mu.eta)^2)))
          z <- eta - offset + (y - mu) / mu.eta
        }
        dev <- sum(family$dev.resids(y, mu, weights))

        res <- dggglm(xdash, w, z)
        coef <- res$x

        if (con$trace) {
          if (iter %% con$printsteps == 0) {
            message(gettextf("pnglm: iter=%d dev=%e", iter, dev))
          }
        }

        if (abs(dev - devold)/(0.1 + abs(dev)) < con$reltol) {
            conv <- TRUE
            break
        }

        devold <- dev
        coefold <- coef
    }

    names(coef) <- xnames

    eta <- x %*% coef + offset
    if (!family$valideta(eta)) {
        stop(gettextf("eta is invalid at %d", iter))
    }
    mu <- as.vector(family$linkinv(eta))
    if (!family$validmu(mu)) {
        stop(gettextf("mu is invalid at %d", iter))
    }
    names(mu) <- ynames

    if (penalized) {
      names(res$y) <- c(ynames, xnames)
      list(coefficients=coef, residual=res$y, fitted.values=mu, K=K, lambda=lambda)
    }
    else {
      names(res$y) <- ynames
      list(coefficients=coef, residual=res$y, fitted.values=mu)
    }
}

#' Options for pnglm.fit
#'
#' Generate a list of option values.
#'
#' @return A list of options.
#' @rdname pnglm.fit
#' @export

pnglm.options <- function() {
  list(maxiter = 2000,
    reltol = sqrt(.Machine$double.eps),
    abstol = 1.0e+200,
    stopcond = NULL,
    trace = FALSE,
    printsteps = 50)
}
