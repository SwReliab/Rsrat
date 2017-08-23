#' truncated logistic distribution
#'
#' Density, distribution function, and quantile function for
#' truncated logistic distribution.
#'
#' @name tlogis
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param location A numeric value of location parameter for the logistic distribution.
#' @param scale A numeric value of scale parameter for the logistic distribution.
#' @param lower.tail A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].
#' @param log,log.p A logical; if TRUE, the probability p is given as log(p).
#' @return 'dtlogis' gives the desity, 'ptlogis' gives the distribution, and
#' 'qtlogis' gives the quantile function.
#' @details
#' The truncated logistic distribution has the cumulative distribution function
#' F(q) = 1-(1-G(q))/(1-G(0)) where G(q) is the cumulative distribution function
#' of the logistic distribution with location and scale.
#' @importFrom stats dlogis plogis qlogis
NULL
#> NULL

#' @rdname tlogis
dtlogis <- function(x, location = 0, scale = 1, log = FALSE) {
  v <- dlogis(x, location=location, scale=scale, log=FALSE) /
    plogis(0, location=location, scale=scale, lower.tail=FALSE)
  if (log) {
    log(v)
  } else {
    v
  }
}

#' @rdname tlogis
ptlogis <- function(q, location = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  if (lower.tail) {
    v <- 1 - plogis(q, location=location, scale=scale, lower.tail=FALSE,
      log.p=FALSE) /
      plogis(0, location=location, scale=scale, lower.tail=FALSE, log.p=FALSE)
  } else {
    v <- plogis(q, location=location, scale=scale, lower.tail=FALSE,
      log.p=FALSE) /
      plogis(0, location=location, scale=scale, lower.tail=FALSE, log.p=FALSE)
  }
  if (log.p) {
    log(v)
  } else {
    v
  }
}

#' @rdname tlogis
qtlogis <- function(p, location = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  pdash <- (1 - p) * plogis(0, location=location, scale=scale,
    lower.tail=FALSE, log.p=FALSE)
  qlogis(p=pdash, location=location, scale=scale, lower.tail=FALSE,
    log.p=FALSE)
}
