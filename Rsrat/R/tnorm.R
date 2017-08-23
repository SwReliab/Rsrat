#' truncated normal distribution
#'
#' Density, distribution function, and quantile function for
#' truncated normal distribution.
#'
#' @name tnorm
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param mean A numeric value of mean parameter for the normal distirbution.
#' @param sd A numeric value of standard deviation for the normal distirbution.
#' @param lower.tail A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].
#' @param log,log.p A logical; if TRUE, the probability p is given as log(p).
#' @return 'dtnorm' gives the desity, 'ptnorm' gives the distribution, and
#' 'qtnorm' gives the quantile function.
#' @details
#' The truncated normal distribution has the cumulative distribution function
#' F(q) = 1-(1-G(q))/(1-G(0)), where G(q) is the cumulative distribution function of
#' normal distribution with mean and sd.
#' @importFrom stats dnorm pnorm qnorm
NULL
#> NULL

#' @rdname tnorm
dtnorm <- function(x, mean = 0, sd = 1, log = FALSE) {
  v <- dnorm(x, mean=mean, sd=sd, log=FALSE) /
    pnorm(0, mean=mean, sd=sd, lower.tail=FALSE)
  if (log) {
    log(v)
  } else {
    v
  }
}

#' @rdname tnorm
ptnorm <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  if (lower.tail) {
    v <- 1 - pnorm(q, mean=mean, sd=sd, lower.tail=FALSE, log.p=FALSE) /
      pnorm(0, mean=mean, sd=sd, lower.tail=FALSE, log.p=FALSE)
  } else {
    v <- pnorm(q, mean=mean, sd=sd, lower.tail=FALSE, log.p=FALSE) /
      pnorm(0, mean=mean, sd=sd, lower.tail=FALSE, log.p=FALSE)
  }
  if (log.p) {
    log(v)
  } else {
    v
  }
}

#' @rdname tnorm
qtnorm <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  pdash <- (1 - p) * pnorm(0, mean=mean, sd=sd, lower.tail=FALSE, log.p=FALSE)
  qnorm(p=pdash, mean=mean, sd=sd, lower.tail=FALSE, log.p=FALSE)
}
