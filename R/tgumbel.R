#' truncated Gumbel distribution (maximum)
#'
#' Density, distribution function, and quantile function for
#' truncated Gumbel distribution (maximum).
#'
#' The truncated Gumbel distribution (maximum) has the cumulative distribution function
#' F(q) = 1-(1-G(q))/(1-G(0)), G(q) = exp(-exp(-z)) where z = (q-loc)/scale
#'
#' @name tgumbel
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param loc A numeric value of location parameter for Gumbel distribution.
#' @param scale A numeric value of scale parameter for Gumbel distribution.
#' @param lower.tail A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].
#' @param log,log.p A logical; if TRUE, the probability p is given as log(p).
#' @return 'dtgumble' gives the desity, 'ptgumbel' gives the distribution, and
#' 'qtgumbel' gives the quantile function.
NULL

#' @rdname tgumbel
dtgumbel <- function(x, loc = 0, scale = 1, log = FALSE) {
  v <- dgumbel(x, loc=loc, scale=scale, log=FALSE) /
    pgumbel(0, loc=loc, scale=scale, lower=FALSE)
  if (log) {
    log(v)
  } else {
    v
  }
}

#' @rdname tgumbel
ptgumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (lower.tail) {
    v <- 1 - pgumbel(q, loc=loc, scale=scale, lower=FALSE, log=FALSE) /
      pgumbel(0, loc=loc, scale=scale, lower=FALSE, log=FALSE)
  } else {
    v <- pgumbel(q, loc=loc, scale=scale, lower=FALSE, log=FALSE) /
      pgumbel(0, loc=loc, scale=scale, lower=FALSE, log=FALSE)
  }
  if (log.p) {
    log(v)
  } else {
    v
  }
}

#' @rdname tgumbel
qtgumbel <- function(p, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  pdash <- (1 - p) * pgumbel(0, loc=loc, scale=scale, lower=FALSE, log=FALSE)
  qgumbel(p=pdash, loc=loc, scale=scale, lower=FALSE, log=FALSE)
}

#' @rdname tgumbel
rtgumbel <- function(n, loc = 0, scale = 1) {
  qtgumbel(runif(n), loc=loc, scale=scale)
}

#' truncated Gumbel distribution (minimum)
#'
#' Density, distribution function, and quantile function for
#' truncated Gumbel distribution (minimum).
#'
#' @name tgumbel.min
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param loc A numeric value of location parameter for Gumbel distribution.
#' @param scale A numeric value of scale parameter for Gumbel distribution.
#' @param lower.tail A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].
#' @param log,log.p A logical; if TRUE, the probability p is given as log(p).
#' @return 'dtgumble.min' gives the desity, 'ptgumbel.min' gives the distribution, and
#' 'qtgumbel.min' gives the quantile function.
#' @details
#' The truncated Gumbel distribution (maximum) has the cumulative distribution function
#' F(q) = 1-(1-G(q))/(1-G(0)), G(q) = exp(-exp(-z)) where z = (q-loc)/scale
NULL
#> NULL

#' @rdname tgumbel.min
dtgumbel.min <- function(x, loc = 0, scale = 1, log = FALSE) {
  v <- dgumbel(x, loc=loc, scale=scale, log=FALSE, min=TRUE) /
    pgumbel(0, loc=loc, scale=scale, lower=FALSE, min=TRUE)
  if (log) {
    log(v)
  } else {
    v
  }
}

#' @rdname tgumbel.min
ptgumbel.min <- function(q, loc = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  if (lower.tail) {
    v <- 1 - pgumbel(q, loc=loc, scale=scale, lower=FALSE, log=FALSE, min=TRUE) /
      pgumbel(0, loc=loc, scale=scale, lower=FALSE, log=FALSE, min=TRUE)
  } else {
    v <- pgumbel(q, loc=loc, scale=scale, lower=FALSE, log=FALSE, min=TRUE) /
      pgumbel(0, loc=loc, scale=scale, lower=FALSE, log=FALSE, min=TRUE)
  }
  if (log.p) {
    log(v)
  } else {
    v
  }
}

#' @rdname tgumbel.min
qtgumbel.min <- function(p, loc = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  pdash <- (1 - p) * pgumbel(0, loc=loc, scale=scale, lower=FALSE, log=FALSE, min=TRUE)
  qgumbel(p=pdash, loc=loc, scale=scale, lower=FALSE, log=FALSE, min=TRUE)
}

#' @rdname tgumbel.min
rtgumbel.min <- function(n, loc = 0, scale = 1) {
  qtgumbel.min(runif(n), loc=loc, scale=scale)
}

