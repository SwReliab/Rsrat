#' Gumbel distribution (maximum)
#'
#' Density, distribution function, and quantile function for
#' Gumbel distribution (maximum).
#'
#' @name gumbel
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param loc A numeric value of location parameter.
#' @param scale A numeric value of scale parameter.
#' @param lower.tail A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].
#' @param log,log.p A logical; if TRUE, the probability p is given as log(p).
#' @return 'dgumble' gives the desity, 'pgumbel' gives the distribution, and
#' 'qgumbel' gives the quantile function.
#' @details
#' The Gumbel distribution (maximum) has the cumulative distribution function
#' F(q) = exp(-exp(-z)) where z = (q-loc)/scale
NULL
#> NULL

#' @rdname gumbel
dgumbel <- function (x, loc = 0, scale = 1, log = FALSE) {
  z <- (x - loc) / scale
  if (log) {
    -log(scale) - z - exp(-z)
  } else {
    exp(-z) * exp(-exp(-z)) / scale
  }
}

#' @rdname gumbel
pgumbel <- function (q, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  z <- (q - loc) / scale
  p <- exp(-exp(-z))
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

#' @rdname gumbel
qgumbel <- function (p, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  z <- -log(-log(p))
  scale * z + loc
}

#' Gumbel distribution (minimum)
#'
#' Density, distribution function, and quantile function for
#' Gumbel distribution (minimum).
#'
#' @name gumbel.min
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param loc A numeric value of location parameter.
#' @param scale A numeric value of scale parameter.
#' @param lower.tail A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].
#' @param log,log.p A logical; if TRUE, the probability p is given as log(p).
#' @return 'dgumble.min' gives the desity, 'pgumbel.min' gives the distribution, and
#' 'qgumbel.min' gives the quantile function.
#' @details
#' The Gumbel distribution (minimum) has the cumulative distribution function
#' F(q) = 1-exp(-exp(-z)) where z = (-q-loc)/scale
NULL
#> NULL

#' @rdname gumbel.min
dgumbel.min <- function (x, loc = 0, scale = 1, log = FALSE) {
  z <- ((-x) - loc) / scale
  if (log) {
    -log(scale) - z - exp(-z)
  } else {
    exp(-z) * exp(-exp(-z)) / scale
  }
}

#' @rdname gumbel.min
pgumbel.min <- function (q, loc = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  z <- ((-q) - loc) / scale
  p <- exp(-exp(-z))
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

#' @rdname gumbel.min
qgumbel.min <- function (p, loc = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  z <- -log(-log(p))
  -(scale * z + loc)
}
