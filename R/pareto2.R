#' Pareto distribution (type 2)
#'
#' Density, distribution function, and quantile function for
#' Pareto (type 2) distribution.
#'
#' The Pareto type 2distribution has the cumulative distribution function
#' F(q) = 1 - (scale / (q + scale))^shape
#'
#' @name pareto2
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param shape A numeric value of shape parameter.
#' @param scale A numeric value of scale parameter.
#' @param lower.tail A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].
#' @param log,log.p A logical; if TRUE, the probability p is given as log(p).
#' @return 'dpareto2' gives the desity, 'ppareto2' gives the distribution, and
#' 'qpareto2' gives the quantile function.
NULL

#' @rdname pareto2
dpareto2 <- function(x, shape = 1, scale = 1, log = FALSE) {
  if (log) {
    log(shape) + shape * log(scale) - (shape+1) * log(x + scale)
  } else {
    (shape / scale) * (scale / (x + scale))^(shape + 1)
  }
}

#' @rdname pareto2
ppareto2 <- function (q, shape = 1, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  p <- (scale / (q + scale))^shape
  if (lower.tail == TRUE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

#' @rdname pareto2
qpareto2 <- function (p, shape = 1, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == TRUE) {
    p <- 1 - p
  }
  scale * (p^(-1/shape) - 1)
}
