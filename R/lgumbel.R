#' log-Gumbel distribution (Frechet)
#'
#' Density, distribution function, and quantile function for
#' log-Gumbel (Frechet) distribution.
#'
#' The log-Gumbel distribution (maximum) has the cumulative distribution function
#' F(q) = exp(-exp(-z)) where z = (log(q)-loc)/scale
#'
#' @name lgumbel
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param loclog A numeric value of location parameter for Gumbel distribution.
#' @param scalelog A numeric value of scale parameter for Gumbel distribution.
#' @param lower.tail A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].
#' @param log,log.p A logical; if TRUE, the probability p is given as log(p).
#' @return 'dlgumble' gives the desity, 'plgumbel' gives the distribution, and
#' 'qlgumbel' gives the quantile function.
NULL

#' @rdname lgumbel
dlgumbel <- function(x, loclog = 0, scalelog = 1, log = FALSE) {
  v <- dgumbel(log(x), loc=loclog, scale=scalelog, log=FALSE)
  if (log == FALSE) {
    r <- v / x
  } else {
    r <- log(v) - log(x)
  }
  r[v == 0] <- 0
  r
}

#' @rdname lgumbel
plgumbel <- function(q, loclog = 0, scalelog = 1, lower.tail = TRUE,
  log.p = FALSE) {
  pgumbel(log(q), loc=loclog, scale=scalelog, lower=lower.tail, log=log.p)
}

#' @rdname lgumbel
qlgumbel <- function(p, loclog = 0, scalelog = 1, lower.tail = TRUE,
  log.p = FALSE) {
  exp(qgumbel(p, loc=loclog, scale=scalelog, lower=lower.tail, log=log.p))
}

#' @rdname lgumbel
rlgumbel <- function(n, loclog = 0, scalelog = 1) {
  exp(rgumbel(n, loc=loclog, scale=scalelog))
}

#' log-Gumbel distribution (Weibull)
#'
#' Density, distribution function, and quantile function for
#' log-Gumbel (Weibull) distribution.
#'
#' The log-Gumbel distribution (minimum) has the cumulative distribution function
#' F(q) = 1-exp(-exp(-z)) where z = (-log(q)-loc)/scale
#'
#' @name lgumbel.min
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param loclog A numeric value of location parameter for Gumbel distribution.
#' @param scalelog A numeric value of scale parameter for Gumbel distribution.
#' @param lower.tail A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].
#' @param log,log.p A logical; if TRUE, the probability p is given as log(p).
#' @return 'dlgumble.min' gives the desity, 'plgumbel.min' gives the distribution, and
#' 'qlgumbel.min' gives the quantile function.
NULL

#' @rdname lgumbel.min
dlgumbel.min <- function(x, loclog = 0, scalelog = 1, log = FALSE) {
  v <- dgumbel(log(x), loc=loclog, scale=scalelog, log=FALSE, min=TRUE)
  if (log == FALSE) {
    r <- v / x
  } else {
    r <- log(v) - log(x)
  }
  r[v == 0] <- 0
  r
}

#' @rdname lgumbel.min
plgumbel.min <- function(q, loclog = 0, scalelog = 1, lower.tail = TRUE,
  log.p = FALSE) {
  pgumbel(log(q), loc=loclog, scale=scalelog, lower=lower.tail, log=log.p, min=TRUE)
}

#' @rdname lgumbel.min
qlgumbel.min <- function(p, loclog = 0, scalelog = 1, lower.tail = TRUE,
  log.p = FALSE) {
  exp(qgumbel(p, loc=loclog, scale=scalelog, lower=lower.tail, log=log.p, min=TRUE))
}

#' @rdname lgumbel.min
rlgumbel.min <- function(n, loclog = 0, scalelog = 1) {
  exp(rgumbel(n, loc=loclog, scale=scalelog, min=TRUE))
}
