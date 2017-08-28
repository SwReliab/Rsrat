#' log-logistic distribution
#'
#' Density, distribution function, and quantile function for
#' log-logistic distribution.
#'
#' The log-logistic distribution has the cumulative distribution function
#' F(q) = G(log(q)) where G(q) is the cumulative distribution function
#' of the logistic distribution with location and scale.
#'
#' @name llogis
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param locationlog A numeric value of location parameter for the logistic distribution.
#' @param scalelog A numeric value of scale parameter for the logistic distribution.
#' @param lower.tail A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].
#' @param log,log.p A logical; if TRUE, the probability p is given as log(p).
#' @return 'dllogis' gives the desity, 'pllogis' gives the distribution, and
#' 'qllogis' gives the quantile function.
#' @importFrom stats dlogis plogis qlogis
NULL

#' @rdname llogis
dllogis <- function(x, locationlog = 0, scalelog = 1, log = FALSE) {
  v <- dlogis(log(x), location=locationlog, scale=scalelog, log=FALSE)
  if (log == FALSE) {
    r <- v / x
  } else {
    r <- log(v) - log(x)
  }
  r[v == 0] <- 0
  r
}

#' @rdname llogis
pllogis <- function(q, locationlog = 0, scalelog = 1, lower.tail = TRUE,
  log.p = FALSE) {
  plogis(log(q), location=locationlog, scale=scalelog, lower.tail=lower.tail,
    log.p=log.p)
}

#' @rdname llogis
qllogis <- function(p, locationlog = 0, scalelog = 1, lower.tail = TRUE,
  log.p = FALSE) {
  exp(qlogis(p, location=locationlog, scale=scalelog, lower.tail=lower.tail,
    log.p=log.p))
}
