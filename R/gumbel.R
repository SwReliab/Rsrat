#' Gumbel distribution
#'
#' Density, distribution function, and quantile function for
#' Gumbel distribution.
#'
#' The Gumbel distribution (maximum) has the cumulative distribution function
#' F(q) = exp(-exp(-z)) where z = (q-loc)/scale
#' The Gumbel distribution (minimum) has the cumulative distribution function
#' F(q) = 1-exp(-exp(-z)) where z = (-q-loc)/scale
#'
#' @name gumbel
#' @param x,q A numeric vector of quantiles.
#' @param p A numeric vector of probabilities.
#' @param loc A numeric value of location parameter.
#' @param scale A numeric value of scale parameter.
#' @param lower A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].
#' @param log A logical; if TRUE, the probability p is given as log(p).
#' @param min A logical; if TRUE, the result of minimum Gumbel is provided.
#' @return 'dgumble' gives the desity, 'pgumbel' gives the distribution, and
#' 'qgumbel' gives the quantile function.
NULL
