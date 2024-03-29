% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lgumbel.R
\name{lgumbel}
\alias{lgumbel}
\alias{dlgumbel}
\alias{plgumbel}
\alias{qlgumbel}
\alias{rlgumbel}
\title{log-Gumbel distribution (Frechet)}
\usage{
dlgumbel(x, loclog = 0, scalelog = 1, log = FALSE)

plgumbel(q, loclog = 0, scalelog = 1, lower.tail = TRUE, log.p = FALSE)

qlgumbel(p, loclog = 0, scalelog = 1, lower.tail = TRUE, log.p = FALSE)

rlgumbel(n, loclog = 0, scalelog = 1)
}
\arguments{
\item{x, q}{A numeric vector of quantiles.}

\item{loclog}{A numeric value of location parameter for Gumbel distribution.}

\item{scalelog}{A numeric value of scale parameter for Gumbel distribution.}

\item{log, log.p}{A logical; if TRUE, the probability p is given as log(p).}

\item{lower.tail}{A logical; if TRUE, the probability is P[X <= x], otherwise, P[X > x].}

\item{p}{A numeric vector of probabilities.}
}
\value{
'dlgumble' gives the desity, 'plgumbel' gives the distribution, and
'qlgumbel' gives the quantile function.
}
\description{
Density, distribution function, and quantile function for
log-Gumbel (Frechet) distribution.
}
\details{
The log-Gumbel distribution (maximum) has the cumulative distribution function
F(q) = exp(-exp(-z)) where z = (log(q)-loc)/scale
}
