## log-gumbel distribution (Frechet type)

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

plgumbel <- function(q, loclog = 0, scalelog = 1, lower.tail = TRUE, log.p = FALSE) {
  pgumbel(log(q), loc=loclog, scale=scalelog, lower.tail=lower.tail, log.p=log.p)
}

qlgumbel <- function(p, loclog = 0, scalelog = 1, lower.tail = TRUE, log.p = FALSE) {
  exp(qgumbel(p, loc=loclog, scale=scalelog, lower.tail=lower.tail, log.p=log.p))
}

## log-gumbel.min distribution (Weibull distribution)

dlgumbel.min <- function(x, loclog = 0, scalelog = 1, log = FALSE) {
  v <- dgumbel.min(log(x), loc=loclog, scale=scalelog, log=FALSE)
  if (log == FALSE) {
    r <- v / x
  } else {
    r <- log(v) - log(x)
  }
  r[v == 0] <- 0
  r
}

plgumbel.min <- function(q, loclog = 0, scalelog = 1, lower.tail = TRUE, log.p = FALSE) {
  pgumbel.min(log(q), loc=loclog, scale=scalelog, lower.tail=lower.tail, log.p=log.p)
}

qlgumbel.min <- function(p, loclog = 0, scalelog = 1, lower.tail = TRUE, log.p = FALSE) {
  exp(qgumbel.min(p, loc=loclog, scale=scalelog, lower.tail=lower.tail, log.p=log.p))
}
