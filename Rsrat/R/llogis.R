## log-logistic distribution

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

pllogis <- function(q, locationlog = 0, scalelog = 1, lower.tail = TRUE,
  log.p = FALSE) {
  plogis(log(q), location=locationlog, scale=scalelog, lower.tail=lower.tail,
    log.p=log.p)
}

qllogis <- function(p, locationlog = 0, scalelog = 1, lower.tail = TRUE,
  log.p = FALSE) {
  exp(qlogis(p, location=locationlog, scale=scalelog, lower.tail=lower.tail,
    log.p=log.p))
}
