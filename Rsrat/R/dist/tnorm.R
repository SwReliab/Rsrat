## truncated normal distribution

dtnorm <- function(x, mean = 0, sd = 1, log = FALSE) {
  v <- dnorm(x, mean=mean, sd=sd, log=FALSE) / pnorm(0, mean=mean, sd=sd, lower.tail=FALSE)
  if (log) {
    log(v)
  } else {
    v
  }
}

ptnorm <- function(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  if (lower.tail) {
    v <- 1 - pnorm(q, mean=mean, sd=sd, lower.tail=FALSE, log.p=FALSE) / pnorm(0, mean=mean, sd=sd, lower.tail=FALSE, log.p=FALSE)
  } else {
    v <- pnorm(q, mean=mean, sd=sd, lower.tail=FALSE, log.p=FALSE) / pnorm(0, mean=mean, sd=sd, lower.tail=FALSE, log.p=FALSE)
  }
  if (log.p) {
    log(v)
  } else {
    v
  }
}

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
