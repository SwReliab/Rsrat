## truncated gumbel distribution

dtgumbel <- function(x, loc = 0, scale = 1, log = FALSE) {
  v <- dgumbel(x, loc=loc, scale=scale, log=FALSE) /
    pgumbel(0, loc=loc, scale=scale, lower.tail=FALSE)
  if (log) {
    log(v)
  } else {
    v
  }
}

ptgumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (lower.tail) {
    v <- 1 - pgumbel(q, loc=loc, scale=scale, lower.tail=FALSE, log.p=FALSE) /
      pgumbel(0, loc=loc, scale=scale, lower.tail=FALSE, log.p=FALSE)
  } else {
    v <- pgumbel(q, loc=loc, scale=scale, lower.tail=FALSE, log.p=FALSE) /
      pgumbel(0, loc=loc, scale=scale, lower.tail=FALSE, log.p=FALSE)
  }
  if (log.p) {
    log(v)
  } else {
    v
  }
}

qtgumbel <- function(p, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  pdash <- (1 - p) * pgumbel(0, loc=loc, scale=scale, lower.tail=FALSE,
    log.p=FALSE)
  qgumbel(p=pdash, loc=loc, scale=scale, lower.tail=FALSE, log.p=FALSE)
}

## truncated gumbel distribution (min)

dtgumbel.min <- function(x, loc = 0, scale = 1, log = FALSE) {
  v <- dgumbel.min(x, loc=loc, scale=scale, log=FALSE) /
    pgumbel.min(0, loc=loc, scale=scale, lower.tail=FALSE)
  if (log) {
    log(v)
  } else {
    v
  }
}

ptgumbel.min <- function(q, loc = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  if (lower.tail) {
    v <- 1 - pgumbel.min(q, loc=loc, scale=scale, lower.tail=FALSE,
      log.p=FALSE) /
      pgumbel.min(0, loc=loc, scale=scale, lower.tail=FALSE, log.p=FALSE)
  } else {
    v <- pgumbel.min(q, loc=loc, scale=scale, lower.tail=FALSE, log.p=FALSE) /
      pgumbel.min(0, loc=loc, scale=scale, lower.tail=FALSE, log.p=FALSE)
  }
  if (log.p) {
    log(v)
  } else {
    v
  }
}

qtgumbel.min <- function(p, loc = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  pdash <- (1 - p) * pgumbel.min(0, loc=loc, scale=scale, lower.tail=FALSE,
    log.p=FALSE)
  qgumbel.min(p=pdash, loc=loc, scale=scale, lower.tail=FALSE, log.p=FALSE)
}
