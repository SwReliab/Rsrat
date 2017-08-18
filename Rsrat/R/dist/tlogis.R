## truncated logistic distribution

dtlogis <- function(x, location = 0, scale = 1, log = FALSE) {
  v <- dlogis(x, location=location, scale=scale, log=FALSE) /
    plogis(0, location=location, scale=scale, lower.tail=FALSE)
  if (log) {
    log(v)
  } else {
    v
  }
}

ptlogis <- function(q, location = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  if (lower.tail) {
    v <- 1 - plogis(q, location=location, scale=scale, lower.tail=FALSE,
      log.p=FALSE) /
      plogis(0, location=location, scale=scale, lower.tail=FALSE, log.p=FALSE)
  } else {
    v <- plogis(q, location=location, scale=scale, lower.tail=FALSE,
      log.p=FALSE) /
      plogis(0, location=location, scale=scale, lower.tail=FALSE, log.p=FALSE)
  }
  if (log.p) {
    log(v)
  } else {
    v
  }
}

qtlogis <- function(p, location = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  pdash <- (1 - p) * plogis(0, location=location, scale=scale,
    lower.tail=FALSE, log.p=FALSE)
  qlogis(p=pdash, location=location, scale=scale, lower.tail=FALSE,
    log.p=FALSE)
}
