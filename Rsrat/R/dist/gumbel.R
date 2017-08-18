### Gumbel distribution

dgumbel <- function (x, loc = 0, scale = 1, log = FALSE) {
  z <- (x - loc) / scale
  if (log) {
    -log(scale) - z - exp(-z)
  } else {
    exp(-z) * exp(-exp(-z)) / scale
  }
}

pgumbel <- function (q, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  z <- (q - loc) / scale
  p <- exp(-exp(-z))
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

qgumbel <- function (p, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  z <- -log(-log(p))
  scale * z + loc
}

### Gumbel distribution (min)

dgumbel.min <- function (x, loc = 0, scale = 1, log = FALSE) {
  z <- ((-x) - loc) / scale
  if (log) {
    -log(scale) - z - exp(-z)
  } else {
    exp(-z) * exp(-exp(-z)) / scale
  }
}

pgumbel.min <- function (q, loc = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  z <- ((-q) - loc) / scale
  p <- exp(-exp(-z))
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

qgumbel.min <- function (p, loc = 0, scale = 1, lower.tail = TRUE,
  log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == FALSE) {
    p <- 1 - p
  }
  z <- -log(-log(p))
  -(scale * z + loc)
}
