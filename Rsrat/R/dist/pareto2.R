### Pareto distribution type 2 (Lomax distribution)

dpareto2 <- function(x, shape = 1, scale = 1, log = FALSE) {
  if (log) {
    log(shape) + shape * log(scale) - (shape+1) * log(x + scale)
  } else {
    (shape / scale) * (scale / (x + scale))^(shape + 1)
  }
}

ppareto2 <- function (q, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  p <- (scale / (q + scale))^shape
  if (lower.tail == TRUE) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}

qpareto2 <- function (p, shape = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) {
    p <- exp(p)
  }
  if (lower.tail == TRUE) {
    p <- 1 - p
  }
  scale * (p^(-1/shape) - 1)
}
