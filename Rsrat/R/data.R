
faultdata.time <- function(time, te) {
  if (missing(te)) {
    len <- length(time)
    fault <- rep(0L, len)
    type <- rep(1L, len)
  } else {
    time <- c(time, te)
    len <- length(time)
    fault <- rep(0L, len)
    type <- rep(1L, len)
    type[len] <- 0L
  }
  faultdata(time, fault, type)
}

faultdata <- function(time, fault = rep.int(0L, length(time)),
  type = rep.int(0L, length(time))) {

  if (length(time) != length(fault))
    stop("Invalid data")
  if (length(time) != length(type))
      stop("Invalid data")
  if (any(time == 0 & fault != 0L & type != 0L))
    stop("Invalid data: zero time exits.")

  tmp <- fault + type
  total <- sum(tmp)
  tmean <- sum(cumsum(time) * tmp) / total
  tmax <- max(cumsum(time)[tmp >= 1L])
  result <- list(
    time=time,
    fault=fault,
    type=type,
    total=total,
    len=length(time),
    mean=tmean,
    max=tmax)
  class(result) <- "Rsrat.faultdata"
  result
}

print.Rsrat.faultdata <- function (x, ..., digits = NULL, quote = FALSE,
  right = TRUE, row.names = TRUE) {
    df <- data.frame(time=x$time, fault=x$fault, type=x$type)
    print.data.frame(df, ..., digits=digits, quote=quote,
      right=right, row.names=row.names)
    invisible(x)
}
