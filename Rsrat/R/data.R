
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
  faultdata.gen(time, fault, type)
}

faultdata.group <- function(time, fault) {
  if (length(time) != length(fault))
    stop("Invalid data")

  type <- rep(0L, length(time))
  faultdata.gen(time, fault, type)
}

faultdata.gen <- function(time, fault, type) {
  if (any(time == 0 & fault != 0L & type != 0L))
    stop("Invalid data: zero time exits.")

  total <- sum(fault)
  tmean <- sum(cumsum(time) * (fault + type)) / total
  tmax <- max(cumsum(time)[(fault + type) >= 1L])
  result <- list(
    fault=fault,
    time=time,
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
