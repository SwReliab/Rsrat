#' Software fault data
#'
#' Function faultdata() creates a list to store the fault data that are used to
#' esiamte model parameters of SRM. This data can represent both fault time and
#' count data.
#'
#' @param time A numeric vector for time intervals.
#' @param fault An integer vector for the number of faults detected in time intervals.
#' The fault detected just at the end of time interal is not counted.
#' @param type Either 0 or 1. If 1, a fault is detected just at the end of corresponding time interval.
#' This is used to represent the fault time data. If 0, no fault is detected at the end of interval.
#' @return A list with the attribute class='Rsrat.faultdata';
#' \item{time}{A vector for time interval.}
#' \item{fault}{A vector for the number of detected faults in intervals.}
#' \item{type}{A vector for the indicator whether a fault is detected at the end of intervals.}
#' \item{total}{An integer for the number of total faults.}
#' \item{len}{An integer for the number of time intervals.}
#' \item{mean}{A numeric value for mean fault detection time from data.}
#' \item{max}{A numeric value for maximum of fault detection time.}
#' @export

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

#' @describeIn faultdata Create fault data from fault time data.
#' @param te A numeric value for the time interval from the last fault to the observation time.
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

print.Rsrat.faultdata <- function (x, ..., digits = NULL, quote = FALSE,
  right = TRUE, row.names = TRUE) {
    df <- data.frame(time=x$time, fault=x$fault, type=x$type)
    print.data.frame(df, ..., digits=digits, quote=quote,
      right=right, row.names=row.names)
    invisible(x)
}
