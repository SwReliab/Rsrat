#' Model name
#'
#' This function provides the model name.
#'
#' @param srm An object of a model
#' @return A string of model name

srmname <- function(srm) {
  if ("srm.nhpp.result" %in% class(srm)) {
    srm <- srm$srm
  }
  srm$name
}

#' Mean value function
#'
#' This function provides the mean value function for a given model.
#'
#' @param t A numeric vector for time
#' @param srm An object of a model
#' @return A vector of mean value function
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(fault=tohma, srm.names = c("exp", "gamma"),
#'             selection = NULL)
#' mvf(c(1,2,3), result$exp)
#' @export

mvf <- function(t, srm) {
  if ("srm.nhpp.result" %in% class(srm)) {
    srm <- srm$srm
  }
  srm$mvf(t)
}

#' Difference of mean value function
#'
#' This function provides the difference of mean value function for a given model.
#'
#' @param t A numeric vector for time
#' @param srm An object of a model
#' @return A vector of mean value function
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(fault=tohma, srm.names = c("exp", "gamma"),
#'             selection = NULL)
#' dmvf(c(1,2,3), result$exp)
#' @export

dmvf <- function(t, srm) {
  if ("srm.nhpp.result" %in% class(srm)) {
    srm <- srm$srm
  }
  srm$dmvf(t)
}

#' Rate function
#'
#' This function provides the rate function for a given model.
#'
#' @param t A numeric vector for time
#' @param srm An object of a model
#' @return A vector of rate function
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(fault=tohma, srm.names = c("exp", "gamma"),
#'             selection = NULL)
#' rate(c(1,2,3), result$exp)
#' @export

rate <- function(t, srm) {
  if ("srm.nhpp.result" %in% class(srm)) {
    srm <- srm$srm
  }
  srm$intensity(t)
}

#' Draw a graph of mean value functions
#'
#' This function draws a graph of mean value function with data and NHPP classes.
#'
#' The return value is an object of ggplot2, and thus a line is added like
#' \code{mvfplot(...) + stat_function(fun=mvf)}.
#'
#' @param time A numeric vector for time intervals.
#' @param fault An integer vector for the number of faults detected in time intervals.
#' The fault detected just at the end of time interal is not counted.
#' @param type Either 0 or 1. If 1, a fault is detected just at the end of corresponding time interval.
#' This is used to represent the fault time data. If 0, no fault is detected at the end of interval.
#' @param te A numeric value for the time interval from the last fault to the observation time.
#' @param data A dataframe. The arguments; time, fault, type, te can also be selected as the columns of dataframe.
#' @param srms A list of estimated results.
#' @param xlab A character string, indicating the label of x-axis.
#' @param ylab A character string, indicating the label of y-axis.
#' @param datalab A character string, indicating the label of data in legend.
#' @param xmax A numeric value, indicating the maximum value of xlim. If NA, the maximum value is automatically determined.
#' @param ymax A numeric value, indicating the maximum value of ylim. If NA, the maximum value is automatically determined.
#' @param colors A vector of strings, indicating the colors in plot.
#' @param ... Other parameters.
#' @return An object of ggplot.
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(fault=tohma, srm.names = c("exp", "gamma"),
#'             selection = NULL)
#' mvfplot(fault=tohma, srms=result)
#' @export

mvfplot <- function(time = NULL, fault = NULL, type = NULL, te = NULL, data = data.frame(),
  srms = list(), xlab = "time", ylab = "# of faults", datalab = "data",
  xmax = NA, ymax = NA, colors = mmcolors, ...) {
  if (class(data) != "Rsrat.faultdata") {
    data <- faultdata.nhpp(substitute(time), substitute(fault),
      substitute(type), substitute(te), data, parent.frame())
  }
  present <- sum(data$time)
  if (is.na(xmax)) {
    xmax <- present * 1.2
  }
  n <- data$fault + data$type
  data <- data.frame(x=cumsum(data$time)[n != 0], y=cumsum(n)[n != 0])
  gp <- ggplot(data, aes_string(x="x", y="y")) + labs(x=xlab, y=ylab) + xlim(c(0,xmax)) + ylim(c(0,ymax))
  gp <- gp + geom_point(aes_(colour=datalab))
  gp <- gp + geom_vline(xintercept=present, linetype="dotted")

  if ("list" %in% class(srms)) {
    for (s in srms) {
      gp <- gp + stat_function(fun=mvf, args=list(srm=s), aes_(colour=srmname(s)))
    }
    gp + scale_colour_manual(NULL, values = colors)
  } else {
    gp <- gp + stat_function(fun=mvf, args=list(srm=srms), aes_(colour=srmname(srms)))
    gp + scale_colour_manual(NULL, values = colors)
  }
}

#' Draw a graph of diff of mean value functions on discrete time domain
#'
#' This function draws a graph of mean value function with data and NHPP classes.
#'
#' The return value is an object of ggplot2, and thus a line is added like
#' \code{dmvfplot(...) + geom_point(stat="identity", position="identity", aes_string(x="x", y="exp", colour=shQuote("exp")))
#'  + geom_line(stat="identity", position="identity", aes_string(x="x", y="exp", colour=shQuote("exp")))}.
#'
#' @param time A numeric vector for time intervals.
#' @param fault An integer vector for the number of faults detected in time intervals.
#' The fault detected just at the end of time interval is not counted.
#' @param type Either 0 or 1. If 1, a fault is detected just at the end of corresponding time interval.
#' This is used to represent the fault time data. If 0, no fault is detected at the end of interval.
#' @param te A numeric value for the time interval from the last fault to the observation time.
#' @param data A dataframe. The arguments; time, fault, type, te can also be selected as the columns of dataframe.
#' @param srms A list of estimated results.
#' @param xlab A character string, indicating the label of x-axis.
#' @param ylab A character string, indicating the label of y-axis.
#' @param datalab A character string, indicating the label of data in legend.
#' @param xmax A numeric value, indicating the maximum value of xlim. If NA, the maximum value is automatically determined.
#' @param ymax A numeric value, indicating the maximum value of ylim. If NA, the maximum value is automatically determined.
#' @param colors A vector of strings, indicating the colors in plot.
#' @param ... Other parameters.
#' @return An object of ggplot.
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(fault=tohma, srm.names = c("exp", "gamma"),
#'             selection = NULL)
#' dmvfplot(fault=tohma, srms=result)
#' @export

dmvfplot <- function(time = NULL, fault = NULL, type = NULL, te = NULL, data = data.frame(),
  srms = list(), xlab = "time", ylab = "# of faults", datalab = "data",
  xmax = NA, ymax = NA, colors = mmcolors, ...) {
  if (class(data) != "Rsrat.faultdata") {
    data <- faultdata.nhpp(substitute(time), substitute(fault),
      substitute(type), substitute(te), data, parent.frame())
  }
  present <- sum(data$time)
  t <- as.numeric(cumsum(data$time))
  n <- as.numeric(data$fault + data$type)

  data <- data.frame(x=t, y=n)
  if ("list" %in% class(srms)) {
    for (s in srms) {
      data <- cbind(data, dmvf(t, s))
    }
    colnames(data) <- c("x", "y", sapply(srms, srmname))
  } else {
    data <- cbind(data, dmvf(t, srms))
    colnames(data) <- c("x", "y", srmname(srms))
  }

  gp <- ggplot(data, aes_string(x="x", y="y")) + labs(x=xlab, y=ylab) + xlim(c(0,xmax)) + ylim(c(0,ymax))
  gp <- gp + geom_bar(stat="identity", position="identity", alpha=0.5, aes_(fill=datalab))
  # gp <- gp + geom_vline(xintercept=present, linetype="dotted")
  if ("list" %in% class(srms)) {
    for (s in srms) {
      # gp <- gp + geom_point(stat="identity", position="identity", aes_string(x="x", y=srmname(s), colour=shQuote(srmname(s))))
      gp <- gp + geom_line(stat="identity", position="identity", aes_string(x="x", y=srmname(s), colour=shQuote(srmname(s))))
    }
  } else {
    # gp <- gp + geom_point(stat="identity", position="identity", aes_string(x="x", y=srmname(srms), colour=shQuote(srmname(srms))))
    gp <- gp + geom_line(stat="identity", position="identity", aes_string(x="x", y=srmname(srms), colour=shQuote(srmname(srms))))
  }
  gp + scale_colour_manual(NULL, values = colors[2:length(colors)]) + scale_fill_manual(NULL, values = colors[1])
}

#' Draw a graph of rate functions on continuous time domain
#'
#' This function draws a graph of rate functions with data and NHPP classes.
#'
#' The return value is an object of ggplot2, and thus a line is added like
#' \code{rateplot(...) + geom_point(stat="identity", position="identity", aes_string(x="x", y="exp", colour=shQuote("exp")))
#'  + geom_line(stat="identity", position="identity", aes_string(x="x", y="exp", colour=shQuote("exp")))}.
#'
#' @param time A numeric vector for time intervals.
#' @param fault An integer vector for the number of faults detected in time
#' intervals. The fault detected just at the end of time interal is not counted.
#' @param type Either 0 or 1. If 1, a fault is detected just at the end of
#' corresponding time interval. This is used to represent the fault time data.
#' If 0, no fault is detected at the end of interval.
#' @param te A numeric value for the time interval from the last fault to the
#' observation time.
#' @param data A dataframe. The arguments; time, fault, type, te can also be
#' selected as the columns of dataframe.
#' @param srms A list of estimated results.
#' @param xlab A character string, indicating the label of x-axis.
#' @param ylab A character string, indicating the label of y-axis.
#' @param datalab A character string, indicating the label of data in legend.
#' @param xmax A numeric value, indicating the maximum value of xlim.
#' If NA, the maximum value is automatically determined.
#' @param ymax A numeric value, indicating the maximum value of ylim.
#' If NA, the maximum value is automatically determined.
#' @param colors A vector of strings, indicating the colors in plot.
#' @param ... Other parameters.
#' @return An object of ggplot.
#' @examples
#' data(dacs)
#' result <- fit.srm.nhpp(fault=tohma, srm.names = c("exp", "gamma"),
#'             selection = NULL)
#' rateplot(fault=tohma, srms=result)
#' @export

rateplot <- function(time = NULL, fault = NULL, type = NULL, te = NULL,
                     data = data.frame(), srms = list(), xlab = "time",
                     ylab = "detection rate", datalab = "data", xmax = NA,
                     ymax = NA, colors = mmcolors, ...) {
  if (class(data) != "Rsrat.faultdata") {
    data <- faultdata.nhpp(substitute(time), substitute(fault),
                           substitute(type), substitute(te), data, parent.frame())
  }
  present <- sum(data$time)
  if (is.na(xmax)) {
    xmax <- present * 1.2
  }

  tt <- c()
  nn <- c()
  tmpt <- 0
  tmpn <- 0
  for (i in 1:length(data$time)) {
    t <- data$time[i]
    n <- data$fault[i] + data$type[i]
    if (t == 0 && n == 0) {
      continue
    }
    if (t == 0) {
      if (tmpt == 0) {
        nn[length(nn)] <- nn[length(nn)] + n
      } else {
        tmpn <- tmpn + n
      }
    } else if (n == 0) {
      tmpt <- tmpt + t
    } else {
      tt <- c(tt, t + tmpt)
      nn <- c(nn, n + tmpn)
      tmpt <- 0
      tmpn <- 0
    }
  }
  data <- data.frame(x=as.numeric(cumsum(tt)), y= nn / tt)
  gp <- ggplot(data) + labs(x=xlab, y=ylab) + xlim(c(0.1,xmax)) + ylim(c(0,ymax))
  gp <- gp + geom_area(stat="identity", position="identity", alpha=0.3, aes_string(x="x", y="y"))
  gp <- gp + geom_vline(xintercept=present, linetype="dotted")
  if ("list" %in% class(srms)) {
    for (s in srms) {
      gp <- gp + stat_function(fun=rate, args=list(srm=s), aes_(colour=srmname(s)))
    }
    gp + scale_colour_manual(NULL, values = colors)
  } else {
    gp <- gp + stat_function(fun=rate, args=list(srm=srms), aes_(colour=srmname(srms)))
    gp + scale_colour_manual(NULL, values = colors)
  }
}

mmcolors <- c(
  "#5E81B5",
  "#E19C24",
  "#8FB032",
  "#EB6235",
  "#8778B3",
  "#C56E1A",
  "#5D9EC7",
  "#FFBF00",
  "#A5609D",
  "#929600",
  "#E95536",
  "#6685D9",
  "#F89F13",
  "#BC5B80",
  "#47B66D"
)
