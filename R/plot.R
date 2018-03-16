#' Draw a graph of mean value functions
#'
#' This function draws a graph of mean value function with data and NHPP classes.
#'
#' The mvf argument should be a list of NHPP classes. In the function, \code{mvf}
#' and \code{name} attributes are used. Therefore, the list with \code{mvf} and
#' \code{name} can be used as an element of the list of mvf argument.
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
#' @param mvf A list of NHPP classes. See details.
#' @param xlab A character string, indicating the label of x-axis.
#' @param ylab A character string, indicating the label of y-axis.
#' @param datalab A character string, indicating the label of data in legend.
#' @param xmax A numeric value, indicating the maximum value of xlim. If NA, the maximum value is automatically determined.
#' @param ymax A numeric value, indicating the maximum value of ylim. If NA, the maximum value is automatically determined.
#' @param colors A vector of strings, indicating the colors in plot.
#' @param ... Other parameters.
#' @return An object of ggplot.
#' @examples
#' data(tomcat5)
#'
#' result <- fit.srm.nhpp(time=tomcat5.catalina$time, fault=tomcat5.catalina$fault,
#'             srm.names = c("exp", "gamma"), selection = "all")
#'
#' mvfplot(time=tomcat5.catalina$time, fault=tomcat5.catalina$fault,
#'      mvf=lapply(result, function(s) s$srm))
#' @export

mvfplot <- function(time = NULL, fault = NULL, type = NULL, te = NULL, data = data.frame(),
  mvf = list(), xlab = "time", ylab = "# of faults", datalab = "data",
  xmax = NA, ymax = NA, colors = mmcolors, ...) {
  if (class(data) != "Rsrat.faultdata") {
    data <- .faultdata.nhpp(substitute(time), substitute(fault),
      substitute(type), substitute(te), data, parent.frame())
  }
  n <- data$fault + data$type
  data <- data.frame(x=cumsum(data$time)[n != 0], y=cumsum(n)[n != 0])
  gp <- ggplot(data, aes_string(x="x", y="y")) + labs(x=xlab, y=ylab) + xlim(c(0,xmax)) + ylim(c(0,ymax))
  gp <- gp + geom_point(aes_(colour=datalab))

  for (s in mvf) {
    gp <- gp + stat_function(fun=s$mvf, aes_(colour=s$name))
  }
  gp + scale_colour_manual(NULL, values = colors)
}

#' Draw a graph of diff of mean value functions on discrete time domain
#'
#' This function draws a graph of mean value function with data and NHPP classes.
#'
#' The dmvf argument should be a list of NHPP classes. In the function, \code{dmvf}
#' and \code{name} attributes are used. Therefore, the list with \code{dmvf} and
#' \code{name} can be used as an element of the list of dmvf argument.
#' The return value is an object of ggplot2, and thus a line is added like
#' \code{dmvfplot(...) + geom_point(stat="identity", position="identity", aes_string(x="x", y="exp", colour=shQuote("exp")))
#'  + geom_line(stat="identity", position="identity", aes_string(x="x", y="exp", colour=shQuote("exp")))}.
#'
#' @param time A numeric vector for time intervals.
#' @param fault An integer vector for the number of faults detected in time intervals.
#' The fault detected just at the end of time interal is not counted.
#' @param type Either 0 or 1. If 1, a fault is detected just at the end of corresponding time interval.
#' This is used to represent the fault time data. If 0, no fault is detected at the end of interval.
#' @param te A numeric value for the time interval from the last fault to the observation time.
#' @param data A dataframe. The arguments; time, fault, type, te can also be selected as the columns of dataframe.
#' @param dmvf A list of NHPP classes. See details.
#' @param xlab A character string, indicating the label of x-axis.
#' @param ylab A character string, indicating the label of y-axis.
#' @param datalab A character string, indicating the label of data in legend.
#' @param xmax A numeric value, indicating the maximum value of xlim. If NA, the maximum value is automatically determined.
#' @param ymax A numeric value, indicating the maximum value of ylim. If NA, the maximum value is automatically determined.
#' @param colors A vector of strings, indicating the colors in plot.
#' @param ... Other parameters.
#' @return An object of ggplot.
#' @examples
#' data(dmet)
#'
#' result <- fit.srm.logit(fault~., data=dmet.ds1)
#'
#' dmvfplot(fault=dmet.ds1$fault, dmvf=list(result$srm))
#' @export

dmvfplot <- function(time = NULL, fault = NULL, type = NULL, te = NULL, data = data.frame(),
  dmvf = list(), xlab = "time", ylab = "# of faults", datalab = "data",
  xmax = NA, ymax = NA, colors = mmcolors, ...) {
  if (class(data) != "Rsrat.faultdata") {
    data <- .faultdata.nhpp(substitute(time), substitute(fault),
      substitute(type), substitute(te), data, parent.frame())
  }
  t <- as.numeric(cumsum(data$time))
  n <- as.numeric(data$fault + data$type)

  data <- data.frame(x=t, y=n)
  for (s in dmvf) {
    data <- cbind(data, s$dmvf(t))
  }
  colnames(data) <- c("x", "y", sapply(dmvf, function(s) s$name))

  gp <- ggplot(data, aes_string(x="x", y="y")) + labs(x=xlab, y=ylab) + xlim(c(0,xmax)) + ylim(c(0,ymax))
  gp <- gp + geom_bar(stat="identity", position="identity", alpha=0.5, aes_(fill=datalab))
  for (s in dmvf) {
    gp <- gp + geom_point(stat="identity", position="identity", aes_string(x="x", y=s$name, colour=shQuote(s$name)))
    gp <- gp + geom_line(stat="identity", position="identity", aes_string(x="x", y=s$name, colour=shQuote(s$name)))
  }
  gp + scale_colour_manual(NULL, values = colors[2:length(colors)]) + scale_fill_manual(NULL, values = colors[1])
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
