#' Software fault data with s-metrics
#'
#' Function faultdata.smet() creates a list to store the fault data and dymanic
#' metrics (d-metrics)that are used to esiamte model parameters of SRM.
#'
#' Since the row.names of data and names of srm.data are used to match the compoments,
#' it had better they all are labelled. If they are not labelled, the pairs are guessed
#' with their order.
#'
#' @param formula An object of class formula. A symbolic description of the
#' model to be fitted. The output variable should be empty.
#' @param data A dataframe for s-metrics.
#' @param names A list for the names of components. This is prioritized to the others.
#' @param offset An integer. This can be used to specify an a priori known
#' component to be included in the linear predictor during fitting. This should
#' be NULL or a numeric vector of length equal to the number of cases.
#' @return A list with the attribute class='Rsrat.faultdata.smet';
#' \item{m}{An integer for the number of components.}
#' \item{names}{A list for the names of components.}
#' \item{metrics}{A matrix for s-metrics.}
#' \item{nmetrics}{An integer for the number of types of metrics.}
#' \item{offset}{A vector for offset.}
#' @examples
#' data(tomcat5)
#' faultdata.smet(~., tomcat5.smet)
#' @export

faultdata.smet <- function(formula, data, names = NULL, offset = NULL) {
  if (is.null(names)) {
    if (is.null(row.names(data))) {
      names <- NULL
    } else {
      names <- row.names(data)
    }
  }
  stopifnot(length(names) <= nrow(data))
  m <- length(names)

  if (is.null(row.names(data))) {
    data <- data[1:m,]
    row.names(data) <- names
  } else {
    stopifnot(all(names %in% row.names(data)))
    data <- data[names,]
  }

  metrics <- model.matrix(formula, data)
  nmetrics <- ncol(metrics)
  if (is.null(offset)) {
    offset <- rep(0, m)
  }
  result <- list(
    names=names,
    metrics=metrics,
    nmetrics=nmetrics,
    m=m,
    offset=offset)
  class(result) <- "Rsrat.faultdata.smet"
  result
}

#' Printing software fault data with s-metrics
#'
#' Print data with s-metrics
#'
#' @param x An object of Rsrat.faultdata.
#' @param ... Other parameters
#' @param digits The minimum number of significant digits.
#' @param quote A logical, indicating whether or not entries are printed with quotes.
#' @param right A logical, indicating whether or not strings are right-aligned.
#' @param row.names A logical or a character vector, indicating whether or not row names are printed.
#' @details This function calls print.data.frame to forms the metrics data.
#' @export

print.Rsrat.faultdata.smet <- function (x, ..., digits = NULL, quote = FALSE,
  right = TRUE, row.names = TRUE) {
    print.data.frame(data.frame(x$metrics), ..., digits=digits, quote=quote,
      right=right, row.names=row.names)
    invisible(x)
}
