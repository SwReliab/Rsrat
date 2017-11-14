#' Software fault data with d-metrics
#'
#' Function faultdata.dmet() creates a list to store the fault data and dymanic
#' metrics (d-metrics)that are used to esiamte model parameters of SRM.
#'
#' @param formula An object of class formula. A symbolic description of the
#' model to be fitted. The output variable should be the column for the number
#' of faults.
#' @param data A dataframe for d-metrics and the number of faults.
#' @param offset An integer. This can be used to specify an a priori known
#' component to be included in the linear predictor during fitting. This should
#' be NULL or a numeric vector of length equal to the number of cases.
#' @return A list with the attribute class='Rsrat.faultdata';
#' \item{fault}{A vector for the number of detected faults.}
#' \item{metrics}{A matrix for metrics.}
#' \item{total}{An integer for the number of total faults.}
#' \item{nmetrics}{An integer for the number of types of metrics.}
#' \item{offset}{A vector for offset.}
#' \item{len}{An integer for the number of data records.}
#' @examples
#' data(dmet)
#' faultdata.dmet(fault~., dmet.ds1)
#' @export

faultdata.dmet <- function(formula, data, offset = NULL) {
  f <- model.frame(formula, data)
  fault <- model.response(f)
  stopifnot(is.vector(fault))
  metrics <- model.matrix(formula, data)
  total <- sum(fault)
  nmetrics <- ncol(metrics)
  if (is.null(offset)) {
    offset <- rep(0, length(fault))
  }

  result <- list(
    fault=fault,
    metrics=metrics,
    total=total,
    nmetrics=nmetrics,
    offset=offset,
    len=length(fault)
  )
  class(result) <- "Rsrat.faultdata.dmet"
  result
}

#' Printing software fault data with d-metrics
#'
#' Print a data frame.
#'
#' @param x An object of Rsrat.faultdata.
#' @param ... Other parameters
#' @param digits The minimum number of significant digits.
#' @param quote A logical, indicating whether or not entries are printed with quotes.
#' @param right A logical, indicating whether or not strings are right-aligned.
#' @param row.names A logical or a character vector, indicating whether or not row names are printed.
#' @details This function calls print.data.frame to forms the metrics data.
#' @export

print.Rsrat.faultdata.dmet <- function (x, ..., digits = NULL, quote = FALSE,
  right = TRUE, row.names = TRUE) {
    df <- data.frame(fault=x$fault, x$metrics)
    print.data.frame(df, ..., digits=digits, quote=quote,
      right=right, row.names=row.names)
    invisible(x)
}
