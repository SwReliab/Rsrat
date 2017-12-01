#' Rsrat: A Software Reliability Assessment Tool
#'
#' This package provides estimation programs for software reliability growth
#' models, metrics-based software reliability growth models with logistic
#' regression for dynamic metrics and Poisson regression for static metrics.
#'
#' @docType package
#' @name Rsrat
#' @import R6 ggplot2
#' @importFrom stats dlogis dnorm plogis pnorm qlogis qnorm nobs model.frame model.response model.matrix
#' @importFrom utils URLencode
#' @importFrom RCurl getURL
#' @importFrom rjson fromJSON
#' @useDynLib Rsrat
NULL
