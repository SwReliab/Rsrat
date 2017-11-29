#' EM step
#'
#' Execute one EM step for software reliability models.
#'
#' @name em
#' @param params A numeric vector for model parameters
#' @param data A faultdata
#' @param divide An integer for the number of integration points.
#' @param eps A numeric value for tolerance error.
#' @param ... Other parameters
#' @return A list with an updated parameter vector (param),
#' absolute difference of parameter vector (pdiff),
#' log-likelihood function for a given parameter vector (llf),
#' the number of total faults (total).
NULL
