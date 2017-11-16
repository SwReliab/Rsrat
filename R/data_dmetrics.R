#' Software fault data
#'
#' Software fault data with d-metrics (dynamic metrics)
#'
#' This data set are software fault data with dynamic metrics such as the number of test cases, and so on.
#' in real embedded software which is control programs for a printer.
#' \describe{
#'  \item{ds1}{The fault data for control programs for a printer.
#'     The number of total faults is 66, and the number of testing periods is 20.
#'     Metrics are;
#'     \describe{
#'       \item{day}{The day of testing period.}
#'       \item{fault}{The number of faults detected at the testing period.}
#'       \item{tc}{The number of test cases executed at the testing period.}
#'       \item{ctc}{The cumulative number of test cases executed by the testing period.}
#'       \item{cov}{The gain of C0 coverage at the testing period.}
#'       \item{ccov}{The C0 coverage by the testing period.}
#'     }}
#'  \item{ds2}{The fault data for control programs for a printer.
#'     The number of total faults is 58, and the number of testing periods is 33.
#'     Metrics are;
#'     \describe{
#'       \item{day}{The day of testing period.}
#'       \item{fault}{The number of faults detected at the testing period.}
#'       \item{tc}{The number of test cases executed at the testing period.}
#'       \item{ctc}{The cumulative number of test cases executed by the testing period.}
#'       \item{cov}{The gain of C0 coverage at the testing period.}
#'       \item{ccov}{The C0 coverage by the testing period.}
#'     }}
#'  \item{ds3}{The fault data for control programs for a printer.
#'     The number of total faults is 52, and the number of testing periods is 30.
#'     Metrics are;
#'     \describe{
#'       \item{day}{The day of testing period.}
#'       \item{fault}{The number of faults detected at the testing period.}
#'       \item{tc}{The number of test cases executed at the testing period.}
#'       \item{ctc}{The cumulative number of test cases executed by the testing period.}
#'       \item{cov}{The gain of C0 coverage at the testing period.}
#'       \item{ccov}{The C0 coverage by the testing period.}
#'     }}
#' }
#'
#' @name dmet
#' @docType data
#' @references H. Okamura, Y. Etani, and T. Dohi (2010) A multi-factor software reliability model
#' based on logistic regression, Proceedings of 21st International Symposium on Software Reliability
#' Engineering (ISSRE2010), 31-40.
NULL

#' @name dmet.ds1
#' @rdname dmet
NULL

#' @name dmet.ds2
#' @rdname dmet
NULL

#' @name dmet.ds3
#' @rdname dmet
NULL
