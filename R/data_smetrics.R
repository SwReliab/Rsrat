#' Software fault data
#'
#' Software fault data with s-metrics (static metrics)
#'
#' This data set are software fault data with static metrics such as LOC and complexity, and so on.
#' in Apache Tomcat5 consisting of the following modules:
#' \describe{
#'  \item{catalina}{The servlet container core.}
#'  \item{connector}{Coyote connectors.}
#'  \item{jasper}{The JSP (JavaServer Pages) compiler and runtime engine.}
#'  \item{servlets}{Servlet API and support programs for CGI, SSI and WebDAV.}
#'  \item{tester}{Unit testing framework.}
#'  \item{webapps}{Web application for administration, documentation and examples.}
#' }
#' The fault data are
#' \describe{
#'  \item{tomcat5.smet}{The table of code metrics;
#'     \describe{
#'       \item{LOC}{Lines of code.}
#'       \item{St}{The number of statements.}
#'       \item{Br}{Percent branch statements.}
#'       \item{Co}{Percent lines with comments.}
#'       \item{Fn}{The number of functions}
#'       \item{Mc}{The maximum McCabe compleixity.}
#'       \item{Ac}{The average McCabe compliexity.}
#'     }}
#'  \item{tomcat5.catalina}{The number of bugs detected for each month in catalina.}
#'  \item{tomcat5.connector}{The number of bugs detected for each month in connector.}
#'  \item{tomcat5.jasper}{The number of bugs detected for each month in jasper.}
#'  \item{tomcat5.servlets}{The number of bugs detected for each month in servlets.}
#'  \item{tomcat5.tester}{The number of bugs detected for each month in tester.}
#'  \item{tomcat5.webapps}{The number of bugs detected for each month in webapps.}
#' }
#'
#' @name tomcat5
#' @docType data
#' @references H. Okamura and T. Dohi (2014), A novel framework of software reliability
#' evaluation with software reliability growth models and software metrics,
#' Proceedings of The 15th IEEE International Symposium on High Assurance Systems
#' Engineering (HASE 2014), 97-104.
NULL

#' @name tomcat5.smet
#' @rdname tomcat5
NULL

#' @name tomcat5.catalina
#' @rdname tomcat5
NULL

#' @name tomcat5.connector
#' @rdname tomcat5
NULL

#' @name tomcat5.jasper
#' @rdname tomcat5
NULL

#' @name tomcat5.servlets
#' @rdname tomcat5
NULL

#' @name tomcat5.tester
#' @rdname tomcat5
NULL

#' @name tomcat5.webapps
#' @rdname tomcat5
NULL
