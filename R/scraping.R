#' Issues in JIRA
#'
#' Collect issues from JIRA with a query
#'
#' @param query A character string for a query
#' @param userpwd A string for user:passwd
#' @param url An URL for JIRA
#' @param startAt An integer to start
#' @param maxResults An integer to the number of issues for a transaction.
#' @examples
#' ## get issuetype, priority, created, resolutiondate, resolution and versions of issues of ACE from Apache JIRA
#' get.issues.jira("project=ACE&fields=id,key,issuetype,priority,created,resolutiondate,resolution,versions",
#'     userpwd = "username:password")
#' @export

get.issues.jira <- function(query, userpwd, url = "https://issues.apache.org/jira/rest/api/2/search",
  startAt = 0, maxResults = 50) {
  url <- paste(url, URLencode(paste("jql=", query, sep=""), reserved = FALSE), sep="?")
  cat(sprintf("URL: %s\n", url))
  issues <- c()
  repeat {
    tmpurl <- paste(url, "&startAt=", startAt, "&maxResults=", maxResults, sep="")
    cat(sprintf("  Getting (%d - %d) ...\n", startAt, startAt + maxResults - 1))
    rec <- fromJSON(getURL(tmpurl, userpwd = userpwd, httpheader = c('Content-Type' = "application/json")))
    startAt <- startAt + maxResults
    issues <- c(issues, rec$issues)
    if (rec$total < startAt) {
      break
    }
  }
  issues
}
