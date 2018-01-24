#' Issues in JIRA
#'
#' Collect issues from JIRA with a query
#'
#' @param query A character string for a query
#' @param fields A character string vector, indicating fields to be downloaded
#' @param user A string for user
#' @param password A string for password
#' @param authtype A string for authentication. The default is \code{basic}.
#' @param url An URL for JIRA
#' @param startAt An integer to start
#' @param maxResults An integer to the number of issues for a transaction.
#' @examples
#' ## get crated date of blocer, critical and major issues
#' ## which affected to 1.0.0 of ACE from Apache JIRA
#' ## query <- "project = ACE AND issuetype = Bug"
#' ## query <- paste(query, "AND priority in (Blocker, Critical, Major)")
#' ## query <- paste(query, "AND affectedVersion = 1.0.0")
#' ## get.issues.jira(query = query, fields="created")
#' @export

get.issues.jira <- function(query, fields = "created", user = NULL, password = "", authtype = "basic",
  url = "https://issues.apache.org/jira/rest/api/2/search", startAt = 0, maxResults = 50) {
  url <- paste(url, "?jql=", URLencode(query, reserved = TRUE), "&fields=", paste(fields, collapse=",", sep=""), sep="")
  cat(sprintf("URL: %s\n", url))
  issues <- c()
  repeat {
    tmpurl <- paste(url, "&startAt=", startAt, "&maxResults=", maxResults, sep="")
    # print(tmpurl)
    cat(sprintf("  Getting (%d - %d) ...\n", startAt, startAt + maxResults - 1))
    if (is.null(user)) {
      rec <- httr::content(httr::GET(tmpurl))
    }
    else {
      rec <- httr::content(httr::GET(tmpurl, httr::authenticate(user, password, authtype)))
    }
    startAt <- startAt + maxResults
    issues <- c(issues, rec$issues)
    if (rec$total < startAt) {
      break
    }
  }
  issues
}

#' get open date from JIRA
#'
#' Collect open date of issues from JIRA with query
#'
#' @param query A character string for a query
#' @param user A string for user
#' @param password A string for password
#' @param authtype A string for authentication. The default is \code{basic}.
#' @param url An URL for JIRA
#' @param startAt An integer to start
#' @param maxResults An integer to the number of issues for a transaction.
#' @param date.min A character string for date, indicating the starting date.
#' @param date.max A character string for date, indicating the ending date.
#' @param by A character string. A increment of date. "1 month", "1 day", etc. are available.
#' @export

get.opendate.jira <- function(query, user = NULL, password = "", authtype = "basic",
  url = "https://issues.apache.org/jira/rest/api/2/search", startAt = 0, maxResults = 50,
  date.min = NULL, date.max = Sys.time(), by = "1 month") {
  data <- get.issues.jira(query, "created", user, password, authtype, url, startAt, maxResults)
  count.date(date=sapply(data, function(x) as.character(as.Date(x$fields$created))), date.min, date.max, by)
}

#' get open and closed date from JIRA
#'
#' Collect open and closed date of issues from JIRA with query
#'
#' @param query A character string for a query
#' @param user A string for user
#' @param password A string for password
#' @param authtype A string for authentication. The default is \code{basic}.
#' @param url An URL for JIRA
#' @param startAt An integer to start
#' @param maxResults An integer to the number of issues for a transaction.
#' @param date.min A character string for date, indicating the starting date.
#' @param date.max A character string for date, indicating the ending date.
#' @param by A character string. A increment of date. "1 month", "1 day", etc. are available.
#' @export

get.openclosedate.jira <- function(query, user = NULL, password = "", authtype = "basic",
  url = "https://issues.apache.org/jira/rest/api/2/search", startAt = 0, maxResults = 50,
  date.min = NULL, date.max = Sys.time(), by = "1 month") {
  data <- get.issues.jira(query, c("created", "resolutiondate"), user, password, authtype, url, startAt, maxResults)
  ddate <- sapply(data, function(x) as.character(as.Date(x$fields$created)))
  cdate <- sapply(data, function(x) if (is.null(x$fields$resolutiondate)) { NA } else { as.character(as.Date(x$fields$resolutiondate)) })
  count.date2(ddate, cdate, date.min, date.max, by)
}

#' Issues in GitHub
#'
#' Collect issues from GitHub with a query
#'
#' @param repository A character string indicating a repository such as "hadley/devtools"
#' @param query A character string for a query
#' @param user A string for user
#' @param password A string for password
#' @param authtype A string for authentication. The default is \code{basic}.
#' @param url An URL for GitHub
#' @param startPage An integer to identify the start page
#' @param perPage An integer to the number of issues for a page
#' @examples
#' ## get crated date of blocer, critical and major issues
#' ## which affected to 1.0.0 of ACE from Apache JIRA
#' ## query <- "project = ACE AND issuetype = Bug"
#' ## query <- paste(query, "AND priority in (Blocker, Critical, Major)")
#' ## query <- paste(query, "AND affectedVersion = 1.0.0")
#' ## get.issues.jira(query = query, fields="created")
#' @export

get.issues.github <- function(repository, query, user = NULL, password = "", authtype = "basic",
  url = "https://api.github.com/repos", startPage = 1L, perPage = 30L) {
  url <- paste(url, "/", repository, "?", URLencode(query, reserved = TRUE), sep="")
  cat(sprintf("URL: %s\n", url))
  issues <- c()
  repeat {
    tmpurl <- paste(url, URLencode(paste("&page=", startPage, "&per_page=", perPage, sep=""), reserved = TRUE), sep="")
    return(tmpurl)
    cat(sprintf("  Getting page %d ...\n", startPage))
    if (is.null(user)) {
      rec <- httr::content(httr::GET(tmpurl))
    }
    else {
      rec <- httr::content(httr::GET(tmpurl, httr::authenticate(user, password, authtype)))
    }
    # if (rec$documentation_url == "https://developer.github.com/v3/#rate-limiting") {
    #   stop(rec$message)
    # }
    if (length(rec) == 0) {
      break
    }
    startPage <- startPage + 1
    issues <- c(issues, rec)
  }
  issues
}
