#' Issues in GitHub
#'
#' Collect issues from GitHub with a query
#'
#' @param repo A character string indicating a repository such as "hadley/devtools"
#' @param query A character string for a query
#' @param user A string for user
#' @param password A string for password or personal access token
#' @param authtype A string for authentication. The default is \code{token} (personal access toke).
#' @param url An URL for GitHub
#' @param startPage An integer to identify the start page
#' @param perPage An integer to the number of issues for a page
#' @export

get.issues.github <- function(repo, query = "", user = NULL, password = "", authtype = "token",
                              url = "https://api.github.com/repos", startPage = 1L, perPage = 100L) {
  url <- paste(url, "/", repo, "/issues?", URLencode(query), sep="")
  cat(sprintf("URL: %s\n", url))
  issues <- c()
  repeat {
    if (query == "") {
      tmpurl <- paste(url, URLencode(paste("page=", startPage, "&per_page=", perPage, sep="")), sep="")
    }
    else {
      tmpurl <- paste(url, URLencode(paste("&page=", startPage, "&per_page=", perPage, sep="")), sep="")
    }
    cat(sprintf("  Getting page %d ...\n", startPage))
    rec <- switch(authtype,
                  "basic" = if (is.null(user)) {
                    httr::content(httr::GET(tmpurl))
                  } else {
                    httr::content(httr::GET(tmpurl, httr::authenticate(user, password, authtype)))
                  },
                  "token" = httr::content(httr::GET(tmpurl, add_headers(Authorization = paste("Bearer ", password, sep="")))),
                  NULL)
    if ("documentation_url" %in% names(rec)) {
      stop(rec$message)
    }
    if (length(rec) == 0) {
      break
    }
    startPage <- startPage + 1
    issues <- c(issues, rec)
  }
  if (is.null(issues)) {
    stop("Cannot get issues. Please check your query.")
  }
  issues
}

#' get open date from GitHub
#'
#' Collect open date of issues from GitHub with query
#'
#' @param repo A character string indicating a repository such as "hadley/devtools"
#' @param query A character string for a query
#' @param user A string for user
#' @param password A string for password or personal access token
#' @param authtype A string for authentication. The default is \code{token} (personal access toke).
#' @param url An URL for GitHub
#' @param startPage An integer to identify the start page
#' @param perPage An integer to the number of issues for a page
#' @param date.min A character string for date, indicating the starting date.
#' @param date.max A character string for date, indicating the ending date.
#' @param by A character string. A increment of date. "1 month", "1 day", etc. are available.
#' @export

get.opendate.github <- function(repo, query = "state=all&labels=bug", user = NULL, password = "", authtype = "token",
  url = "https://api.github.com/repos", startPage = 1L, perPage = 100L,
  date.min = NULL, date.max = Sys.time(), by = "1 month") {
  data <- get.issues.github(repo, query, user, password, authtype, url, startPage, perPage)
  ddate <- sapply(data, function(x) as.character(as.Date(x$created_at)))
  count.date(ddate, date.min, date.max, by)
}

#' get open and closed date from GitHub
#'
#' Collect open and closed date of issues from GitHub with query
#'
#' @param repo A character string indicating a repository such as "hadley/devtools"
#' @param query A character string for a query
#' @param user A string for user
#' @param password A string for password or personal access token
#' @param authtype A string for authentication. The default is \code{token} (personal access toke).
#' @param url An URL for GitHub
#' @param startPage An integer to identify the start page
#' @param perPage An integer to the number of issues for a page
#' @param date.min A character string for date, indicating the starting date.
#' @param date.max A character string for date, indicating the ending date.
#' @param by A character string. A increment of date. "1 month", "1 day", etc. are available.
#' @export

get.openclosedate.github <- function(repo, query = "state=all&labels=bug", user = NULL, password = "", authtype = "token",
  url = "https://api.github.com/repos", startPage = 1L, perPage = 100L,
  date.min = NULL, date.max = Sys.time(), by = "1 month") {
  data <- get.issues.github(repo, query, user, password, authtype, url, startPage, perPage)
  ddate <- sapply(data, function(x) as.character(as.Date(x$created_at)))
  cdate <- sapply(data, function(x) if (is.null(x$closed_at)) { NA } else { as.character(as.Date(x$closed_at)) })
  count.date2(ddate, cdate, date.min, date.max, by)
}
