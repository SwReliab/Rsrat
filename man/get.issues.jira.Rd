% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scraping-jira.R
\name{get.issues.jira}
\alias{get.issues.jira}
\title{Issues in JIRA}
\usage{
get.issues.jira(
  query,
  fields = "created",
  user = NULL,
  password = "",
  authtype = "basic",
  url = "https://issues.apache.org/jira/rest/api/2/search",
  startAt = 0,
  maxResults = 50
)
}
\arguments{
\item{query}{A character string for a query}

\item{fields}{A character string vector, indicating fields to be downloaded}

\item{user}{A string for user}

\item{password}{A string for password}

\item{authtype}{A string for authentication. The default is \code{basic}.}

\item{url}{An URL for JIRA}

\item{startAt}{An integer to start}

\item{maxResults}{An integer to the number of issues for a transaction.}
}
\description{
Collect issues from JIRA with a query
}
