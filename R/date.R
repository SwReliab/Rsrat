#' count date
#'
#' Make grouped data for detection data.
#'
#' @param date A vector of character strings for dates.
#' @param date.min A character string for date, indicating the starting date.
#' @param date.max A character string for date, indicating the ending date.
#' @param by A character string. A increment of date. "1 month", "1 day", etc. are available.
#' @export

count.date <- function(date, date.min = NULL, date.max = Sys.time(), by = "1 month") {
 dd <- as.Date(date)
 date.max <- as.Date(date.max)
 if (is.null(date.min)) {
   date.min <- as.Date(cut(min(dd), breaks=by))
 } else {
   date.min <- as.Date(date.min)
   date.min <- as.Date(cut(date.min, breaks=by))
 }
 v <- seq(from=date.min, to=date.max, by=by)
 if (v[length(v)] != date.max) {
   v <- c(v, date.max+1)
 }
 dl <- numeric(0)
 tt <- numeric(0)
 res <- numeric(length(v)-1)
 for (i in 1:(length(v)-1)) {
   l <- sum(sapply(dd, function(x) (v[i] <= x) && (x < v[i+1])))
   dl <- c(dl, format(v[i], "%Y-%m-%d"))
   tt <- c(tt, difftime(v[i+1], v[i]))
   res[i] <- l
 }
 list(date=dl, time=tt, counts=res, data.min=date.min, date.max=date.max)
}

#' count date 2
#'
#' Make grouped data for detection and correction.
#'
#' @param ddate A vector of character strings for detection dates.
#' @param cdate A vector of character strings for correction dates.
#' @param date.min A character string for date, indicating the starting date.
#' @param date.max A character string for date, indicating the ending date.
#' @param by A character string. A increment of date. "1 month", "1 day", etc. are available.
#' @note
#' The format of date is yyyy-mm-dd. NA means that the correction has not been done by the max date.
#' @export

count.date2 <- function(ddate, cdate, date.min = NULL, date.max = Sys.time(), by = "1 month") {
 dd <- as.Date(ddate)
 dc <- as.Date(cdate)
 da <- c(dd, dc[!is.na(dc)])
 date.max <- as.Date(date.max)
 if (is.null(date.min)) {
   date.min <- as.Date(cut(min(da), breaks=by))
 } else {
   date.min <- as.Date(date.min)
   date.min <- as.Date(cut(date.min, breaks=by))
 }
 v <- seq(from=date.min, to=date.max, by=by)
 if (v[length(v)] != date.max) {
   v <- c(v, date.max+1)
 }
 dl <- numeric(0)
 tt <- numeric(0)
 res <- Matrix::Matrix(0, length(v)-1, length(v))
 for (i in 1:(length(v)-1)) {
   df <- sapply(dd, function(x) (v[i] <= x) && (x < v[i+1]))
   tmp <- dc[df]
   res[i,length(v)] <- sum(is.na(tmp))
   tmp <- tmp[!is.na(tmp)]
   for (j in 1:(length(v)-1)) {
     if (length(tmp) >= 1) {
       res[i,j] <- sum(sapply(tmp, function(x) (v[j] <= x) && (x < v[j+1])))
     }
   }
   dl <- c(dl, format(v[i], "%Y-%m-%d"))
   tt <- c(tt, difftime(v[i+1], v[i]))
 }
 list(date=dl, time=tt, counts=res, data.min=date.min, date.max=date.max)
}
