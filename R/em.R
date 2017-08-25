#' EM step
#'
#' Execute one EM step for software reliability models.
#'
#' @name em
#' @param params A numeric vector for model parameters
#' @param data A faultdata
#' @param divide An integer for the number of integration points.
#' @param ... Other parameters
#' @return A list with an updated parameter vector (param),
#' absolute difference of parameter vector (pdiff),
#' log-likelihood function for a given parameter vector (llf),
#' the number of total faults (total).
NULL
#> NULL

#' @rdname em
em.exp <- function(params, data, ...) {
  stopifnot(all(length(data$time) == data$len,
    length(data$type) == data$len, length(data$type) == data$len))
  stopifnot(length(params) == 2)
  res <- .C("em_exp_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,length(params))),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

#' @rdname em
em.gamma <- function(params, data, divide, ...) {
  stopifnot(all(length(data$time) == data$len,
    length(data$type) == data$len, length(data$type) == data$len))
  stopifnot(length(params) == 3)
  res <- .C("em_gamma_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            divide=as.integer(divide),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,length(params))),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

#' @rdname em
em.pareto <- function(params, data, ...) {
  stopifnot(all(length(data$time) == data$len,
    length(data$type) == data$len, length(data$type) == data$len))
  stopifnot(length(params) == 3)
  res <- .C("em_pareto_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,length(params))),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

#' @rdname em
em.tnorm <- function(params, data, ...) {
  stopifnot(all(length(data$time) == data$len,
    length(data$type) == data$len, length(data$type) == data$len))
  stopifnot(length(params) == 3)
  res <- .C("em_tnorm_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,length(params))),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

#' @rdname em
em.lnorm <- function(params, data, ...) {
  stopifnot(all(length(data$time) == data$len,
    length(data$type) == data$len, length(data$type) == data$len))
  stopifnot(length(params) == 3)
  res <- .C("em_lnorm_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,length(params))),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

#' @rdname em
em.tlogist <- function(params, data, ...) {
  stopifnot(all(length(data$time) == data$len,
    length(data$type) == data$len, length(data$type) == data$len))
  stopifnot(length(params) == 3)
  res <- .C("em_tlogist_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,length(params))),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

#' @rdname em
em.llogist <- function(params, data, ...) {
  stopifnot(all(length(data$time) == data$len,
    length(data$type) == data$len, length(data$type) == data$len))
  stopifnot(length(params) == 3)
  res <- .C("em_llogist_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,length(params))),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

#' @rdname em
em.txvmax <- function(params, data, ...) {
  stopifnot(all(length(data$time) == data$len,
    length(data$type) == data$len, length(data$type) == data$len))
  stopifnot(length(params) == 3)
  res <- .C("em_txvmax_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,length(params))),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

#' @rdname em
em.lxvmax <- function(params, data, ...) {
  stopifnot(all(length(data$time) == data$len,
    length(data$type) == data$len, length(data$type) == data$len))
  stopifnot(length(params) == 3)
  res <- .C("em_lxvmax_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,length(params))),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

#' @rdname em
em.txvmin <- function(params, data, ...) {
  stopifnot(all(length(data$time) == data$len,
    length(data$type) == data$len, length(data$type) == data$len))
  stopifnot(length(params) == 3)
  res <- .C("em_txvmin_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,length(params))),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

#' @rdname em
em.lxvmin <- function(params, data, ...) {
  stopifnot(all(length(data$time) == data$len,
    length(data$type) == data$len, length(data$type) == data$len))
  stopifnot(length(params) == 3)
  res <- .C("em_lxvmin_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,length(params))),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}
