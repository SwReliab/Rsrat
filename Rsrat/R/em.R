

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
