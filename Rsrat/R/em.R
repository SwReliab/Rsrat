
emfit.srm <- function(srm, data,
  initialize = TRUE,
  maxiter = 1000, rtol = 1.0e-6, atol = 1.0e-3,
  termination = "llf") {
    ## init
    if (initialize) {
      srm$init_params(data)
    }

    iter <- 1
    conv <- FALSE
    param <- srm$params

    ## termination
    if (termination == "parameter" || termination == "param") {
      term.fn <- function(res0, res1) {
        sdiff <- res1$llf - res0$llf
        para0 <- c(res0$param)
        para1 <- c(res1$param)
        aerror <- max(abs(para1 - para0))
        rerror <- aerror / max(abs(para0))
        return(c(aerror, rerror, sdiff))
      }
    } else if (termination == "llf") {
      term.fn <- function(res0, res1) {
        sdiff <- res1$llf - res0$llf
        aerror <- abs(res1$llf - res0$llf)
        rerror <- aerror / abs(res0$llf)
        return(c(aerror, rerror, sdiff))
      }
    } else {
      stop("wrong termination condition.")
    }

    ## repeat emstep
    res0 <- list(param=param, llf=-Inf)

    repeat {
      res1 <- srm$em(res0$param, data)
      error <- term.fn(res0, res1)

      if (!is.finite(res1$llf) || !all(is.finite(res1$param))) {
        warning(sprintf("LLF or param becomes +-Inf, NaN or NA: %s %d",
          srm$name, iter))
        res1 <- res0
        break
      }
      if (error[3] < 0) {
        warning(sprintf("LLF decreses: %s %d %e", srm$name,
          iter, error[3]))
      }
      if ((error[1] < atol) && (error[2] < rtol)) {
        conv <- TRUE
        break
      }
      if (iter >= maxiter) {
        warning("Did not converge to MLE by max iteration.")
        break
      }
      iter <- iter + 1
      res0 <- res1
    }

    result <- list(
      initial = param,
      srm = srm,
      llf = res1$llf,
      convergence = conv,
      iter = iter,
      aerror = error[1],
      rerror = error[2]
    )
    class(result) <- "Rsrat.result"
    result
}

em.exp <- function(params, data, ...) {
  res <- .C("em_exp_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,2)),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

em.gamma <- function(params, data, divide, ...) {
  res <- .C("em_gamma_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            divide=as.integer(divide),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,3)),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

em.pareto <- function(params, data, ...) {
  res <- .C("em_pareto_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,3)),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

em.tnorm <- function(params, data, ...) {
  res <- .C("em_tnorm_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,3)),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

em.lnorm <- function(params, data, ...) {
  res <- .C("em_lnorm_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,3)),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

em.tlogist <- function(params, data, ...) {
  res <- .C("em_tlogist_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,3)),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

em.llogist <- function(params, data, ...) {
  res <- .C("em_llogist_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,3)),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

em.txvmax <- function(params, data, ...) {
  res <- .C("em_txvmax_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,3)),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

em.lxvmax <- function(params, data, ...) {
  res <- .C("em_lxvmax_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,3)),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

em.txvmin <- function(params, data, ...) {
  res <- .C("em_txvmin_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,3)),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}

em.lxvmin <- function(params, data, ...) {
  res <- .C("em_lxvmin_emstep", PACKAGE="Rsrat",
            dsize=as.integer(data$len),
            time=as.double(data$time),
            num=as.double(data$fault),
            type=as.integer(data$type),
            npara=as.integer(length(params)),
            para=as.double(params),
            pdiff=as.double(array(0,3)),
            retllf=as.double(array(0,1)),
            total=as.double(array(0,1))
  )
  list(param=res$para, pdiff=res$pdiff, llf=res$retllf, total=res$total)
}
