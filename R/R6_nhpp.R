#' Class for NHPP-based software reliability model
#'
#' @docType class
#' @name NHPP
#' @format \code{\link{R6Class}} object.
#' @seealso \code{\link{srm}}
#' @export
NHPP <- R6::R6Class("NHPP",
  private = list(
    Ft = function(t, lower.tail = TRUE) { NA },
    invFt = function(p) { NA },
    ft = function(t) { NA },
    rf = function(n) { NA }
  ),
  public = list(
    #' @field name A character string for the name of model.
    name = NA,

    #' @field params A numeric vector for the model parameters.
    params = NA,

    #' @field df An integer for the degrees of freedom of the model.
    df = NA,

    #' @field data Data to esimate parameters.
    data = NA,

    #' @description
    #' Print model parameters.
    #' @param digits An integer to determine the number of digits for print.
    #' @param ... Others
    print = function(digits = max(3, getOption("digits") - 3), ...) {
      cat(gettextf("Model name: %s\n", self$name))
      print.default(format(self$params, digits = digits), print.gap = 2, quote = FALSE)
    },

    #' @description
    #' The number of total faults.
    #' @return The number of total faults.
    omega = function() { self$params[1L] },

    #' @description 
    #' The mean value function.
    #' @param t Time
    #' @return The expected number of faults detected before time t. 
    mvf = function(t) { self$omega() * private$Ft(t) },

    #' @description 
    #' The difference of mean value function on discrete time domain.
    #' @param t A vector for time series.
    #' @return A vector of the expected number of faults in each time interval.
    dmvf = function(t) {
      t <- c(0, t)
      diff(self$omega() * private$Ft(t))
    },

    #' @description 
    #' The inverse of mean value function.
    #' @param x The number of faults.
    #' @return The first passage time at which the mean value function exceeds x.
    inv_mvf = function(x) {
      p <- x / self$omega()
      sapply(p, function(pp) {
        if (pp > 0 && pp < 1) {
          private$invFt(pp)
        } else {
          NA
        }
      })
    },

    #' @description 
    #' The intensity function.
    #' @param t Time
    #' @return The intensity (the first derivative of mean value function) at time t.
    intensity = function(t) { self$omega() * private$ft(t) },

    #' @description 
    #' Generate samples of fault-detection times following the NHPP
    #' over [lower,upper). If diff is true, the method provides
    #' time interval of samples.
    #' @param diff A logical meaning the generated samples are time intervals or not.
    #' @param lower A value indicating the minimum of time domain.
    #' @param upper A value indicating the maximum of time domain.
    #' @return A vector of samples following NHPP.
    sample = function(diff = TRUE, lower = 0, upper = +Inf) {
      n <- rpois(1, self$omega())
      x <- private$rf(n)
      x <- sort(x[x>=lower & x<=upper])
      if (diff == TRUE) {
        diff(c(lower, x))
      } else {
        x
      }
    },

    #' @description 
    #' The reliability function.
    #' @param t A value of time from s.
    #' @param s A value of current time.
    #' @return The probability that no fault is detected during [s, t+s)
    reliab = function(t, s) { exp(-(self$mvf(t+s) - self$mvf(s))) },

    #' @description 
    #' The expected number of residual faults at time t.
    #' @param t A value of time.
    #' @return The expected residual number of faults at time t.
    residual = function(t) { self$omega() * private$Ft(t, lower.tail=FALSE) },

    #' @description 
    #' The fault-free probability at time t.
    #' @param t A value of time.
    #' @return The probability that software does not have any fault at time t.
    ffp = function(t) { exp(-self$residual(t)) },

    #' @description 
    #' The instantaneous MTBF at time t.
    #' @param t A value of time.
    #' @return The inverse of intensity at time t.
    imtbf = function(t) { 1 / self$intensity(t) },

    #' @description
    #' The cumulative MTBF at time t.
    #' @param t A value of time.
    #' @return the cumulative MTBF.
    cmtbf = function(t) { t / self$mvf(t) },

    #' @description 
    #' The time at which the software reliability attains the probability p
    #' from the orign s.
    #' @param s A value of origin time.
    #' @param p A value of probability.
    #' @return the time at which the software reiability become p.
    median = function(s, p = 0.5) {
      if (p > self$ffp(s)) {
        private$invFt((self$mvf(s) - log(p)) / self$omega())
      } else {
        NA
      }
    },

    #' @description 
    #' Set initial parameters from a given data.
    #' This is used to set the initial value for
    #' the fitting algorithm.
    #' @param data Data.
    init_params = function(data) { NA },

    #' @description 
    #' Set model parameters.
    #' @param params Parameters.
    set_params = function(params) { self$params <- params },

    #' @description
    #' Set data to be used in the fitting algorithm.
    #' @param data Data.
    set_data = function(data) { self$data <- data },

    #' @description 
    #' Execute an EM step.
    #' @param params Parameters.
    #' @param data Data.
    #' @return
    #' A list with the following
    #' \describe{
    #' \item{param}{Updated parameters.}
    #' \item{pdiff}{Absolute difference of parameter vector.}
    #' \item{llf}{Log-likelihood function for a given parameter vector.}
    #' \item{total}{The number of total faults.}
    #' }
    em = function(params, data) { NA },

    #' @description 
    #' The log-likelihood function for a given data.
    #' @param data Data.
    #' @return Log-likelihood function.
    llf = function(data) {
      n <- data$len
      te <- sum(data$time)
      ct1 <- cumsum(data$time)
      ct0 <- c(0,ct1)[1:n]

      tt <- ct1[data$type != 0]
      ct1 <- ct1[data$fault != 0]
      ct0 <- ct0[data$fault != 0]
      fn <- data$fault[data$fault != 0]

      retval <- data$total * log(self$omega())
      if (length(tt) != 0) {
        retval <- retval + sum(log(private$ft(tt)))
      }
      retval <- retval +
        sum(apply(cbind(fn, ct1, ct0), 1,
                  function(v) v[1] * log(private$Ft(v[2]) - private$Ft(v[3]))
                        - lgamma(v[1] + 1)))
      retval <- retval - self$omega() * private$Ft(te)
      retval
    },

    #' @description 
    #' Provides the absolute error and relative error from the
    #' the two results. The default is the difference of log-likelihood.
    #' @param res0 One result
    #' @param res1 Another result
    #' @return A vector of abusolute error, relative error and the difference.
    comp_error = function(res0, res1) {
      sdiff <- res1$llf - res0$llf
      aerror <- abs(res1$llf - res0$llf)
      rerror <- aerror / abs(res0$llf)
      c(aerror, rerror, sdiff)
    }
  )
)

#' @rdname NHPP
#' @export
ExpSRM <- R6::R6Class("ExpSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      pexp(t, rate=self$rate(), lower.tail=lower.tail)
    },
    invFt = function(p) { qexp(p, rate=self$rate()) },
    ft = function(t) { dexp(t, rate=self$rate()) },
    rf = function(n) { rexp(n, rate=self$rate()) }
  ),
  public = list(
    name = "exp",
    df = 2,
    rate = function() { self$params[2L] },
    initialize = function(omega = 1, rate = 1) {
      self$params <- c(omega, rate)
      names(self$params) <- c("omega", "rate")
    },
    init_params = function(data) {
      self$params <- c(data$total, 1.0/data$mean)
    },
    em = function(params, data) {
      em_exp_emstep(params, data)
    }
  )
)

#' @rdname NHPP
#' @export
GammaSRM <- R6::R6Class("GammaSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      pgamma(t, shape=self$shape(), rate=self$rate(), lower.tail=lower.tail)
    },
    invFt = function(p) { qgamma(p, shape=self$shape(), rate=self$rate()) },
    ft = function(t) { dgamma(t, shape=self$shape(), rate=self$rate()) },
    rf = function(n) { rgamma(n, shape=self$shape(), rate=self$rate()) }
  ),
  public = list(
    name = "gamma",
    df = 3,
    shape = function() { self$params[2L] },
    rate = function() { self$params[3L] },
    initialize = function(omega = 1, shape = 1, rate = 1) {
      self$params <- c(omega, shape, rate)
      names(self$params) <- c("omega", "shape", "rate")
    },
    init_params = function(data) {
      self$params <- c(data$total, 1.0, 1.0/data$mean)
    },
    em = function(params, data, divide = 15) {
      em_gamma_emstep(params, data, divide)
    }
  )
)

#' @rdname NHPP
#' @export
ParetoSRM <- R6::R6Class("ParetoSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      ppareto2(t, shape=self$shape(), scale=self$scale(),
        lower.tail=lower.tail)
    },
    invFt = function(p) {
      qpareto2(p, shape=self$shape(), scale=self$scale())
    },
    ft = function(t) {
      dpareto2(t, shape=self$shape(), scale=self$scale())
    },
    rf = function(n) {
      rpareto2(n, shape=self$shape(), scale=self$scale())
    }
  ),
  public = list(
    name = "pareto",
    df = 3,
    shape = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, shape = 1, scale = 1) {
      self$params <- c(omega, shape, scale)
      names(self$params) <- c("omega", "shape", "scale")
    },
    init_params = function(data) {
      self$params <- c(1.0, 1.0, 1.0)
    },
    em = function(params, data) {
      em_pareto_emstep(params, data)
    }
  )
)

#' @rdname NHPP
#' @export
TNormSRM <- R6::R6Class("TNormSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      ptnorm(t, mean=self$mean(), sd=self$sd(), lower.tail=lower.tail)
    },
    invFt = function(p) { qtnorm(p, mean=self$mean(), sd=self$sd()) },
    ft = function(t) { dtnorm(t, mean=self$mean(), sd=self$sd()) },
    rf = function(n) { rtnorm(n, mean=self$mean(), sd=self$sd()) }
  ),
  public = list(
    name = "tnorm",
    df = 3,
    mean = function() { self$params[2L] },
    sd = function() { self$params[3L] },
    initialize = function(omega = 1, mean = 0, sd = 1) {
      self$params <- c(omega, mean, sd)
      names(self$params) <- c("omega", "mean", "sd")
    },
    init_params = function(data) {
      self$params <- c(1.0, 0.0, data$mean)
    },
    em = function(params, data) {
     em_tnorm_emstep(params, data)
    }
  )
)

#' @rdname NHPP
#' @export
LNormSRM <- R6::R6Class("LNormSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      plnorm(t, meanlog=self$meanlog(), sdlog=self$sdlog(),
        lower.tail=lower.tail)
      },
    invFt = function(p) {
      qlnorm(p, meanlog=self$meanlog(), sdlog=self$sdlog())
    },
    ft = function(t) {
      dlnorm(t, meanlog=self$meanlog(), sdlog=self$sdlog())
    },
    rf = function(n) {
      rlnorm(n, meanlog=self$meanlog(), sdlog=self$sdlog())
    }
  ),
  public = list(
    name = "lnorm",
    df = 3,
    meanlog = function() { self$params[2L] },
    sdlog = function() { self$params[3L] },
    initialize = function(omega = 1, meanlog = 0, sdlog = 1) {
      self$params <- c(omega, meanlog, sdlog)
      names(self$params) <- c("omega", "meanlog", "sdlog")
    },
    init_params = function(data) {
      self$params <- c(1.0, 1.0, max(log(data$mean), 1.0))
    },
    em = function(params, data) {
      em_lnorm_emstep(params, data)
    }
  )
)

#' @rdname NHPP
#' @export
TLogisSRM <- R6::R6Class("TLogisSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      ptlogis(t, location=self$location(), scale=self$scale(),
        lower.tail=lower.tail)
    },
    invFt = function(p) {
      qtlogis(p, location=self$location(), scale=self$scale())
    },
    ft = function(t) {
      dtlogis(t, location=self$location(), scale=self$scale())
    },
    rf = function(n) {
      rtlogis(n, location=self$location(), scale=self$scale())
    }
  ),
  public = list(
    name = "tlogis",
    df = 3,
    location = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, location = 0, scale = 1) {
      self$params <- c(omega, location, scale)
      names(self$params) <- c("omega", "location", "scale")
    },
    init_params = function(data) {
      self$params <- c(1.0, 0.0, data$mean)
    },
    em = function(params, data) {
      # em_tlogis_emstep_mo(params, data)
      eres <- em_tlogis_estep(params, data)
      res <- optim(par=params[2:3], fn=em_tlogis_pllf, data=data,
        w0=eres$w0, w1=eres$w1, control=list(fnscale=-1))
      list(
        llf=eres$llf,
        param=c(eres$omega, res$par),
        pdiff=c(eres$omega, res$par) - params,
        total=eres$total
      )
    }
  )
)

#' @rdname NHPP
#' @export
LLogisSRM <- R6::R6Class("LLogisSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      pllogis(t, locationlog=self$locationlog(), scalelog=self$scalelog(),
        lower.tail=lower.tail)
    },
    invFt = function(p) {
      qllogis(p, locationlog=self$locationlog(), scalelog=self$scalelog())
    },
    ft = function(t) {
      dllogis(t, locationlog=self$locationlog(), scalelog=self$scalelog())
    },
    rf = function(n) {
      rllogis(n, locationlog=self$locationlog(), scalelog=self$scalelog())
    }
  ),
  public = list(
    name = "llogis",
    df = 3,
    locationlog = function() { self$params[2L] },
    scalelog = function() { self$params[3L] },
    initialize = function(omega = 1, locationlog = 0, scalelog = 1) {
      self$params <- c(omega, locationlog, scalelog)
      names(self$params) <- c("omega", "locationlog", "scalelog")
    },
    init_params = function(data) {
      self$params <- c(1.0, 1.0, max(log(data$mean), 1.0))
    },
    em = function(params, data) {
      eres <- em_llogis_estep(params, data)
      res <- optim(par=params[2:3], fn=em_llogis_pllf, data=data, w1=eres$w1,
        control=list(fnscale=-1))
      list(
        llf=eres$llf,
        param=c(eres$omega, res$par),
        pdiff=c(eres$omega, res$par) - params,
        total=eres$total
      )
    }
  )
)

#' @rdname NHPP
#' @export
TXVMaxSRM <- R6::R6Class("TXVMaxSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      ptgumbel(t, loc=self$loc(), scale=self$scale(), lower.tail=lower.tail)
    },
    invFt = function(p) { qtgumbel(p, loc=self$loc(), scale=self$scale()) },
    ft = function(t) { dtgumbel(t, loc=self$loc(), scale=self$scale()) },
    rf = function(n) { rtgumbel(n, loc=self$loc(), scale=self$scale()) }
  ),
  public = list(
    name = "txvmax",
    df = 3,
    loc = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, loc = 0, scale = 1) {
      self$params <- c(omega, loc, scale)
      names(self$params) <- c("omega", "loc", "scale")
    },
    init_params = function(data) {
      self$params <- c(1.0, 0.0, data$max/3)
    },
    em = function(params, data) {
      # em_txvmax_emstep_mo(params, data)
      eres <- em_txvmax_estep(params, data)
      res <- optim(par=params[2:3], fn=em_txvmax_pllf, data=data,
        w0=eres$w0, w1=eres$w1, control=list(fnscale=-1))
      list(
        llf=eres$llf,
        param=c(eres$omega, res$par),
        pdiff=c(eres$omega, res$par) - params,
        total=eres$total
      )
    }
  )
)

#' @rdname NHPP
#' @export
LXVMaxSRM <- R6::R6Class("LXVMaxSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      plgumbel(t, loclog=self$loclog(), scalelog=self$scalelog(),
        lower.tail=lower.tail)
    },
    invFt = function(p) {
      qlgumbel(p, loclog=self$loclog(), scalelog=self$scalelog())
    },
    ft = function(t) {
      dlgumbel(t, loclog=self$loclog(), scalelog=self$scalelog())
    },
    rf = function(n) {
      rlgumbel(n, loclog=self$loclog(), scalelog=self$scalelog())
    }
  ),
  public = list(
    name = "lxvmax",
    df = 3,
    loclog = function() { self$params[2L] },
    scalelog = function() { self$params[3L] },
    initialize = function(omega = 1, loclog = 0, scalelog = 1) {
      self$params <- c(omega, loclog, scalelog)
      names(self$params) <- c("omega", "loclog", "scalelog")
    },
    init_params = function(data) {
      self$params <- c(1.0, 1.0, max(log(data$max), 1.0))
    },
    em = function(params, data) {
      eres <- em_lxvmax_estep(params, data)
      res <- optim(par=params[2:3], fn=em_lxvmax_pllf, data=data, w1=eres$w1,
        control=list(fnscale=-1))
      list(
        llf=eres$llf,
        param=c(eres$omega, res$par),
        pdiff=c(eres$omega, res$par) - params,
        total=eres$total
      )
    }
  )
)

#' @rdname NHPP
#' @export
TXVMinSRM <- R6::R6Class("TXVMinSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      ptgumbel.min(t, loc=self$loc(), scale=self$scale(),
        lower.tail=lower.tail)
    },
    invFt = function(p) {
      qtgumbel.min(p, loc=self$loc(), scale=self$scale())
    },
    ft = function(t) {
      dtgumbel.min(t, loc=self$loc(), scale=self$scale())
    },
    rf = function(n) {
      rtgumbel.min(n, loc=self$loc(), scale=self$scale())
    }
  ),
  public = list(
    name = "txvmin",
    df = 3,
    loc = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, loc = 0, scale = 1) {
      self$params <- c(omega, loc, scale)
      names(self$params) <- c("omega", "loc", "scale")
    },
    init_params = function(data) {
      self$params <- c(data$total, -data$mean, data$max/3)
    },
    em = function(params, data) {
      eres <- em_txvmin_estep(params, data)
      res <- optim(par=params[2:3], fn=em_txvmin_pllf, data=data,
        w0=eres$w0, w1=eres$w1, control=list(fnscale=-1))
      list(
        llf=eres$llf,
        param=c(eres$omega, res$par),
        pdiff=c(eres$omega, res$par) - params,
        total=eres$total
      )
    }
  )
)

#' @rdname NHPP
#' @export
LXVMinSRM <- R6::R6Class("LXVMinSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      plgumbel.min(t, loclog=self$loclog(), scalelog=self$scalelog(),
        lower.tail=lower.tail)
    },
    invFt = function(p) {
      qlgumbel.min(p, loclog=self$loclog(), scalelog=self$scalelog())
    },
    ft = function(t) {
      dlgumbel.min(t, loclog=self$loclog(), scalelog=self$scalelog())
    },
    rf = function(n) {
      rlgumbel.min(n, loclog=self$loclog(), scalelog=self$scalelog())
    }
  ),
  public = list(
    name = "lxvmin",
    df = 3,
    loclog = function() { self$params[2L] },
    scalelog = function() { self$params[3L] },
    initialize = function(omega = 1, loclog = 0, scalelog = 1) {
      self$params <- c(omega, loclog, scalelog)
      names(self$params) <- c("omega", "loclog", "scalelog")
    },
    init_params = function(data) {
      self$params <- c(1.0, 0.0, max(log(data$max), 1.0))
    },
    em = function(params, data) {
      eres <- em_lxvmin_estep(params, data)
      res <- optim(par=params[2:3], fn=em_lxvmin_pllf, data=data, w1=eres$w1,
        control=list(fnscale=-1))
      list(
        llf=eres$llf,
        param=c(eres$omega, res$par),
        pdiff=c(eres$omega, res$par) - params,
        total=eres$total
      )
    }
  )
)
