# R6 classes

NHPP <- R6Class("NHPP",
  private = list(
    Ft = function(t, lower.tail = TRUE) { NA },
    invFt = function(p) { NA },
    ft = function(t) { NA }
  ),
  public = list(
    name = NA,
    params = NA,
    omega = function() { self$params[1L] },
    mvf = function(t) { self$omega() * private$Ft(t) },
    intensity = function(t) { self$omega() * private$ft(t) },
    reliab = function(t, s) { exp(-(self$mvf(t+s) - self$mvf(s))) },
    residual = function(t) { self$omega() * private$Ft(t, lower.tail=FALSE) },
    ffp = function(t) { exp(-self$residual(t)) },
    median = function(s, p = 0.5) {
      if (p > self$ffp(s)) {
        private$invFt((self$mvf(s) - log(p)) / self$omega())
      } else {
        NA
      }
    }
  )
)

ExpSRM <- R6Class("ExpSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      pexp(t, rate=self$rate(), lower.tail=lower.tail)
    },
    invFt = function(p) { qexp(p, rate=self$rate()) },
    ft = function(t) { dexp(t, rate=self$rate()) }
  ),
  public = list(
    rate = function() { self$params[2L] },
    initialize = function(omega = 1, rate = 1) {
      self$params = c(omega, rate)
    },
    set_params = function(omega, rate) { self$params = c(omega, rate) }
  )
)

GammaSRM <- R6Class("GammaSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      pgamma(t, shape=self$shape(), rate=self$rate(), lower.tail=lower.tail)
    },
    invFt = function(p) { qgamma(p, shape=self$shape(), rate=self$rate()) },
    ft = function(t) { dgamma(x=t, shape=self$shape(), rate=self$rate()) }
  ),
  public = list(
    shape = function() { self$params[2L] },
    rate = function() { self$params[3L] },
    initialize = function(omega = 1, shape = 1, rate = 1) {
      self$params = c(omega, shape, rate)
    },
    set_params = function(omega, shape, rate) {
      self$params = c(omega, shape, rate)
    }
  )
)

ParetoSRM <- R6Class("ParetoSRM",
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
      dpareto2(x=t, shape=self$shape(), scale=self$scale())
    }
  ),
  public = list(
    shape = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, shape = 1, scale = 1) {
      self$params = c(omega, shape, scale)
    },
    set_params = function(omega, shape, scale) {
      self$params = c(omega, shape, scale)
    }
  )
)

TNormSRM <- R6Class("TNormSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      ptnorm(t, mean=self$mean(), sd=self$sd(), lower.tail=lower.tail)
    },
    invFt = function(p) { qtnorm(p, mean=self$mean(), sd=self$sd()) },
    ft = function(t) { dtnorm(x=t, mean=self$mean(), sd=self$sd()) }
  ),
  public = list(
    mean = function() { self$params[2L] },
    sd = function() { self$params[3L] },
    initialize = function(omega = 1, mean = 0, sd = 1) {
      self$params = c(omega, mean, sd)
    },
    set_params = function(omega, mean, sd) {
      self$params = c(omega, mean, sd)
    }
  )
)

LNormSRM <- R6Class("LNormSRM",
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
      dlnorm(x=t, meanlog=self$meanlog(), sdlog=self$sdlog())
    }
  ),
  public = list(
    meanlog = function() { self$params[2L] },
    sdlog = function() { self$params[3L] },
    initialize = function(omega = 1, meanlog = 0, sdlog = 1) {
      self$params = c(omega, meanlog, sdlog)
    },
    set_params = function(omega, meanlog, sdlog) {
      self$params = c(omega, meanlog, sdlog)
    }
  )
)

TLogisSRM <- R6Class("TLogisSRM",
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
      dtlogis(x=t, location=self$location(), scale=self$scale())
    }
  ),
  public = list(
    location = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, location = 0, scale = 1) {
      self$params = c(omega, location, scale)
    },
    set_params = function(omega, location, scale) {
      self$params = c(omega, location, scale)
    }
  )
)

LLogisSRM <- R6Class("LLogisSRM",
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
      dllogis(x=t, locationlog=self$locationlog(), scalelog=self$scalelog())
    }
  ),
  public = list(
    locationlog = function() { self$params[2L] },
    scalelog = function() { self$params[3L] },
    initialize = function(omega = 1, locationlog = 0, scalelog = 1) {
      self$params = c(omega, locationlog, scalelog)
    },
    set_params = function(omega, locationlog, scalelog) {
      self$params = c(omega, locationlog, scalelog)
    }
  )
)

TXVMaxSRM <- R6Class("TXVMaxSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      ptgumbel(t, loc=self$loc(), scale=self$scale(), lower.tail=lower.tail)
    },
    invFt = function(p) { qtgumbel(p, loc=self$loc(), scale=self$scale()) },
    ft = function(t) { dtgumbel(x=t, loc=self$loc(), scale=self$scale()) }
  ),
  public = list(
    loc = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, loc = 0, scale = 1) {
      self$params = c(omega, loc, scale)
    },
    set_params = function(omega, loc, scale) {
      self$params = c(omega, loc, scale)
    }
  )
)

LXVMaxSRM <- R6Class("LXVMaxSRM",
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
      dlgumbel(x=t, loclog=self$loclog(), scalelog=self$scalelog())
    }
  ),
  public = list(
    loclog = function() { self$params[2L] },
    scalelog = function() { self$params[3L] },
    initialize = function(omega = 1, loclog = 0, scalelog = 1) {
      self$params = c(omega, loclog, scalelog)
    },
    set_params = function(omega, loclog, scalelog) {
      self$params = c(omega, loclog, scalelog)
    }
  )
)

TXVMinSRM <- R6Class("TXVMinSRM",
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
      dtgumbel.min(x=t, loc=self$loc(), scale=self$scale())
    }
  ),
  public = list(
    loc = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, loc = 0, scale = 1) {
      self$params = c(omega, loc, scale)
    },
    set_params = function(omega, loc, scale) {
      self$params = c(omega, loc, scale)
    }
  )
)

LXVMinSRM <- R6Class("LXVMinSRM",
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
      dlgumbel.min(x=t, loclog=self$loclog(), scalelog=self$scalelog())
    }
  ),
  public = list(
    loclog = function() { self$params[2L] },
    scalelog = function() { self$params[3L] },
    initialize = function(omega = 1, loclog = 0, scalelog = 1) {
      self$params = c(omega, loclog, scalelog)
    },
    set_params = function(omega, loclog, scalelog) {
      self$params = c(omega, loclog, scalelog)
    }
  )
)
