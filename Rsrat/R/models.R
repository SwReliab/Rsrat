# R6 classes

srm.names <- c(
  "exp",
  "gamma",
  "pareto",
  "tnorm",
  "lnorm",
  "tlogis",
  "llogis",
  "txvmax",
  "lxvmax",
  "txvmin",
  "lxvmin"
)

srm <- function(names) {
  if (length(names) == 1L) {
    create.srm.model(names)
  }
  else {
    result <- lapply(names, create.srm.model)
    names(result) <- names
    result
  }
}

create.srm.model <- function(name) {
  switch(name,
    "exp"=ExpSRM$new(),
    "gamma"=GammaSRM$new(),
    "pareto"=ParetoSRM$new(),
    "tnorm"=TNormSRM$new(),
    "lnorm"=LNormSRM$new(),
    "tlogis"=TLogisSRM$new(),
    "llogis"=LLogisSRM$new(),
    "txvmax"=TXVMaxSRM$new(),
    "lxvmax"=LXVMaxSRM$new(),
    "txvmin"=TXVMinSRM$new(),
    "lxvmin"=LXVMinSRM$new(),
    NA
  )
}

NHPP <- R6Class("NHPP",
  private = list(
    Ft = function(t, lower.tail = TRUE) { NA },
    invFt = function(p) { NA },
    ft = function(t) { NA }
  ),
  public = list(
    name = NA,
    params = NA,
    df = NA,
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
    },
    init_params = function(data) { NA },
    set_params = function(params) { self$params <- params },
    em = function(params, data) { NA }
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
    name = "ExpSRM",
    df = 2,
    rate = function() { self$params[2L] },
    initialize = function(omega = 1, rate = 1) {
      self$params <- c(omega, rate)
    },
    init_params = function(data) {
      self$params <- c(data$total, 1.0/data$mean)
    },
    em = function(params, data) {
      em.exp(params, data)
    }
  )
)

GammaSRM <- R6Class("GammaSRM",
  inherit = NHPP,
  private = list(
    Ft = function(t, lower.tail = TRUE) {
      pgamma(t, shape=self$shape(), rate=self$rate(), lower.tail=lower.tail)
    },
    invFt = function(p) { qgamma(p, shape=self$shape(), rate=self$rate()) },
    ft = function(t) { dgamma(t, shape=self$shape(), rate=self$rate()) }
  ),
  public = list(
    name = "GammaSRM",
    df = 3,
    shape = function() { self$params[2L] },
    rate = function() { self$params[3L] },
    initialize = function(omega = 1, shape = 1, rate = 1) {
      self$params <- c(omega, shape, rate)
    },
    init_params = function(data) {
      self$params <- c(data$total, 1.0, 1.0/data$mean)
    },
    em = function(params, data, divide = 15) {
      em.gamma(params, data, divide)
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
      dpareto2(t, shape=self$shape(), scale=self$scale())
    }
  ),
  public = list(
    name = "ParetoSRM",
    df = 3,
    shape = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, shape = 1, scale = 1) {
      self$params <- c(omega, shape, scale)
    },
    init_params = function(data) {
      self$params <- c(1.0, 1.0, 1.0)
    },
    em = function(params, data) {
      em.pareto(params, data)
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
    ft = function(t) { dtnorm(t, mean=self$mean(), sd=self$sd()) }
  ),
  public = list(
    name = "TNormSRM",
    df = 3,
    mean = function() { self$params[2L] },
    sd = function() { self$params[3L] },
    initialize = function(omega = 1, mean = 0, sd = 1) {
      self$params <- c(omega, mean, sd)
    },
    init_params = function(data) {
      self$params <- c(1.0, 0.0, data$mean)
    },
    em = function(params, data) {
      em.tnorm(params, data)
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
      dlnorm(t, meanlog=self$meanlog(), sdlog=self$sdlog())
    }
  ),
  public = list(
    name = "LNormSRM",
    df = 3,
    meanlog = function() { self$params[2L] },
    sdlog = function() { self$params[3L] },
    initialize = function(omega = 1, meanlog = 0, sdlog = 1) {
      self$params <- c(omega, meanlog, sdlog)
    },
    init_params = function(data) {
      self$params <- c(1.0, 1.0, max(log(data$mean), 1.0))
    },
    em = function(params, data) {
      em.lnorm(params, data)
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
      dtlogis(t, location=self$location(), scale=self$scale())
    }
  ),
  public = list(
    name = "TLogisSRM",
    df = 3,
    location = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, location = 0, scale = 1) {
      self$params <- c(omega, location, scale)
    },
    init_params = function(data) {
      self$params <- c(1.0, 0.0, data$mean)
    },
    em = function(params, data) {
      em.tlogist(params, data)
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
      dllogis(t, locationlog=self$locationlog(), scalelog=self$scalelog())
    }
  ),
  public = list(
    name = "LLogisSRM",
    df = 3,
    locationlog = function() { self$params[2L] },
    scalelog = function() { self$params[3L] },
    initialize = function(omega = 1, locationlog = 0, scalelog = 1) {
      self$params <- c(omega, locationlog, scalelog)
    },
    init_params = function(data) {
      self$params <- c(1.0, 1.0, max(log(data$mean), 1.0))
    },
    em = function(params, data) {
      em.llogist(params, data)
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
    ft = function(t) { dtgumbel(t, loc=self$loc(), scale=self$scale()) }
  ),
  public = list(
    name = "TXVMaxSRM",
    df = 3,
    loc = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, loc = 0, scale = 1) {
      self$params <- c(omega, loc, scale)
    },
    init_params = function(data) {
      self$params <- c(1.0, 0.0, data$max/3)
    },
    em = function(params, data) {
      em.txvmax(params, data)
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
      dlgumbel(t, loclog=self$loclog(), scalelog=self$scalelog())
    }
  ),
  public = list(
    name = "LXVMaxSRM",
    df = 3,
    loclog = function() { self$params[2L] },
    scalelog = function() { self$params[3L] },
    initialize = function(omega = 1, loclog = 0, scalelog = 1) {
      self$params <- c(omega, loclog, scalelog)
    },
    init_params = function(data) {
      self$params <- c(1.0, 1.0, max(log(data$max), 1.0))
    },
    em = function(params, data) {
      em.lxvmax(params, data)
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
      dtgumbel.min(t, loc=self$loc(), scale=self$scale())
    }
  ),
  public = list(
    name = "TXVMinSRM",
    df = 3,
    loc = function() { self$params[2L] },
    scale = function() { self$params[3L] },
    initialize = function(omega = 1, loc = 0, scale = 1) {
      self$params <- c(omega, loc, scale)
    },
    init_params = function(data) {
      self$params <- c(data$total, -data$mean, data$max/3)
    },
    em = function(params, data) {
      em.txvmin(params, data)
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
      dlgumbel.min(t, loclog=self$loclog(), scalelog=self$scalelog())
    }
  ),
  public = list(
    name = "LXVMinSRM",
    df = 3,
    loclog = function() { self$params[2L] },
    scalelog = function() { self$params[3L] },
    initialize = function(omega = 1, loclog = 0, scalelog = 1) {
      self$params <- c(omega, loclog, scalelog)
    },
    init_params = function(data) {
      self$params <- c(1.0, 0.0, max(log(data$max), 1.0))
    },
    em = function(params, data) {
      em.lxvmin(params, data)
    }
  )
)
