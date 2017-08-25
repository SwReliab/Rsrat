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
