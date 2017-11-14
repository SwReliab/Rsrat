#' Model names
#'
#' All software reliablity models.
#'
#' Return a list of all models.
#' \describe{
#' \item{exp}{Exponential model (Goel-Okumoto model).}
#' \item{gamma}{Gamma model.}
#' \item{pareto}{Pareto model (Littlewood model).}
#' \item{tnorm}{Truncated normal model.}
#' \item{lnorm}{Log-normal model.}
#' \item{tlogis}{Truncated logisitc model (Inflection S-shpaed model, Logistic curve).}
#' \item{txvmax}{Truncated exttreme-value max model (Gompertz curve).}
#' \item{lxvmax}{Log-exttreme-value max model (Frechet typed).}
#' \item{txvmin}{Truncated exttreme-value min model.}
#' \item{lxvmin}{Log-exttreme-value min model (Goel model, Weibull model).}
#' }
#'
#' @return A list of all software reliability models (exp, gamma,
#' pareto, tnorm, lnorm, tlogis, llogis, txvmax, lxvmax, txvmin, lxvmin).
#' @export

srm.models <- c(
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

#' Software reliability model
#'
#' Generate a software reliability model.
#'
#' Return a class of NHPP.
#'
#' @param names A character vector, indicating the model (\code{\link{srm.models}}).
#' @return A class of NHPP (\code{\link{NHPP}}).
#' @export

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
    "dglm.logit"=dGLM.logit$new(),
    "dglm.probit"=dGLM.probit$new(),
    "dglm.cloglog"=dGLM.cloglog$new(),
    NA
  )
}
