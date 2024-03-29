% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emfit.R
\name{emfit}
\alias{emfit}
\title{Estimate SRM parameters.}
\usage{
emfit(
  srm,
  data,
  initialize = TRUE,
  maxiter = 2000,
  reltol = 1e-06,
  abstol = 0.001,
  trace = FALSE,
  printsteps = 50,
  ...
)
}
\arguments{
\item{srm}{A model.}

\item{data}{A faultdata.}

\item{initialize}{Either TRUE or FALSE. If TRUE, the model parameters are
initilized with a given data before executing the fitting algorithm.}

\item{maxiter}{An integer for the maximum number of iterations in the fitting algorithm.}

\item{reltol}{A numeric value. The algorithm stops if the relative error is
less than \emph{reltol} and the absolute error is less than \emph{abstol}.}

\item{abstol}{A numeric value. The algorithm stops if the relative error is
less than \emph{reltol} and the absolute error is less than \emph{abstol}.}

\item{trace}{A logical. If TRUE, the intermediate parameters are printed.}

\item{printsteps}{An integer for print.}

\item{...}{A list for other parameters which are sent to the \code{em} method of srm.}
}
\value{
A list with components;
\item{initial}{A vector for initial parameters.}
\item{srm}{A class of NHPP. The SRM with the estiamted parameters.}
\item{llf}{A numeric value for the maximum log-likelihood function.}
\item{df}{An integer for degrees of freedom.}
\item{convergence}{A boolean meaning the alorigthm is converged or not.}
\item{iter}{An integer for the number of iterations.}
\item{aerror}{A numeric value for absolute error.}
\item{rerror}{A numeric value for relative error.}
}
\description{
This function provides the maximum likelihood (ML) estiamtes for a given SRM 
with a given data.bThe ML estimates are computed with the EM algorithm.
The initial parameters for the EM algorithm are automatically decided 
if the flag \emph{initialize} is TRUE.
}
\examples{
data(dacs)
data <- faultdata(fault=sys1g)
emfit(srm("exp"), data)
}
