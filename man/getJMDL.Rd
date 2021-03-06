% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilite.R
\name{getJMDL}
\alias{getJMDL}
\alias{getJMDL.jmdlMod}
\title{Extract or Get Generalized Components from a Fitted Joint Mean
Correlation Model}
\usage{
getJMDL(object, name)

\method{getJMDL}{jmdlMod}(object, name = c("m", "Y", "X", "W", "offset",
  "theta", "beta", "gamma", "delta", "loglik", "family", "q", "time", "stdbeta",
  "stdgamma", "stddelta"))
}
\arguments{
\item{object}{a fitted joint mean correlation model of class "JmdlMod", i.e.,
typically the result of jmdl().}

\item{name}{a character vector specifying the name(s) of the "component".

possible values are:
\describe{
  \item{\code{"m"}}{a vector of number of measurement for each subject}
  \item{\code{"Y"}}{response matrix}
  \item{\code{"X"}}{model matrix for mean structure}
  \item{\code{"W"}}{model matrix for correlation structure (the lower
  triangular matrix)}
  \item{\code{"offset"}}{a vecter to be added to a linear predictor}
  \item{\code{"theta"}}{parameter estimates of joint mean correlation model}
  \item{\code{"beta"}}{parameter estimates for mean structure model}
  \item{\code{"delta"}}{parameter estimates for mean structure model (for
  Nbinom model)}
  \item{\code{"gamma"}}{parameter estimates for correlation structure (the
  lower triangular matrix)}
  \item{\code{"stdbeta"}}{standard error for parameter beta}
  \item{\code{"stddelta"}}{standard error for parameter delta}
  \item{\code{"stdgamma"}}{standard error for parameter gamma}
  \item{\code{"loglik"}}{log-likelihood, except for a constant}
  \item{\code{"family"}}{the marginal distributions of the discrete variables}
  \item{\code{"q"}}{degree of polynomial of the time lag to model the lower
  triangular matrix}
  \item{\code{"time"}}{a vector of time from the data}
}}
}
\description{
Extract (or "get") "components" - in a generalized sense - from
a fitted joint mean correlation model from an object of class "JmdlMod".
}
\section{Methods (by class)}{
\itemize{
\item \code{jmdlMod}: Extract or Get Generalized Components from a Fitted Joint
Mean Correlation Model
}}

\examples{

mydat <- toydata
fit <- jmdl(Y|id|time ~ X, data = mydat, q = 2, family ='Bernoulli')
beta <- getJMDL(fit, "beta")
beta
loglik  <- getJMDL(fit, "loglik")
loglik

}
