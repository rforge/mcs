\name{methods}


\alias{methods}

\alias{model_response}
\alias{log_lik}
\alias{aic}
\alias{bic}

\alias{variable.names.lmSubsets}
\alias{formula.lmSubsets}
\alias{model.frame.lmSubsets}
\alias{model.matrix.lmSubsets}
\alias{model_response.lmSubsets}
\alias{refit.lmSubsets}
\alias{deviance.lmSubsets}
\alias{sigma.lmSubsets}
\alias{log_lik.lmSubsets}
\alias{logLik.lmSubsets}
\alias{aic.lmSubsets}
\alias{AIC.lmSubsets}
\alias{bic.lmSubsets}
\alias{BIC.lmSubsets}
\alias{coef.lmSubsets}
\alias{vcov.lmSubsets}
\alias{fitted.lmSubsets}
\alias{residuals.lmSubsets}

\alias{variable.names.lmSelect}
\alias{formula.lmSelect}
\alias{model.frame.lmSelect}
\alias{model.matrix.lmSelect}
\alias{model_response.lmSelect}
\alias{refit.lmSelect}
\alias{deviance.lmSelect}
\alias{sigma.lmSelect}
\alias{log_lik.lmSelect}
\alias{logLik.lmSelect}
\alias{aic.lmSelect}
\alias{AIC.lmSelect}
\alias{bic.lmSelect}
\alias{BIC.lmSelect}
\alias{coef.lmSelect}
\alias{vcov.lmSelect}
\alias{fitted.lmSelect}
\alias{residuals.lmSelect}


\title{Methods for 'lmSubsets' and 'lmSelect' objects}


\description{
  Extractor methods for \code{"lmSubsets"} and \code{"lmSelect"}
  objects.
}


\usage{
\method{variable.names}{lmSubsets}(object, size, best = 1, \ldots)
\method{formula}{lmSubsets}(x, size, best = 1, \ldots)
\method{model.frame}{lmSubsets}(formula, size, best = 1, \ldots)
\method{model.matrix}{lmSubsets}(object, size, best = 1, \ldots)
\method{model_response}{lmSubsets}(data, \ldots)
\method{refit}{lmSubsets}(object, size, best = 1, \ldots)
\method{deviance}{lmSubsets}(object, size, best = 1, \ldots, drop = TRUE)
\method{sigma}{lmSubsets}(object, size, best = 1, \ldots)
\method{log_lik}{lmSubsets}(object, size, best = 1, \ldots, drop = TRUE)
\method{logLik}{lmSubsets}(object, size, best = 1, \ldots)
\method{aic}{lmSubsets}(object, size, best = 1, \ldots, k = 2, drop = TRUE)
\method{AIC}{lmSubsets}(object, size, best = 1, \ldots, k = 2)
\method{bic}{lmSubsets}(object, size, best = 1, \ldots, drop = TRUE)
\method{BIC}{lmSubsets}(object, size, best = 1, \ldots)
\method{coef}{lmSubsets}(object, size, best = 1, \ldots)
\method{vcov}{lmSubsets}(object, size, best = 1, \ldots)
\method{fitted}{lmSubsets}(object, size, best = 1, \ldots)
\method{residuals}{lmSubsets}(object, size, best = 1, \ldots)

\method{variable.names}{lmSelect}(object, best = 1, \ldots)
\method{formula}{lmSelect}(x, best, \ldots)
\method{model.frame}{lmSelect}(formula, best, \ldots)
\method{model.matrix}{lmSelect}(object, best, \ldots)
\method{model_response}{lmSelect}(data, \ldots)
\method{refit}{lmSelect}(object, best = 1, \ldots)
\method{deviance}{lmSelect}(object, best = 1, \ldots)
\method{sigma}{lmSelect}(object, best = 1, \ldots)
\method{log_lik}{lmSelect}(object, best = 1, \ldots)
\method{logLik}{lmSelect}(object, best = 1, \ldots)
\method{aic}{lmSelect}(object, best = 1, \ldots, k = 2)
\method{AIC}{lmSelect}(object, best = 1, \ldots, k = 2)
\method{bic}{lmSelect}(object, best = 1, \ldots)
\method{BIC}{lmSelect}(object, best = 1, \ldots)
\method{coef}{lmSelect}(object, best = 1, \ldots)
\method{vcov}{lmSelect}(object, best = 1, \ldots)
\method{fitted}{lmSelect}(object, best = 1, \ldots)
\method{residuals}{lmSelect}(object, best = 1, \ldots)
}


\arguments{
  \item{object, formula, data, x}{An object of class \code{"lmSubsets"} or
    \code{"lmSelect"}.}
  \item{size}{The submodel size.}
  \item{best}{The submodel ranking.}
  \item{\dots}{Forwarded arguments.}
  \item{k}{AIC penalty.}
  \item{drop}{Reduce dimensionality of returned object.}
}


\details{
  The extractor methods work for \code{"lmSubsets"} and
  \code{"lmSelect"} objects, i.e., that have been generated using the
  formula interface.

  The information is extracted from a subset model of a given size.  If
  a method is not available, \code{refit} can be called explicitly to
  obtain the corresponding \code{"lm"} object.
}


\seealso{\code{\link{lmSubsets}}, \code{\link{lmSelect}},
  \code{\link{refit}}.}


\examples{
## load data
data("AirPollution", package = "lmSubsets")

## fit subsets (5 best subsets per size)
lm_all <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 5)

## extract information (for best submodel of size 3)
coef(lm_all, size = 3)
vcov(lm_all, size = 3)
residuals(lm_all, size = 3)
fitted(lm_all, size = 3)
model.matrix(lm_all, size = 3)

## select best (BIC) submodels
lm_best <- lmSelect(lm_all)

## extract information
deviance(lm_best)
logLik(lm_best)
AIC(lm_best)
BIC(lm_best, best = 1:5)

## refit model
lm5 <- refit(lm_all, size = 5)
summary(lm5)
## (Note that the p-values are not valid due to model selection.)
}


\keyword{regression}
