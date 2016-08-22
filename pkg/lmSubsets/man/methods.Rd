\name{methods}


\alias{methods}

\alias{model.response}

\alias{variable.names.lmSubsets}
\alias{formula.lmSubsets}
\alias{model.frame.lmSubsets}
\alias{model.matrix.lmSubsets}
\alias{model.response.lmSubsets}
\alias{refit.lmSubsets}
\alias{deviance.lmSubsets}
\alias{logLik.lmSubsets}
\alias{AIC.lmSubsets}
\alias{BIC.lmSubsets}
\alias{coef.lmSubsets}
\alias{vcov.lmSubsets}
\alias{fitted.lmSubsets}
\alias{residuals.lmSubsets}

\alias{variable.names.lmSelect}
\alias{formula.lmSelect}
\alias{model.frame.lmSelect}
\alias{model.matrix.lmSelect}
\alias{model.response.lmSelect}
\alias{refit.lmSelect}
\alias{deviance.lmSelect}
\alias{logLik.lmSelect}
\alias{AIC.lmSelect}
\alias{BIC.lmSelect}
\alias{coef.lmSelect}
\alias{vcov.lmSelect}
\alias{fitted.lmSelect}
\alias{residuals.lmSelect}


\title{Methods for 'lmSubsets' and 'lmSelect' objects}


\description{
  Extractor methods for \code{lmSubsets} and \code{lmSelect}
  objects.
}


\usage{
\method{variable.names}{lmSubsets}(object, size, best = 1, \ldots)
\method{formula}{lmSubsets}(x, size, best = 1, \ldots)
\method{model.frame}{lmSubsets}(formula, size, best = 1, \ldots)
\method{model.matrix}{lmSubsets}(object, size, best = 1, \ldots)
\method{model.response}{lmSubsets}(data, \ldots)
\method{refit}{lmSubsets}(object, size, best = 1, \ldots)
\method{deviance}{lmSubsets}(object, size, best = 1, \ldots, drop = TRUE)
\method{logLik}{lmSubsets}(object, size, best = 1, \ldots, drop = TRUE)
\method{AIC}{lmSubsets}(object, size, best = 1, \ldots, k = 2, drop = TRUE)
\method{BIC}{lmSubsets}(object, size, best = 1, \ldots, drop = TRUE)
\method{coef}{lmSubsets}(object, size, best = 1, \ldots)
\method{vcov}{lmSubsets}(object, size, best = 1, \ldots)
\method{fitted}{lmSubsets}(object, size, best = 1, \ldots)
\method{residuals}{lmSubsets}(object, size, best = 1, \ldots)

\method{variable.names}{lmSelect}(object, best = 1, \ldots, drop = TRUE)
\method{formula}{lmSelect}(x, best, \ldots)
\method{model.frame}{lmSelect}(formula, best, \ldots)
\method{model.matrix}{lmSelect}(object, best, \ldots)
\method{model.response}{lmSelect}(data, \ldots)
\method{refit}{lmSelect}(object, best = 1, \ldots)
\method{deviance}{lmSelect}(object, best = 1, \ldots)
\method{logLik}{lmSelect}(object, best = 1, \ldots, drop = TRUE)
\method{AIC}{lmSelect}(object, best = 1, \ldots, k = 2, drop = TRUE)
\method{BIC}{lmSelect}(object, best = 1, \ldots, drop = TRUE)
\method{coef}{lmSelect}(object, best = 1, \ldots)
\method{vcov}{lmSelect}(object, best = 1, \ldots)
\method{fitted}{lmSelect}(object, best = 1, \ldots)
\method{residuals}{lmSelect}(object, best = 1, \ldots)
}


\arguments{
  \item{object, formula, data, x}{An object of class \code{lmSubsets} or
    \code{lmSelect}.}
  \item{size}{The subset size.}
  \item{best}{The subset rank.}
  \item{\dots}{Forwarded arguments.}
  \item{k}{AIC penalty.}
  \item{drop}{Reduce dimensionality of returned object.}
}


\details{
  The extractor methods work for \code{lmSubsets} and \code{lmSelect}
  objects that have been generated using the formula interface.

  The information is extracted from the model refitted to a given size.
  If a method is not available, \code{refit} can be called explicitly to
  obtain the \code{lm} object.

  For convenience, the submodel size can be inferred from the name of an
  information criterion passed as the \code{size} argument.  Currently,
  only \code{"AIC"} and \code{"BIC"} are recognized.

  The method \code{refit} returns an \code{lm} object fitted to the
  desired size and rank.
}


\seealso{\code{\link{lmSubsets}}, \code{\link{lmSelect}},
  \code{\link{refit}}.}


\examples{
## load data
data("AirPollution", package = "lmSubsets")

## fit subsets (5 best subsets per size)
all.AirPoll <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 5)

## extract information (for best subset of size 3)
coef(all.AirPoll, size = 3)
vcov(all.AirPoll, size = 3)
residuals(all.AirPoll, size = 3)
fitted(all.AirPoll, size = 3)
model.matrix(all.AirPoll, size = 3)

## select best subsets
best.AirPoll <- lmSelect(all.AirPoll)

## extract information (for best BIC subset)
deviance(best.AirPoll)
logLik(best.AirPoll)
AIC(best.AirPoll)
BIC(best.AirPoll, best = 1:5)

## refit model (inferred size)
lm5 <- refit(all.AirPoll, size = "BIC")
summary(lm5)
## (Note that the p-values are not valid due to model selection.)
}


\keyword{regression}
