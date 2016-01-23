\name{methods}


\alias{methods}

\alias{model.response}

\alias{variable.names.lmSubsets}
\alias{model.frame.lmSubsets}
\alias{model.response.lmSubsets}
\alias{model.matrix.lmSubsets}
\alias{formula.lmSubsets}
\alias{refit.lmSubsets}
\alias{coef.lmSubsets}
\alias{vcov.lmSubsets}
\alias{fitted.lmSubsets}
\alias{residuals.lmSubsets}
\alias{deviance.lmSubsets}
\alias{logLik.lmSubsets}
\alias{AIC.lmSubsets}
\alias{BIC.lmSubsets}

\alias{variable.names.lmSelect}
\alias{model.frame.lmSelect}
\alias{model.response.lmSelect}
\alias{model.matrix.lmSelect}
\alias{formula.lmSelect}
\alias{refit.lmSelect}
\alias{coef.lmSelect}
\alias{vcov.lmSelect}
\alias{fitted.lmSelect}
\alias{residuals.lmSelect}
\alias{deviance.lmSelect}
\alias{logLik.lmSelect}
\alias{AIC.lmSelect}
\alias{BIC.lmSelect}


\title{Methods for 'lmSubsets' and 'lmSelect' objects}


\description{
  Extractor methods for \code{lmSubsets} and \code{lmSelect}
  objects.
}


\usage{
\method{variable.names}{lmSubsets}(object, size, best = 1, \ldots)
\method{model.frame}{lmSubsets}(formula, \ldots)
\method{model.response}{lmSubsets}(data, \ldots)
\method{model.matrix}{lmSubsets}(object, size, best = 1, \ldots)
\method{formula}{lmSubsets}(x, size, best = 1, \ldots)
\method{refit}{lmSubsets}(object, size, best = 1, \ldots)
\method{coef}{lmSubsets}(object, size, best = 1, \ldots)
\method{vcov}{lmSubsets}(object, size, best = 1, \ldots)
\method{fitted}{lmSubsets}(object, size, best = 1, \ldots)
\method{residuals}{lmSubsets}(object, size, best = 1, \ldots)
\method{deviance}{lmSubsets}(object, size, best = 1, \ldots)
\method{logLik}{lmSubsets}(object, size, best = 1, \ldots)
\method{AIC}{lmSubsets}(object, size, best = 1, \ldots, k = 2)
\method{BIC}{lmSubsets}(object, size, best = 1, \ldots)

\method{variable.names}{lmSelect}(object, best = 1, \ldots)
\method{model.frame}{lmSelect}(formula, \ldots)
\method{model.response}{lmSelect}(data, \ldots)
\method{model.matrix}{lmSelect}(object, best = 1, \ldots)
\method{formula}{lmSelect}(x, best = 1, \ldots)
\method{refit}{lmSelect}(object, best = 1, \ldots)
\method{coef}{lmSelect}(object, best = 1, \ldots)
\method{vcov}{lmSelect}(object, best = 1, \ldots)
\method{fitted}{lmSelect}(object, best = 1, \ldots)
\method{residuals}{lmSelect}(object, best = 1, \ldots)
\method{deviance}{lmSelect}(object, best = 1, \ldots)
\method{logLik}{lmSelect}(object, best = 1, \ldots)
\method{AIC}{lmSelect}(object, best = 1, \ldots, k = 2)
\method{BIC}{lmSelect}(object, best = 1, \ldots)
}


\arguments{
  \item{object, formula, data, x}{An object of class \code{lmSubsets} or
    \code{lmSelect}.}
  \item{size}{The subset size.}
  \item{best}{The subset rank.}
  \item{k}{AIC penalty.}
  \item{\dots}{Forwarded arguments.}
}


\details{
  The extractor methods work for \code{lmSubsets} and \code{lmSelect}
  objects that have been generated using the formula interface.

  The information is extracted from the model refitted to a given size.
  If a method is not available, \code{refit} can always be called
  explicitly to obtain a new \code{lm} object of the desired size or
  rank.

  The method \code{refit} returns an \code{lm} object fitted to the
  desired size and rank.
}


\seealso{\code{\link{lmSubsets}}, \code{\link{lmSelect}},
  \code{\link{refit}}.}


\examples{
## load data
data("AirPollution", package = "lmSubsets")

## fit subsets (5 best subsets per size)
lmSub.AP <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 5)

## extract information (for subset of size 3)
coef(lmSub.AP, size = 3)
vcov(lmSub.AP, size = 3)
residuals(lmSub.AP, size = 3)
fitted(lmSub.AP, size = 3)
model.matrix(lmSub.AP, size = 3)

## select best subsets
lmSel.AP <- lmSelect(lmSub.AP)

## extract information (for best BIC subset)
deviance(lmSel.AP)
logLik(lmSel.AP)
AIC(lmSel.AP)
BIC(lmSel.AP, best = 1:5)

## refit model
lm5 <- refit(lmSub.AP, size=5)
summary(lm5)
## (Note that the p-values are not valid due to model selection.)
}


\keyword{regression}
