\name{methods}


\alias{methods}

\alias{variable.names.lmSubsets}
\alias{formula.lmSubsets}
\alias{model.frame.lmSubsets}
\alias{model.matrix.lmSubsets}
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
\alias{formula.lmSelect}
\alias{model.frame.lmSelect}
\alias{model.matrix.lmSelect}
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
\method{variable.names}{lmSubsets}(object, size, best = 1, \ldots,
.full = FALSE, .cmpl = FALSE)
\method{formula}{lmSubsets}(x, \ldots)
\method{model.frame}{lmSubsets}(formula, \ldots)
\method{model.matrix}{lmSubsets}(object, size, best = 1, \ldots)
\method{refit}{lmSubsets}(object, \ldots)
\method{coef}{lmSubsets}(object, \ldots)
\method{vcov}{lmSubsets}(object, \ldots)
\method{fitted}{lmSubsets}(object, \ldots)
\method{residuals}{lmSubsets}(object, \ldots)
\method{deviance}{lmSubsets}(object, size, best = 1, \ldots)
\method{logLik}{lmSubsets}(object, size, best = 1, \ldots)
\method{AIC}{lmSubsets}(object, size, best = 1, \ldots, k = 2)
\method{BIC}{lmSubsets}(object, size, best = 1, \ldots)

\method{variable.names}{lmSelect}(object, best = 1, \ldots,
.full = FALSE, .cmpl = FALSE)
\method{formula}{lmSelect}(x, \ldots)
\method{model.frame}{lmSelect}(formula, \ldots)
\method{model.matrix}{lmSelect}(object, best = 1, \ldots)
\method{refit}{lmSelect}(object, \ldots)
\method{coef}{lmSelect}(object, \ldots)
\method{vcov}{lmSelect}(object, \ldots)
\method{fitted}{lmSelect}(object, \ldots)
\method{residuals}{lmSelect}(object, \ldots)
\method{deviance}{lmSelect}(object, best = 1, \ldots)
\method{logLik}{lmSelect}(object, best = 1, \ldots)
\method{AIC}{lmSelect}(object, best = 1, \ldots, k = 2)
\method{BIC}{lmSelect}(object, best = 1, \ldots)
}


\arguments{
  \item{object, formula, x}{An object of class \code{lmSubsets}
    or \code{lmSelect}.}
  \item{size}{The subset size.}
  \item{best}{The subset rank.}
  \item{k}{AIC penalty.}
  \item{\dots}{Forwarded arguments.}
  \item{.full, .cmpl}{For internal use.}
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
sub1 <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 5)

## extract information (for subset of size 3)
coef(sub1, size = 3)
vcov(sub1, size = 3)
residuals(sub1, size = 3)
fitted(sub1, size = 3)
model.matrix(sub1, size = 3)

## select best subsets
sel1 <- lmSelect(sub1)

## extract information (for best BIC subset)
deviance(sel1)
logLik(sel1)
AIC(sel1)
BIC(sel1, best = 1:5)

## refit model
lm5 <- refit(sub1, size=5)
summary(lm5)
## (Note that the p-values are not valid due to model selection.)
}


\keyword{regression}
