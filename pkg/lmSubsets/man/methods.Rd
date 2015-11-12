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


\description{Extractor methods for \code{\link{mcsSubset}} and
\code{\link{lmSelect}} objects.}


\usage{
\method{variable.names}{lmSubsets}(object, size = NULL, best = 1,
\ldots, .full = FALSE, .neg = FALSE)
\method{formula}{lmSubsets}(x, \ldots)
\method{model.frame}{lmSubsets}(formula, \ldots)
\method{model.matrix}{lmSubsets}(object, size = NULL, best = 1, \ldots)
\method{refit}{lmSubsets}(object, \ldots)
\method{coef}{lmSubsets}(object, \ldots)
\method{vcov}{lmSubsets}(object, \ldots)
\method{fitted}{lmSubsets}(object, \ldots)
\method{residuals}{lmSubsets}(object, \ldots)
\method{deviance}{lmSubsets}(object, size = NULL, best = 1, \ldots)
\method{logLik}{lmSubsets}(object, size = NULL, best = 1, \ldots, df)
\method{AIC}{lmSubsets}(object, size = NULL, best = 1, \ldots, k = NULL)
\method{BIC}{lmSubsets}(object, size = NULL, best = 1, \ldots)

\method{variable.names}{lmSelect}(object, size = NULL, best = 1,
\ldots, .full = FALSE, .neg = FALSE)
\method{formula}{lmSelect}(x, \ldots)
\method{model.frame}{lmSelect}(formula, \ldots)
\method{model.matrix}{lmSelect}(object, size = NULL, best = 1, \ldots)
\method{refit}{lmSelect}(object, \ldots)
\method{coef}{lmSelect}(object, \ldots)
\method{vcov}{lmSelect}(object, \ldots)
\method{fitted}{lmSelect}(object, \ldots)
\method{residuals}{lmSelect}(object, \ldots)
\method{deviance}{lmSelect}(object, size = NULL, best = 1, \ldots)
\method{logLik}{lmSelect}(object, size = NULL, best = 1, \ldots, df)
\method{AIC}{lmSelect}(object, size = NULL, best = 1, \ldots, k = NULL)
\method{BIC}{lmSelect}(object, size = NULL, best = 1, \ldots)
}


\arguments{
  \item{object, formula, x}{An object of class \code{mcsSubset}.}
  \item{size}{The subset size.}
  \item{best}{The subset rank.}
  \item{df}{Degrees of freedom.}
  \item{k}{AIC penalty.}
  \item{\dots}{Forwarded arguments.}
  \item{.full}{For internal use.}
  \item{.neg}{For internal use.}
}


\details{
  The extractor methods work for \code{lmSubsets} and \code{lmSelect}
  objects that have been generated using the formula interface.

  The information is extracted from the model refitted to a given size.
  If a method is not available, \code{refit} can always be called
  explicitly to obtain a new \code{lm} object of the desired size or
  rank.  In some cases, the \code{size} and \code{best} arguments may be
  a vector specifying multiple values.

  The method \code{refit} returns an \code{lm} object fitted to the
  desired size and rank.
}


\seealso{\code{\link{lmSubsets}}, \code{\link{lmSelect}}, \code{\link{refit}}.}


\examples{
## load data
data("AirPollution", package = "lmSubsets")

## fit subsets (5 best subsets per size)
xs <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 5)

## extract information (for subset of size 3)
coef(xs, size = 3)
vcov(xs, size = 3)
residuals(xs, size = 3)
fitted(xs, size = 3)
model.matrix(xs, size = 3)

## summarize (BIC)
sx <- summary(xs, penalty = log(nrow(AirPollution)))

## extract information (for best BIC subset)
deviance(sx)
logLik(sx)
AIC(sx)

## refit model
lm5 <- refit(sx, size=5)
summary(lm5)
## (Note that the p-values are not valid due to model selection.)

## select 5 best subsets using AIC
xs <- lmSelect(mortality ~ ., data = AirPollution, penalty = 2, nbest = 5)

## summarize
summary(xs)

## extract deviance for best subset
deviance(xs)

## extract BIC for all subsets
BIC(xs, best = 1:5)
}


\keyword{regression}
