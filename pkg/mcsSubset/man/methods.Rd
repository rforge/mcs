\name{methods}

\alias{methods}

\alias{variable.names.mcsSubset}
\alias{formula.mcsSubset}
\alias{model.frame.mcsSubset}
\alias{model.matrix.mcsSubset}
\alias{refit.mcsSubset}
\alias{coef.mcsSubset}
\alias{vcov.mcsSubset}
\alias{fitted.mcsSubset}
\alias{residuals.mcsSubset}
\alias{deviance.mcsSubset}
\alias{logLik.mcsSubset}
\alias{AIC.mcsSubset}


\title{Methods for 'mcsSubset' objects}


\description{Extractor methods for \code{\link{mcsSubset}} objects.}


\usage{
\method{variable.names}{mcsSubset}(object, size = NULL, best = 1,
\ldots, .full = FALSE, .neg = FALSE)
\method{formula}{mcsSubset}(x, \ldots)
\method{model.frame}{mcsSubset}(formula, \ldots)
\method{model.matrix}{mcsSubset}(object, size = NULL, best = 1, \ldots)
\method{refit}{mcsSubset}(object, \ldots)
\method{coef}{mcsSubset}(object, \ldots)
\method{vcov}{mcsSubset}(object, \ldots)
\method{fitted}{mcsSubset}(object, \ldots)
\method{residuals}{mcsSubset}(object, \ldots)
\method{deviance}{mcsSubset}(object, size = NULL, best = 1, \ldots)
\method{logLik}{mcsSubset}(object, size = NULL, best = 1, \ldots, df)
\method{AIC}{mcsSubset}(object, size = NULL, best = 1, \ldots, k = NULL)
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
  The extractor methods work for \code{mcsSubset} objects that
  have been generated using the formula interface.

  The information is extracted from the model refitted to a given size.
  If a method is not available, \code{refit} can always be called
  explicitly to obtain a new \code{lm} object of the desired size or
  rank.  In some cases, the \code{size} and \code{best} arguments may be
  a vector specifying multiple values.

  The method \code{refit} returns an \code{lm} object fitted to the
  desired size and rank.
}


\seealso{\code{\link{mcsSubset}}, \code{\link{refit}}.}


\examples{
## load data
data("AirPollution", package = "mcsSubset")

## fit subsets (5 best subsets per size)
xs <- mcsSubset(mortality ~ ., data = AirPollution, nbest = 5)

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
xs <- mcsSubset(mortality ~ ., data = AirPollution, penalty = 2, nbest = 5)

## summarize
summary(xs)

## extract deviance for best subset
deviance(xs)

## extract BIC for all subsets
AIC(xs, best = 1:5, k = log(nrow(AirPollution)))
}


\keyword{regression}
