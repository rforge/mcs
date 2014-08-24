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

\alias{variable.names.mcsSubset2}
\alias{formula.mcsSubset2}
\alias{model.frame.mcsSubset2}
\alias{model.matrix.mcsSubset2}
\alias{refit.mcsSubset2}
\alias{coef.mcsSubset2}
\alias{vcov.mcsSubset2}
\alias{fitted.mcsSubset2}
\alias{residuals.mcsSubset2}
\alias{deviance.mcsSubset2}
\alias{logLik.mcsSubset2}
\alias{AIC.mcsSubset2}


\title{Methods for 'mcsSubset' objects}


\description{Extractor methods for \code{\link{mcsSubset}} and
  \code{\link{mcsSubset2}} objects.}


\usage{
\method{variable.names}{mcsSubset}(object, size, rank = 1, \ldots)
\method{variable.names}{mcsSubset2}(object, rank = 1, \ldots)

\method{formula}{mcsSubset}(x, size, rank = 1, \ldots)
\method{formula}{mcsSubset2}(x, rank = 1, \ldots)

\method{model.frame}{mcsSubset}(formula, \ldots)
\method{model.frame}{mcsSubset2}(formula, \ldots)

\method{model.matrix}{mcsSubset}(object, size, rank = 1, \ldots)
\method{model.matrix}{mcsSubset2}(object, rank = 1, \ldots)

\method{refit}{mcsSubset}(object, size, rank = 1, \ldots, mask.call = TRUE)
\method{refit}{mcsSubset2}(object, rank = 1, \ldots, mask.call = TRUE)

\method{coef}{mcsSubset}(object, size, rank = 1, \ldots)
\method{coef}{mcsSubset2}(object, rank = 1, \ldots)

\method{vcov}{mcsSubset}(object, size, rank = 1, \ldots)
\method{vcov}{mcsSubset2}(object, rank = 1, \ldots)

\method{fitted}{mcsSubset}(object, size, rank = 1, \ldots)
\method{fitted}{mcsSubset2}(object, rank = 1, \ldots)

\method{residuals}{mcsSubset}(object, size, rank = 1, \ldots)
\method{residuals}{mcsSubset2}(object, rank = 1, \ldots)

\method{deviance}{mcsSubset}(object, size = NULL, rank = 1, \ldots)
\method{deviance}{mcsSubset2}(object, rank = 1, \ldots)

\method{logLik}{mcsSubset}(object, size = NULL, rank = 1, \ldots)
\method{logLik}{mcsSubset2}(object, rank = 1, \ldots)

\method{AIC}{mcsSubset}(object, size = NULL, rank = 1, \ldots, k = 2)
\method{AIC}{mcsSubset2}(object, rank = 1, \ldots, k = NULL)
}


\arguments{
  \item{object, formula, x}{An object of class \code{mcsSubset} or \code{mcsSubset2}.}
  \item{size}{The subset size.}
  \item{rank}{The subset rank.}
  \item{mask.call}{Used internally.}
  \item{k}{AIC penalty.}
  \item{\dots}{Forwarded arguments.}
}


\details{
  The extractor methods work for \code{mcsSubset} and \code{mcsSubset2}
  objects that have been generated using the formula interface.

  The information is extracted from the model refitted to a given size.
  If a method is not available, \code{refit} can always be called
  explicitly to obtain a new \code{lm} object of the desired size or
  rank.  Most methods for \code{mcsSubset} objects require a \code{size}
  argument and an optional \code{rank} argument, whereas most methods
  for \code{mcsSubset2} objects take an optional \code{rank} argument.  In
  some cases, the \code{size} and \code{rank} arguments may be a vector
  specifying multiple values.

  The method \code{refit} returns an \code{lm} object fitted to the
  desired size and rank.
}


\seealso{\code{\link{mcsSubset}}, \code{\link{refit}}.}


\examples{
## load data
data("AirPollution", package = "mcsSubset")


#############
## mcsSubset ##
#############

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
lm5 <- refit(sx)
summary(lm5)
## (Note that the p-values are not valid due to model selection.)


##############
## mcsSubset2 ##
##############

## select 5 best subsets using AIC
xs <- mcsSubset(mortality ~ ., data = AirPollution, penalty = 2, nbest = 5)

## summarize
summary(xs)

## extract deviance for best subset
deviance(xs)

## extract BIC for all subsets
AIC(xs, rank = 1:5, k = log(nrow(AirPollution)))
}


\keyword{regression}
