\name{methods}

\alias{methods}

\alias{variable.names.xsubset}
\alias{formula.xsubset}
\alias{model.frame.xsubset}
\alias{model.matrix.xsubset}
\alias{refit.xsubset}
\alias{coef.xsubset}
\alias{vcov.xsubset}
\alias{fitted.xsubset}
\alias{residuals.xsubset}
\alias{deviance.xsubset}
\alias{logLik.xsubset}
\alias{AIC.xsubset}

\alias{variable.names.xsubset2}
\alias{formula.xsubset2}
\alias{model.frame.xsubset2}
\alias{model.matrix.xsubset2}
\alias{refit.xsubset2}
\alias{coef.xsubset2}
\alias{vcov.xsubset2}
\alias{fitted.xsubset2}
\alias{residuals.xsubset2}
\alias{deviance.xsubset2}
\alias{logLik.xsubset2}
\alias{AIC.xsubset2}


\title{Methods for 'xsubset' objects}


\description{Extractor methods for \code{\link{xsubset}} and
  \code{\link{xsubset2}} objects.}


\usage{
\method{variable.names}{xsubset}(object, size, rank = 1, \ldots)
\method{variable.names}{xsubset2}(object, rank = 1, \ldots)

\method{formula}{xsubset}(x, size, rank = 1, \ldots)
\method{formula}{xsubset2}(x, rank = 1, \ldots)

\method{model.frame}{xsubset}(formula, \ldots)
\method{model.frame}{xsubset2}(formula, \ldots)

\method{model.matrix}{xsubset}(object, size, rank = 1, \ldots)
\method{model.matrix}{xsubset2}(object, rank = 1, \ldots)

\method{refit}{xsubset}(object, size, rank = 1, \ldots, mask.call = TRUE)
\method{refit}{xsubset2}(object, rank = 1, \ldots, mask.call = TRUE)

\method{coef}{xsubset}(object, size, rank = 1, \ldots)
\method{coef}{xsubset2}(object, rank = 1, \ldots)

\method{vcov}{xsubset}(object, size, rank = 1, \ldots)
\method{vcov}{xsubset2}(object, rank = 1, \ldots)

\method{fitted}{xsubset}(object, size, rank = 1, \ldots)
\method{fitted}{xsubset2}(object, rank = 1, \ldots)

\method{residuals}{xsubset}(object, size, rank = 1, \ldots)
\method{residuals}{xsubset2}(object, rank = 1, \ldots)

\method{deviance}{xsubset}(object, size = NULL, rank = 1, \ldots)
\method{deviance}{xsubset2}(object, rank = 1, \ldots)

\method{logLik}{xsubset}(object, size = NULL, rank = 1, \ldots)
\method{logLik}{xsubset2}(object, rank = 1, \ldots)

\method{AIC}{xsubset}(object, size = NULL, rank = 1, \ldots, k = 2)
\method{AIC}{xsubset2}(object, rank = 1, \ldots, k = NULL)
}


\arguments{
  \item{object, formula, x}{An object of class \code{xsubset} or \code{xsubset2}.}
  \item{size}{The subset size.}
  \item{rank}{The subset rank.}
  \item{mask.call}{Used internally.}
  \item{k}{AIC penalty.}
  \item{\dots}{Forwarded arguments.}
}


\details{
  The extractor methods work for \code{xsubset} and \code{xsubset2}
  objects that have been generated using the formula interface.

  The information is extracted from the model refitted to a given size.
  If a method is not available, \code{refit} can always be called
  explicitly to obtain a new \code{lm} object of the desired size or
  rank.  Most methods for \code{xsubset} objects require a \code{size}
  argument and an optional \code{rank} argument, whereas most methods
  for \code{xsubset2} objects take an optional \code{rank} argument.  In
  some cases, the \code{size} and \code{rank} arguments may be a vector
  specifying multiple values.

  The method \code{refit} returns an \code{lm} object fitted to the
  desired size and rank.
}


\seealso{\code{\link{xsubset}}, \code{\link{refit}}.}


\examples{
## load data
data("AirPollution", package = "xsubset")


#############
## xsubset ##
#############

## fit subsets (5 best subsets per size)
xs <- xsubset(mortality ~ ., data = AirPollution, nbest = 5)

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
## xsubset2 ##
##############

## select 5 best subsets using AIC
xs <- xsubset(mortality ~ ., data = AirPollution, penalty = 2, nbest = 5)

## summarize
summary(xs)

## extract deviance for best subset
deviance(xs)

## extract BIC for all subsets
AIC(xs, rank = 1:5, k = log(nrow(AirPollution)))
}


\keyword{regression}
