\name{methods}

\alias{model.response}
\alias{model.weights}
\alias{model.offset}
\alias{environment}

\alias{environment.xsubset}
\alias{variable.names.xsubset}
\alias{formula.xsubset}
\alias{model.frame.xsubset}
\alias{model.matrix.xsubset}
\alias{model.response.xsubset}
\alias{model.weights.xsubset}
\alias{model.offset.xsubset}
\alias{terms.xsubset}
\alias{refit.xsubset}
\alias{coef.xsubset}
\alias{vcov.xsubset}
\alias{fitted.xsubset}
\alias{residuals.xsubset}
\alias{deviance.xsubset}
\alias{logLik.xsubset}
\alias{AIC.xsubset}


\title{Methods for xsubset objects}

\description{

  Methods \code{\link{xsubset}} objects.

}

\usage{

\method{refit}{xsubset}(object, size, \ldots, mask.call = TRUE)

\method{refit}{xselect}(object, \ldots, mask.call = TRUE)


model.response(object, \ldots)

model.weights(object, \ldots)

model.offset(object, \ldots)

environment(object, \ldots)


\ldots

}

\arguments{

  \item{object}{An object of class \code{xsubset}, \code{xselect}, or
    \code{lm}.}

  \item{size}{The requested size of the refitted model.}

  \item{mask.call}{Used internally.}

  \item{\dots}{Ignored.}

}

\details{

  The method \code{refit.xsubset} refits the associated \code{lm} object
  to the desired size.  The method \code{refit.xselect} refits to the
  best size computed by \code{xselect} (see the \code{best} component of
  the return value of \code{\link{xselect.default}}).

  Many methods that work for \code{lm} objects also work for
  \code{xsubset} and \code{xselect} objects.  Currently, the methods
  include \code{\link[stats]{variable.names}},
  \code{\link[stats]{formula}}, \code{\link[stats]{model.frame}},
  \code{\link[stats]{model.matrix}}, \code{\link[stats]{terms}},
  \code{\link[stats]{coef}}, \code{\link[stats]{vcov}},
  \code{\link[stats]{fitted}}, \code{\link[stats]{residuals}},
  \code{\link[stats]{deviance}}, \code{\link[stats]{logLik}}, and
  \code{\link[stats]{AIC}}.  Additionally, new generics have been
  defined.  These are \code{model.response}, \code{model.weights},
  \code{model.offset}, and \code{environment}.

  The information is extracted from the model refitted to a given size.
  If a method is not available, \code{refit} can always be called
  explicitly to obtain a new \code{lm} object of the desired size.  Most
  methods for \code{xsubset} objects require a \code{size} argument,
  whereas methods for \code{xselect} objects determine the size from the
  best subset size computed by the \code{xselect} method.  In some
  cases, the \code{size} argument may be a vector specifying multiple
  sizes.

}

\seealso{\code{\link{xsubset}}, \code{\link{xselect}},
  \code{\link{summary}}.}

\examples{

## load data (with logs for relative potentials)
data("AirPollution", package = "xsubset")
for(i in 12:14) AirPollution[[i]] <- log(AirPollution[[i]])


#############
## xsubset ##
#############

## fit best subsets
xs <- xsubset(mortality ~ ., data = AirPollution)

## summarize BIC
sx <- summary(xs, aic.penalty = log(nrow(AirPollution)))

## extract information
deviance(xs, size = sx$best)
logLik(xs, size = sx$best)
AIC(xs, size = sx$best)
AIC(xs, size = sx$best, k = log(nrow(AirPollution)))

## extract information for all subsets fitted
deviance(xs)
AIC(xs)
AIC(xs, k = log(nrow(AirPollution)))

## refit model
lm5 <- refit(xs, size = sx$best)
summary(lm5)
## (Note that the p-values are not valid due to model selection.)

## extract model information
coef(xs, size = 3)
vcov(xs, size = 3)
residuals(xs, size = 3)
fitted(xs, size = 3)
model.matrix(xs, size = 3)


#############
## xselect ##
#############

## select best subset
xs <- xselect(mortality ~ ., data = AirPollution, criterion = "Cp")

## summarize
summary(xs)

## extract informaion for best model
deviance(xs)
logLik(xs)
AIC(xs)
AIC(xs, k = log(nrow(AirPollution)))

}

\keyword{regression}
