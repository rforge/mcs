\name{summary}

\alias{summary.xsubset}
\alias{print.summary.xsubset}
\alias{plot.summary.xsubset}

\alias{summary.xselect}
\alias{print.summary.xselect}

\title{Summarizing xsubset and xselect objects}

\description{

  Summary methods for objects of class \code{"xsubset"} and
  \code{"xselect"}.

}

\usage{

\method{summary}{xsubset}(object, size = NULL, aic.penalty = 2, \dots)

\method{summary}{xselect}(object, \ldots)

}

\arguments{

  \item{object}{An object of class \code{"xsubset"} or
    \code{"xselect"}.}

  \item{size}{Vector of subset sizes, defaulting to all subset sizes
    evaluated in \code{object}.}

  \item{aic.penalty}{The penalty used by the modified AIC.}

  \item{\dots}{Ignored.}

}

\details{

  The \code{summary} method for \code{xsubset} ojbects computes the
  modified AIC for all evaluated subset sizes.  The desired penalty is
  specified by the \code{aic.penalty} argument.

}

\value{

  An object of class \code{"summary.xsubset"} or
  \code{"summary.xselect"}.

  Objects of class \code{"summary.xsubset"} are of type \code{"list"}
  and contain the following components:

  \item{rss}{The residual sum of squares of the evaluated subset
    models.}

  \item{aic}{The modified AIC of the evaluated subset models.}

  \item{best}{The size of the best subset model.}

  \item{lm}{The summary of the associated \code{lm} object.}

}

\examples{

## load data (with logs for relative potentials)
data("AirPollution", package = "xsubset")
for(i in 12:14) AirPollution[[i]] <- log(AirPollution[[i]])

## fit best subsets
xs <- xsubset(mortality ~ ., data = AirPollution)

## visualize/summarize RSS and BIC
plot(xs)
sx <- summary(xs, aic.penalty = log(nrow(AirPollution)))
print(sx)
plot(sx)

## refit model (best BIC by default)
lm5 <- refit(xs, size = sx$best)
summary(lm5)

## using 'xselect'
xs <- xselect(mortality ~ ., data = AirPollution, criterion = "BIC")
summary(xs)

}

\keyword{regression}
