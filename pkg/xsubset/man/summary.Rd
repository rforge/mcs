\name{summary}

\alias{summary}

\alias{summary.xsubset}
\alias{print.summary.xsubset}
\alias{plot.summary.xsubset}

\alias{summary.xsubset2}


\title{Summarize xsubset objects}

\description{
  Summary methods for objects of class \code{"xsubset"}.
}

\usage{
\method{summary}{xsubset}(object, penalty = 2, \dots)
\method{summary}{xsubset2}(object, penalty = NULL, \ldots)
}

\arguments{
  \item{object}{An object of class \code{"xsubset"}.}
  \item{penalty}{The AIC penalty.}
  \item{\dots}{Ignored.}
}

\details{
  The \code{summary} method for \code{xsubset} objects computes the AIC
  for the desired \code{penalty}.
}

\value{
  An object of class \code{"summary.xsubset"}.  The subset models
  contained in the summary object are sorted in order of increasing AIC.

  Extractor methods may be called on summary objects to obtain
  information about the best subset models.  See \code{\link{methods}}
  for more details.
}

\seealso{\code{\link{xsubset}}, \link{methods}.}

\examples{
## load data
data("AirPollution", package = "xsubset")

## fit best subsets
xs <- xsubset(mortality ~ ., data = AirPollution)

## summarize (BIC)
sx <- summary(xs, penalty = log(nrow(AirPollution)))
sx

## refit (best BIC)
lm <- refit(sx)
summary(lm)

## refit (second best BIC)
lm <- refit(sx, rank = 2)
summary(lm)
}

\keyword{regression}
