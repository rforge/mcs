\name{summary}

\alias{summary}

\alias{summary.mcsSubset}
\alias{print.summary.mcsSubset}
\alias{plot.summary.mcsSubset}


\title{Summarize mcsSubset objects}

\description{
  Summary methods for objects of class \code{mcsSubset}.
}

\usage{
\method{summary}{mcsSubset}(object, penalty = 2, \dots)
}

\arguments{
  \item{object}{An object of class \code{mcsSubset}.}
  \item{penalty}{The AIC penalty.}
  \item{\dots}{Ignored.}
}

\details{
  The \code{summary} method for \code{mcsSubset} objects computes the AIC
  for the desired \code{penalty}.
}

\value{
  An object of class \code{summary.mcsSubset}.  The subset models
  contained in the summary object are sorted in order of increasing AIC.

  Extractor methods may be called on summary objects to obtain
  information about the best subset models.  See \code{\link{methods}}
  for more details.
}

\seealso{\code{\link{mcsSubset}}, \link{methods}.}

\examples{
## load data
data("AirPollution", package = "mcsSubset")

## fit best subsets
xs <- mcsSubset(mortality ~ ., data = AirPollution)

## summarize (BIC)
sx <- summary(xs, penalty = log(nrow(AirPollution)))
sx

## refit (best BIC)
lm <- refit(sx, rank=1)
summary(lm)

## refit (second best BIC)
lm <- refit(sx, rank=2)
summary(lm)
}

\keyword{regression}
