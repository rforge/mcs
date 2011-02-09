\name{summary}

\alias{summary}

\alias{summary.xsubset}
\alias{print.summary.xsubset}
\alias{plot.summary.xsubset}
\alias{refit.summary.xsubset}
\alias{summary.summary.xsubset}

\alias{summary.xsubset2}
\alias{print.summary.xsubset2}
\alias{plot.summary.xsubset2}
\alias{refit.summary.xsubset2}
\alias{summary.summary.xsubset2}


\title{Summarize xsubset and xselect objects}

\description{
  Summary methods for objects of class \code{"xsubset"} and
  \code{"xselect"}.
}

\usage{
\method{summary}{xsubset}(object, rank = 1, penalty = 2, \dots)
\method{summary}{xsubset2}(object, \ldots)

\method{refit}{summary.xsubset}(object, \ldots, mask.call = TRUE)
\method{refit}{summary.xsubset2}(object, \ldots, mask.call = TRUE)
}

\arguments{
  \item{object}{An object of class \code{"xsubset"} or \code{"xselect"}.}
  \item{rank}{The subset rank.}
  \item{penalty}{The AIC penalty.}
  \item{mask.call}{Rebuild call to \code{summary} for pretty printing.}
  \item{\dots}{Ignored.}
}

\details{
  The \code{summary} method for \code{xsubset} objects computes the AIC
  for the desired \code{rank}.  The desired penalty is specified by
  the \code{penalty} argument.

  The \code{refit} may be called on an object of class
  \code{"summary.xsubset"} or \code{"summary.xsubset2"} to obtain an
  \code{"lm"} object that corresponds to the optimal subset size.
}

\value{
  An object of class \code{"summary.xsubset"} or \code{"summary.xsubset2"}.

  Objects of class \code{"summary.xsubset"} are objects of class
  \code{"xsubset"} extended with the following components:
  \item{sum.rank}{Summary rank.}
  \item{sum.penalty}{Summary AIC penalty.}
  \item{sum.rss}{Summary RSS.}
  \item{sum.aic}{Summary AIC.}
  \item{sum.best}{Index of best subset model.}

  Objects of class \code{"summary.xsubset2"} are identical to objects of
  class \code{"xsubset2"}.

  The \code{refit} method returns an object of class \code{"lm"} fitted
  to the best subset size.
}

\seealso{\code{\link{xsubset}}, \link{methods}, \code{\link{refit}}.}

\examples{
## load data (with logs for relative potentials)
data("AirPollution", package = "xsubset")
for(i in 12:14) AirPollution[[i]] <- log(AirPollution[[i]])

## fit best subsets
xs <- xsubset(mortality ~ ., data = AirPollution)

## visualize/summarize RSS and BIC
sx <- summary(xs, penalty = log(nrow(AirPollution)))
print(sx)
plot(sx)

## refit model (best BIC)
lm <- refit(sx)
summary(lm)
}

\keyword{regression}
