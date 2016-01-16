\name{summary}


\alias{summary}

\alias{summary.lmSubsets}
\alias{print.summary.lmSubsets}
\alias{plot.summary.lmSubsets}

\alias{summary.lmSelect}
\alias{print.summary.lmSelect}
\alias{plot.summary.lmSelect}


\title{Summarize 'lmSubsets' and 'lmSelect' objects}


\description{
  Summary methods for objects of class \code{lmSubsets} and
  \code{lmSelect}.
}


\usage{
\method{summary}{lmSubsets}(object, penalty = "BIC", \dots)
\method{summary}{lmSelect}(object, penalty = "AIC", \dots)
}


\arguments{
  \item{object}{An object of class \code{lmSubsets} or \code{lmSelect}.}
  \item{penalty}{Penalty per parameter.}
  \item{\dots}{Ignored.}
}


\details{
  The \code{summary} method for \code{lmSubsets} and \code{lmSelect}
  objects computes the AIC given the desired \code{penalty}.
}


\value{
  An object of class \code{summary.lmSubsets} or
  \code{summary.lmSelect}.
}


\seealso{\code{\link{lmSubsets}}, \link{methods}.}


\examples{
## load data
data("AirPollution", package = "lmSubsets")

## fit all subsets
lmSub.AP <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 10)

## summarize (AIC)
summary(lmSub.AP, penalty = "AIC")

## fit best subsets (BIC)
lmSel.AP <- lmSelect(mortality ~ ., data = AirPollution, nbest = 10)

## summarize
summary(lmSel.AP, penalty = list("AIC", 3, 4, 5))
}

\keyword{regression}
