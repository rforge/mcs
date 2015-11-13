\name{summary}


\alias{summary}

\alias{summary.lmSubsets}
\alias{print.summary.lmSubsets}
\alias{plot.summary.lmSubsets}

%% \alias{summary.lmSelect}
%% \alias{print.summary.lmSelect}
%% \alias{plot.summary.lmSelect}


\title{Summarize 'lmSubsets' and 'lmSelect' objects}


\description{
  Summary methods for objects of class \code{lmSubsets} and
  \code{lmSelect}.
}


\usage{
\method{summary}{lmSubsets}(object, penalty = 2, \dots)
%% \method{summary}{lmSelect}(object, penalty = 2, \dots)
}


\arguments{
  \item{object}{An object of class \code{lmSubsets} or \code{lmSelect}.}
  \item{penalty}{The AIC penalty.}
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

## fit best subsets
xs <- lmSubsets(mortality ~ ., data = AirPollution)

## summarize (BIC)
sx <- summary(xs, penalty = log(nrow(AirPollution)))
sx

## or
sx <- summary(xs, penalty = "BIC")
sx
}

\keyword{regression}
