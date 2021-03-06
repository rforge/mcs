\name{refit}

\alias{refit}

\title{Refitting Models}

\description{
  Generic function for refitting a models on various subsets
  or reweighted data sets.
}

\usage{
  refit(object, \dots)
}

\arguments{
  \item{object}{An object.}
  \item{\dots}{Forwarded arguments.}
}

\details{
  \code{refit} is a new generic function for refitting a certain model
  object on multiple versions of a data set (and is hence different from
  \code{update}).  Applications refit models after some kind of
  model selection, e.g., variable subset selection, partitioning,
  reweighting etc.

  The generic is similar to the one provided in \pkg{modeltools} and
  \pkg{fxregime} (and should fulfill the same purpose).  To avoid
  dependence on this package, it is also provided here.
}

\seealso{\code{\link{methods}}}

\keyword{regression}
