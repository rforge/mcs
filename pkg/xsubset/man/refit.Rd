\name{refit}
\alias{refit}

\title{Re-fitting Models}

\description{
Generic function for re-fitting a models on various subsets or
re-weighted data sets.
}

\usage{
  refit(object, \dots)
}

\arguments{
  \item{object}{an object.}
  \item{\dots}{arguments passed to methods.}
}

\details{
  \code{refit} is a new generic function for re-fitting a certain
  model object on multiple versions of a data set (and is hence
  different from \code{update}). Applications are re-fitting models
  after some kind of model selection, e.g., variable subset selection,
  partitioning, re-weighting etc.
    
  The generic is similar to the one provided in \pkg{modeltools}
  and \pkg{fxregime}
  (and should fulfill the same purpose). To avoid dependence on
  this package, it is also provided here.
}

\seealso{\code{\link[xsubset]{xsubset}}}

\keyword{regression}
