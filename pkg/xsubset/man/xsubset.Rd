\name{xsubset}

\alias{xsubset}
\alias{xsubset.model}
\alias{xsubset.print}

\title{Variable subset selection}

\description{Compute the best variable subsets.}

\usage{
  xsubset(a, y, icpt = TRUE, prad = 1, tau = 0, vmin = NULL, vmax = NULL)
  xsubset.model(x, size)
  xsubset.print(x)
}

\arguments{
  \item{a}{Data matrix (independent variables).}
  \item{y}{Dependent variable.}
  \item{icpt}{Intercept.}
  \item{prad}{Preordering radius.}
  \item{tau}{Tolerance parameter.}
  \item{vmin}{Minimal subset size (number of variables).}
  \item{vmax}{Maximal subset size (number of variables).}
  \item{size}{Subset model size.}
  \item{x}{Value returned by \code{xsubset}.}
}

\details{

  The function \code{xsubset} computes the best variable-subset models
  of the model specified by \code{a} and \code{y}, where \code{a} is a
  \eqn{m} by \eqn{n} (\eqn{m>n}) matrix and \code{y} an \eqn{m}-element
  vector.  The best subset model is derived for each number of
  variables.  The residual sum of squares (RSS) and the variable subset
  (indices) is returned for each solution subset model.

  By default, an intercept term is added to the model.  It can be
  disabled by setting \code{icpt=FALSE}.

  The function will preorder the variables to reduce execution time if
  \code{prad}\eqn{>0} (\code{prad}\eqn{\le n}).  Best execution times
  are usually attained for approximately \code{prad}\eqn{=n/3}.

  If the exact solutions are not needed, and to prevent the function
  from performing an exhaustive search, a tolerance \code{tau}\eqn{>0}
  can be specified.  If \code{tau} is a vector of length \eqn{n}, then
  \code{tau[i]} is the tolerance for subset models with \eqn{i}
  variables.

  If only subsets of a certain size are of interest, \code{vmin} and
  \code{vmax} can be used to reduce the search space
  (\eqn{1\le}\code{vmin}\eqn{\le}\code{vmax}\eqn{\le n}).

  The function \code{xsubset.model} extracts the subset model (RSS,
  variable names) for the specified subset size.

  The function \code{xsubset.print} prints all derived subset models.

}

\value{

  The return value \code{x} is a list with two entries:
  \item{\code{rsel}}{A list such that \code{rsel[[i]]}\eqn{=}RSS\eqn{(Si)}.}
  \item{\code{isel}}{A list such that \code{isel[[i]]}\eqn{=Si}.}

  Here, \eqn{Si} denotes the solution subset with \eqn{i} variables
  (\code{vmin}\eqn{\le i\le}\code{vmax}).  Note that \code{x$isel[[i]]}
  is an integer vector that can be used to index the original matrix
  \code{a} to extract the solution subset model.

}

\examples{

  ## load pollute dataset
  data(pollute)

  ## index of dependent variable
  idep <- 16

  ## compute all subset models (with intercept)
  x <- xsubset(pollute[-idep], pollute[idep])
  xsubset.print(x)

  ## double check RSS of solution subset 5
  size <- 5
  ## extract subset model
  s <- xsubset.model(x, size)
  print(s$xnames)
  print(s$yname)
  print(s$rss)
  ## data matrix
  xy <- pollute[c(s$xnames, s$yname)]
  ## add intercept
  ixy <- cbind(1, xy)
  ## compute QR factorization
  irz <- qr.R(qr(ixy))
  ## compute RSS
  rss <- irz[size + 2, size + 2]^2
  print(rss)

  ## shortcut
  xy <- pollute[c(x$isel[[size]], idep)]
  irz <- qr.R(qr(cbind(1, xy)))
  rss <- irz[size + 2, size + 2]^2
  print(rss)

}

\references{

  Hofmann, M. and Gatu, C. and Kontoghiorghes, E. J. (2007).  Efficient
  algorithms for computing the best subset regression models for
  large-scale problems.  Computational Statistics and Data Analysis, 52,
  16-29.

  Gatu, C. and Kontoghiorghes, E. J.  (2006).  Branch-and-bound
  algorithms for computing the best subset regression models.  Journal
  of Computational and Graphical Statistics, 15, 139-156.

}
