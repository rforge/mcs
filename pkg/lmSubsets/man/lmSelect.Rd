\name{lmSelect}


\alias{lmSelect}
\alias{lmSelect.lmSubsets}
\alias{lmSelect.lm}
\alias{lmSelect.formula}
\alias{lmSelect.default}

\alias{print.lmSelect}
\alias{plot.lmSelect}


\title{All-Subsets Regression}


\description{Best-subsets regression for ordinary linear models.}

\usage{
lmSelect(object, \dots)

\method{lmSelect}{lmSubsets}(object, \dots, penalty = "AIC")

\method{lmSelect}{formula}(formula, \ldots, lm = FALSE)

\method{lmSelect}{lm}(object, \dots, penalty = "AIC")

\method{lmSelect}{default}(object, y, weights = NULL, offset = NULL,
include = NULL, exclude = NULL, penalty = "AIC", tolerance = 0,
pradius = NULL, nbest = 1, \ldots, .algo = "hbba")
}


\arguments{
  \item{formula, object}{An object of class \code{lm}, \code{formula} or
    \code{matrix}.}
  \item{y}{The response variable.}
  \item{weights, offset}{Part of the standard \code{formula} interface.}
  \item{include, exclude}{Index vectors designating variables that are
    forced in or out of the model, respectively.}
  \item{penalty}{Penalty per parameter (see \code{\link[stats]{AIC}}).}
  \item{tolerance}{Heuristic tolerance.}
  \item{pradius}{Preordering radius.}
  \item{nbest}{Number of best subsets.}
  \item{\dots}{Ignored.}
  \item{lm}{If \code{true}, compute \code{lm} component.}
  \item{.algo}{Internal use.}
}


\details{
  The function \code{lmSelect} computes best-variable-subsets regression
  for ordinary linear models: The \code{nbest} best subset models are
  computed according to an information criterion of the AIC family.

  See \code{\link{lmSubsets}} for further information.
}


\value{
  An object of class \code{"lmSelect"}, i.e. a list with the following
  components:
  \item{nobs}{Number of observations.}
  \item{nvar}{Number of variables (not including intercept, if any).}
  \item{weights}{Weights vector.}
  \item{offset}{Offset component.}
  \item{include}{Included variables.}
  \item{exclude}{Excluded variables.}
  \item{size}{Subset sizes.}
  \item{intercept}{\code{TRUE} if regression has an intercept term;
    \code{FALSE} otherwise.}
  \item{penalty}{Penalty per parameter.}
  \item{tolerance}{Heuristic tolerance.}
  \item{nbest}{Number of best subsets.}
  \item{df}{Degrees of freedom.}
  \item{rss}{Residual sum of squares.}
  \item{aic}{Subset values.}
  \item{which}{Selected variables.}
}


\references{
  Hofmann, M. and Gatu, C. and Kontoghiorghes, E. J. (2007).  Efficient
  Algorithms for Computing the Best Subset Regression Models for
  Large-Scale Problems. \emph{Computational Statistics \& Data Analysis},
  \bold{52}, 16--29.

  Gatu, C. and Kontoghiorghes, E. J.  (2006).  Branch-and-Bound
  Algorithms for Computing the Best Subset Regression Models.
  \emph{Journal of Computational and Graphical Statistics},
  \bold{15}, 139--156.
}


\seealso{\code{\link{lmSubsets}}, \code{\link{summary}},
  \link{methods}.}


\examples{
## load data (with logs for relative potentials)
data("AirPollution", package = "lmSubsets")

#################
## basic usage ##
#################

## canonical example: fit best subsets
xz <- lmSelect(mortality ~ ., data = AirPollution, nbest = 20)

## equivalent to:
# xs <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 20)
# xz <- lmSelect(xs)

## visualize RSS
plot(xz)

## summarize
summary(xz)

## plot summary
plot(summary(xz))
}

\keyword{regression}
