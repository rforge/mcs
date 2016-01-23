\name{lmSelect}


\alias{lmSelect}
\alias{lmSelect.fit}
\alias{lmSubsets.select}

\alias{print.lmSelect}
\alias{plot.lmSelect}


\title{All-Subsets Regression}


\description{Best-subsets regression for ordinary linear models.}

\usage{
lmSelect(formula, data, subset, weights, na.action,
model = TRUE, x = FALSE, y = FALSE, contrasts = NULL, offset, \ldots)

lmSelect.fit(x, y, weights = NULL, offset = NULL,
include = NULL, exclude = NULL, penalty = "BIC", tolerance = 0,
pradius = NULL, nbest = 1, \ldots, .algo = "hpbba")

lmSubsets.select(object, penalty = "BIC", ...)

lmSubsets.select(object, penalty = "BIC", ...)
}


\arguments{
  \item{object}{An object of class \code{lmSubsets}.}
  \item{formula, data, subset, weights, na.action, model, contrasts,
    offset}{Standard formula interface.}
  \item{x, y}{The model matrix and response.}
  \item{include, exclude}{Index vectors designating variables that are
    forced in or out of the model, respectively.}
  \item{penalty}{Penalty per parameter (see \code{\link[stats]{AIC}}).}
  \item{tolerance}{Heuristic tolerance.}
  \item{pradius}{Preordering radius.}
  \item{nbest}{Number of best subsets.}
  \item{\dots}{Ignored.}
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
  \item{intercept}{\code{TRUE} if model has intercept term;
    \code{FALSE} otherwise.}
  \item{include}{Included variables.}
  \item{exclude}{Excluded variables.}
  \item{size}{Subset sizes.}
  \item{penalty}{Penalty per parameter.}
  \item{tolerance}{Heuristic tolerance.}
  \item{nbest}{Number of best subsets.}
  \item{df}{Degrees of freedom.}
  \item{rss}{Residual sum of squares.}
  \item{val}{Subset values.}
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
lmSel.AP <- lmSelect(mortality ~ ., data = AirPollution, nbest = 20)

## equivalent to:
\dontrun{
lmSub.AP <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 20)
lmSel.AP <- lmSelect(lmSub.AP)
}

## visualize RSS
plot(lmSel.AP)

## summarize
summary(lmSel.AP)
}

\keyword{regression}
