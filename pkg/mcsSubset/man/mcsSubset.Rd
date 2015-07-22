\name{mcsSubset}


\alias{mcsSubset}
\alias{mcsSubset.lm}
\alias{mcsSubset.formula}
\alias{mcsSubset.default}

\alias{print.mcsSubset}
\alias{plot.mcsSubset}


\title{All-Subsets Regression}


\description{All-subsets regression for ordinary linear models.}

\usage{
mcsSubset(object, \dots)

\method{mcsSubset}{formula}(formula, \ldots, lm = FALSE)

\method{mcsSubset}{lm}(object, \dots, penalty = 0)

\method{mcsSubset}{default}(object, y, include = NULL, exclude = NULL,
size = NULL, penalty = 0, tolerance = 0, pradius = NULL, nbest = 1,
\ldots, .algo = "hbba")
}


\arguments{
  \item{formula, object}{An object of class \code{lm}, \code{formula} or
    \code{matrix}.}
  \item{y}{The response variable.}
  \item{include, exclude}{Index vectors designating variables that are
    forced in or out of the model, respectively.  The vectors may
    consist of (integer) indexes, (character) names, or (logical) bits
    selecting the desired columns.  The integer indexes correspond to
    the position of the variables in the model matrix; the intercept, if
    any, has index \code{1}.  By default, all variables are included.}
  \item{size}{Vector of subset sizes (not counting the intercept, if any).
    By default, the best subsets are computed for each subset size (as
    determined by way of \code{include} and \code{exclude}).  Ignored if
    \code{penalty != 0}.}
  \item{penalty}{Penalty per parameter (see \code{\link[stats]{AIC}}).
    If \code{penalty == 0}, determine subsets with lowest RSS for each
    subset size; otherwise, determine subset(s) with overall lowest
    AIC.}
  \item{tolerance}{If \code{penalty == 0}, a numeric vector (expanded if
    necessary), where \code{tolerance[n]} is the tolerance employed for
    subsets of size \code{n}; otherwise, a single value indicating the
    overall tolerance.}
  \item{pradius}{Preordering radius.}
  \item{nbest}{Number of best subsets to report.}
  \item{\dots}{Ignored.}
  \item{lm}{If \code{true}, compute \code{lm} component.}
  \item{.algo}{Internal use.}
}


\details{
  The function \code{mcsSubset} computes all variable-subsets regression
  for ordinary linear models.  The function is generic and provides
  various methods to conveniently specify the regressor and response
  variables.  The standard \code{formula} interface (see
  \code{\link[stats]{lm}}) can be used, or the information can be
  extracted from an already fitted \code{lm} object.  The regressor
  matrix and response variable can also be passed in directly.

  By default (i.e. \code{penalty == 0}), the method computes the
  \code{nbest} best subset models for every subset size, where the
  "best" models are the models with the lowest residual sum of squares
  (RSS).  The scope of the search can be limited to certain subset sizes
  by setting \code{size}.  A tolerance vector (expanded if necessary)
  may be specified to speed up the algorithm, where \code{tolerance[n]}
  is the tolerance applied to subset models of size \code{n}.

  Alternatively (\code{penalty > 0}), the overall (i.e. over all sizes)
  \code{nbest} best subset models may be computed according to an
  information criterion of the AIC family.  A single tolerance value may
  be specified to speed up the search.

  By way of \code{include} and \code{exclude}, variables may be forced
  into or out of the regression, respectively.

  The function will preorder the variables to reduce execution time if
  \code{pradius > 0}.  Good execution times are usually attained for
  approximately \code{pradius = n/3} (default value), where \code{n} is
  the number of regressors after evaluation \code{include} and
  \code{exclude}.

  A set of standard extractor functions for fitted model objects is
  available for objects of class \code{"mcsSubset"}.  See
  \code{\link{methods}} for more details.
}


\value{
  An object of class \code{"mcsSubset"}, i.e. a list
  with the following components:
  \item{weights}{Weights.}
  \item{offset}{Offset.}
  \item{nobs}{Number of observations.}
  \item{nvar}{Number of variables (not including intercept, if any).}
  \item{x.names}{Names of all design variables.}
  \item{y.name}{Name of response variable.}
  \item{include}{Indexes of variables forced in.}
  \item{exclude}{Indexes of variables forced out.}
  \item{intercept}{\code{TRUE} if regression has an intercept term;
    \code{FALSE} otherwise.}
  \item{penalty}{AIC penalty.}
  \item{nbest}{Number of best subsets.}

  When \code{penalty == 0}:
  \item{size}{Subset sizes.}
  \item{tolerance}{Tolerance vector.}
  \item{rss}{A two dimensional numeric \code{nbest x nvar} array.}
  \item{which}{A three dimensional logical \code{nvar x nbest x nvar}
    array.}

  The entry \code{rss[i, n]} corresponds to the RSS of the \code{i}-th
  best subset model of size \code{n}.  The entry \code{which[j, i, n]}
  has value \code{TRUE} if the \code{i}-th best subset model of size
  \code{n} contains the \code{j}-th variable.

  When \code{penalty != 0}:
  \item{tolerance}{Tolerance value.}
  \item{rss}{A one dimensional numeric array of length \code{nbest}.}
  \item{aic}{A one dimensional numeric array of length \code{nbest}.}
  \item{which}{A two dimensional logical \code{nvar x nbest} array.}

  The entries \code{rss[i]} and \code{aic[i]} correspond to the RSS and
  AIC of the \code{i}-th best subset model, respectively.  The entry
  \code{which[j, i]} is \code{TRUE} if the \code{i}-th best subset model
  contains variable \code{j}.
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

\seealso{\code{\link{summary}}, \link{methods}.}


\examples{
## load data (with logs for relative potentials)
data("AirPollution", package = "mcsSubset")

#################
## basic usage ##
#################

## canonical example: fit best subsets
xs <- mcsSubset(mortality ~ ., data = AirPollution)

## visualize RSS
plot(xs)

## summarize
summary(xs)

## plot summary
plot(summary(xs))

## forced inclusion/exclusion of variables
xs <- mcsSubset(mortality ~ ., data = AirPollution,
                include = "noncauc", exclude = "whitecollar")

## or equivalently
xs <- mcsSubset(mortality ~ ., data = AirPollution,
                include = 10, exclude = 11)
summary(xs)

##########################
## find best BIC models ##
##########################

## find 10 best subset models
xs <- mcsSubset(mortality ~ ., data = AirPollution,
                penalty = "BIC", nbest = 10)

## summarize
summary(xs)

## visualize BIC and RSS
plot(summary(xs))
}

\keyword{regression}
