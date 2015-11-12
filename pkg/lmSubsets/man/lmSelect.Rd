\name{lmSelect}


\alias{lmSelect}
\alias{lmSelect.lm}
\alias{lmSelect.formula}
\alias{lmSelect.default}
\alias{lmSelect.lmSubsets}

\alias{print.lmSelect}
\alias{plot.lmSelect}


\title{Best-Subset Regression}


\description{Best-subset regression for ordinary linear models.}

\usage{
lmSelect(object, \dots)

\method{lmSelect}{formula}(formula, \ldots, lm = FALSE)

\method{lmSelect}{lm}(object, \dots)

\method{lmSelect}{default}(object, y, include = NULL, exclude = NULL,
penalty = 0, tolerance = 0, pradius = NULL, nbest = 1, \ldots, .algo =
"hbba")

\method{lmSelect}{lmSubsets}(object, \ldots)
}


\arguments{
  \item{formula, object}{An object of class \code{lm}, \code{formula},
    \code{matrix} or \code{lmSubsets}.}
  \item{y}{The response variable.}
  \item{include, exclude}{Index vectors designating variables that are
    forced in or out of the model, respectively.  The vectors may
    consist of (integer) indexes, (character) names, or (logical) bits
    selecting the desired columns.  The integer indexes correspond to
    the position of the variables in the model matrix; the intercept, if
    any, has index \code{1}.  By default, all variables are included.}
  \item{penalty}{Penalty per parameter (see \code{\link[stats]{AIC}}).}
  \item{tolerance}{A single value indicating the overall tolerance.}
  \item{pradius}{Preordering radius.}
  \item{nbest}{Number of best subsets to report.}
  \item{\dots}{Ignored.}
  \item{lm}{If \code{true}, compute \code{lm} component.}
  \item{.algo}{Internal use.}
}


\details{
  The function \code{lmSelect} computes all variable-subsets regression
  for ordinary linear models.  The function is generic and provides
  various methods to conveniently specify the regressor and response
  variables.  The standard \code{formula} interface (see
  \code{\link[stats]{lm}}) can be used, or the information can be
  extracted from an already fitted \code{lm} object.  The regressor
  matrix and response variable can also be passed in directly.

  The overall (i.e. over all sizes) \code{nbest} best subset models may
  be computed according to an information criterion of the AIC family.
  A single tolerance value may be specified to speed up the search.

  By way of \code{include} and \code{exclude}, variables may be forced
  into or out of the regression, respectively.

  The function will preorder the variables to reduce execution time if
  \code{pradius > 0}.  Good execution times are usually attained for
  approximately \code{pradius = n/3} (default value), where \code{n} is
  the number of regressors after evaluation \code{include} and
  \code{exclude}.

  A set of standard extractor functions for fitted model objects is
  available for objects of class \code{"lmSelect"}.  See
  \code{\link{methods}} for more details.
}


\value{
  An object of class \code{"lmSelect"}, i.e. a list
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
  \item{tolerance}{Tolerance value.}
  \item{rss}{A one dimensional numeric array of length \code{nbest}.}
  \item{aic}{A one dimensional numeric array of length \code{nbest}.}
  \item{which}{A two dimensional logical \code{nvar x nbest} array.}

  The entries \code{rss[i]} and \code{aic[i]} correspond to the RSS and
  AIC of the \code{i}-th best subset model, respectively.  The entry
  \code{which[j, i]} is \code{TRUE} if the \code{i}-th best subset model
  contains variable \code{j}.
}


\seealso{\code{\link{lmSubsets}}, \code{\link{summary}}, \code{\link{methods}}.}


\examples{
## load data (with logs for relative potentials)
data("AirPollution", package = "lmSubsets")

##########################
## find best BIC models ##
##########################

## find 10 best subset models
xs <- lmSelect(mortality ~ ., data = AirPollution,
               penalty = "BIC", nbest = 10)

## summarize
summary(xs)

## visualize BIC and RSS
plot(summary(xs))
}

\keyword{regression}
