\name{xselect}

\alias{xselect}
\alias{xselect.matrix}
\alias{xselect.data.frame}
\alias{xselect.formula}
\alias{xselect.lm}
\alias{xselect.default}

\title{Best-Subset Selection}

\description{

  Best-subset selection via all-subsets regression in linear
  regression models.

}

\usage{

  xselect(object, \dots)

  \method{xselect}{matrix}(object, y, \ldots)

  \method{xselect}{data.frame}(object, y = 1, \dots)

  \method{xselect}{formula}(object, data = NULL, row.subset = NULL,
    weights = NULL, na.action = na.omit, model = TRUE, x = FALSE,
    y = FALSE, contrasts = NULL, offset = NULL, \ldots)

  \method{xselect}{lm}(object, \dots)

  \method{xselect}{default}(object, include = NULL, exclude = NULL,
    tolerance = 0, criterion = NULL, \ldots, pradius = NULL, precision =
    "double")

}

\arguments{

  \item{object}{An R object.}

  \item{data, row.subset, weights, na.action, model, x, contrasts,
    offset}{Arguments forwarded to \code{\link[stats]{lm}}.}

  \item{y}{A numeric vector representing the response variable \eqn{y}
    for \code{xsubset.matrix}, the index (numeric or character) of the
    response variable in the data frame for \code{xsubset.data.frame},
    or a logical for \code{xsubset.formula} (see \code{\link[stats]{lm}}
    for details).}

  \item{include, exclude}{Index vectors designating variables that are
    forced into or out of the model, respectively.  The vectors can
    consist of (integer) indexes, (character) names, or (logical) bits
    selecting the desired columns.  All variables are included by
    default.}

  \item{pradius}{Preordering radius, by default \code{round(nreg/3)},
    where \code{nreg} is the number of regressors.}

  \item{tolerance}{The tolerance.}

  \item{criterion}{Criterion name (\code{"AIC"}, \code{"AICc"},
    \code{"AICc2"}, \code{"AICu"}, \code{"BIC"}, or \code{"Cp"}) or
    penalty for modified AIC.}

  \item{precision}{Use either single (\code{"single"}) or double
    (\code{"double"}) precision arithmetics.}

  \item{\dots}{Ignored.}

}

\details{

  These methods serve the same purpose as the \code{xsubset} family of
  methods.  However, a specified criterion is employed during the
  search.  The result of the search is the best subset model over all
  sizes.

  The criterion can be one of "AIC", "AICc", "AICc2", "AICu", "BIC", or
  "mAIC".  The latter is given by the numeric value of the penalty to be
  employed by the modified AIC criterion.

  See \code{\link[stats]{lm}} for more information.

}


\value{

  An object of class \code{"xselect"}, i.e. a list with the following
  components:

  \item{call}{The original function call.}

  \item{lm}{An object of type \code{lm} representing the full model.}
  
  \item{nreg}{The (original) number of regressors.}

  \item{include}{The (computed) index of included variables.}

  \item{exclude}{The (computed) index of excluded variables.}

  \item{best}{The size of the best subset.}

  \item{which}{A list indicating the selected variables.  The component
    \code{which[[i]]} is a numeric vector of length \code{size[i]}
    holding the indexes of the variables that form the selected subset
    of size \code{size[i]}.}

  \item{value}{The value of the modified AIC for the best subset.}

  \item{intercept}{Whether an intercept is included in the regression.}

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

\seealso{\code{\link{summary.xselect}}}

\examples{

## data
data(AirPollution, package = "xsubset")

## AIC model selection
xs_aic <- xselect(mortality ~ ., data = AirPollution, criterion = "AIC")

## or equivalently
xs_aic <- xselect(mortality ~ ., data = AirPollution, criterion = 2)
print(xs_aic)

## summarize
sx_aic <- summary(xs_aic)
print(sx_aic)

}

\keyword{regression}
