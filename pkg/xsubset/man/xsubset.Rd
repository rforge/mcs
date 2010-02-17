\name{xsubset}

\alias{xsubset}
\alias{xsubset.matrix}
\alias{xsubset.data.frame}
\alias{xsubset.formula}
\alias{xsubset.lm}
\alias{xsubset.default}

\alias{print.xsubset}
\alias{plot.xsubset}

\title{All-Subsets Regression}

\description{All-subsets regression for linear models.}

\usage{

  xsubset(object, \dots)

  \method{xsubset}{matrix}(object, y, \ldots)

  \method{xsubset}{data.frame}(object, y = 1, \dots)

  \method{xsubset}{formula}(object, data = NULL, row.subset = NULL,
    weights = NULL, na.action = na.omit, model = TRUE, x = FALSE,
    y = FALSE, contrasts = NULL, offset = NULL, \ldots)

  \method{xsubset}{lm}(object, \dots)

  \method{xsubset}{default}(object, include = NULL, exclude = NULL,
    size = NULL, tolerance = 0, \ldots, pradius = NULL, precision = "double")

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

  \item{size}{Vector of subset sizes (number of variables including the
    intercept, if any) to be evaluated.  The default is to include all
    possible sizes.  If not specified, the value of \code{size} is
    determined from \code{include} and \code{exclude}.}

  \item{pradius}{Preordering radius, by default \code{round(nreg/3)},
    where \code{nreg} is the number of regressors.}

  \item{tolerance}{The tolerance, where \code{tolerance[i]} is the
    tolerance employed for subset size \code{size[i]}.}

  \item{precision}{Use either single (\code{"single"}) or double
    (\code{"double"}) precision arithmetics.}

  \item{\dots}{Ignored.}

}


\details{

  The function \code{xsubset} computes all variable-subsets for linear
  regression models.  The function is generic and provides various
  methods to conveniently specify the regressor and response variable.
  Either the standard \code{formula} interface as in
  \code{\link[stats]{lm}} can be used, or the information can be
  extracted from an already fitted \code{lm} object.  The regressor
  matrix and response variable can also be passed in directly.  In all
  cases, the default method is eventually called.

  The number of observations (rows) has to be larger than the number of
  variables (columns).  By default, the best subset model is derived for
  each number of variables, but the scope of the search can be limited
  to certain subset sizes by setting \code{size}.

  \bold{Preordering:} The function will preorder the variables to reduce
  execution times if \code{pradius}\eqn{>0} (\code{pradius}\eqn{\le n}),
  where \eqn{n} is the number of regressors.  Best execution times are
  usually attained for approximately \code{pradius}\eqn{=n/3}, which is
  the default value.

  \bold{Heuristic:} If the exact solution is not needed, a tolerance
  \code{tolerance}\eqn{>0} can be specified for each subset size
  (expanded if necessary) to speed up the procedure.  The solution
  satisfies RSS\eqn{(k)\le(1+\tau_k)}RSS\eqn{^*(k)}, where RSS\eqn{(k)}
  is the residual sum of squares (RSS) of the subset of size \eqn{k},
  RSS\eqn{^*} is the (theoretical) RSS of the exact solution, and
  \eqn{\tau_k} the tolerance for subset size \eqn{k}.  Specifically,
  \code{tolerance[i]} is the tolerance for subsets of size
  \code{size[i]}.
  
  A set of standard extractor functions for fitted model objects is
  available for objects of class \code{"xsubset"}.  See
  \code{\link{summary.xsubset}} for more details.

}

\value{

  An object of class \code{"xsubset"}, i.e. a list with the following
  components:

  \item{call}{The original function call.}

  \item{lm}{An object of type \code{lm} representing the full model.}
  
  \item{nreg}{The (original) number of regressors.}

  \item{include}{The (computed) index of included variables.}

  \item{exclude}{The (computed) index of excluded variables.}

  \item{size}{A vector giving the evaluated subset sizes (computed).}

  \item{which}{A list indicating the selected variables.  The component
    \code{which[[i]]} is a numeric vector of length \code{size[i]}
    holding the indexes of the variables that form the selected subset
    of size \code{size[i]}.}

  \item{rss}{A vector indicating the residual sum of squares (RSS) of
    the selected subsets.  The element \code{rss[i]} is the RSS of the
    subset of size \code{size[i]}.}

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

\seealso{\code{\link{summary.xsubset}}}

\examples{

## load data (with logs for relative potentials)
data("AirPollution", package = "xsubset")

#################
## basic usage ##
#################

## simplest form
xs <- xsubset(AirPollution, y = "mortality")

## which, in this case, is equivalent to
## (by default, 'y' is last column)
xs <- xsubset(as.matrix(AirPollution))

## canonical example: fit best subsets
xs <- xsubset(mortality ~ ., data = AirPollution)

## visualize RSS
plot(xs)

## visualize/summarize RSS and AIC
sx_aic <- summary(xs)
print(sx_aic)
plot(sx_aic)

## extract information (subset size 7)
logLik(xs, 7)
coef(xs, 7)

## extract information for all subsets fitted
deviance(xs)
AIC(xs)
AIC(xs, k = log(nrow(AirPollution)))

## refit best model
lm10 <- refit(xs, summary(xs)$best)
summary(lm10)
## (Note that the p-values are not valid due to model selection.)


#####################
## further options ##
#####################

## summarize BIC
sx_bic <- summary(xs, aic.penalty = log(nrow(AirPollution)))
print(sx_bic)
plot(sx_bic)

## forced inclusion/exclusion of variables in the model
xs_fix <- xsubset(mortality ~ ., data = AirPollution,
  include = c("(Intercept)", "noncauc"), exclude = "whitecollar")
## or equivalently
xs_fix <- xsubset(mortality ~ ., data = AirPollution,
                  include = c(1, 10), exclude = 11)
plot(xs_fix)
plot(xs, size = 2:15)

## weights (artificial example)
set.seed(123)
w <- runif(nrow(AirPollution), max = 3)
w[w < 0.3] <- 0
xs_wls <- xsubset(mortality ~ precipitation + temperature1 + age,
                  weights = w, data = AirPollution)

}

\keyword{regression}
