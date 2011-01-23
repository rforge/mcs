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

\method{xsubset}{matrix}(object, y = NULL, \ldots)

\method{xsubset}{data.frame}(object, y.which = 1, \dots)

\method{xsubset}{formula}(object, data = NULL, row.subset = NULL,
weights = NULL, na.action = na.omit, model = TRUE, x.return = FALSE,
y.return = FALSE, contrasts = NULL, offset = NULL, \ldots)

\method{xsubset}{lm}(object, \dots)

\method{xsubset}{default}(object, include = NULL, exclude = NULL, size =
NULL, criterion = 0, tolerance = 0, pradius = NULL, nbest = 1, \ldots)

}

\arguments{

  \item{object}{An R object.}

  \item{data, row.subset, weights, na.action, model, x, contrasts,
    offset}{Arguments forwarded to \code{\link[stats]{lm}}.}

  \item{y}{A numeric vector representing the response variable.}

  \item{y.which}{The index (numeric or character) of the response
    variable in the data frame.}

  \item{x.return, y.return}{Logical (see \code{\link[stats]{lm}} for
    details).}

  \item{include, exclude}{Index vectors designating variables that are
    forced in or out of the model, respectively.  The vectors may
    consist of (integer) indexes, (character) names, or (logical) bits
    selecting the desired columns.  By default, all variables are
    included.  The intercept is included if present in the model.}

  \item{size}{Vector of subset sizes (counting the intercept, if any) to
    be evaluated.  By default, all possible sizes (determined by
    \code{include} and \code{exclude}) are computed.}

  \item{criterion}{The selection criterion.  Available criteria are
  \code{"RSS"}, \code{"AIC"}, and \code{"BIC"}.  If \code{criterion} is
  numeric, then it is interpreted as a penalty per parameter (see
  \code{\link[stats]{AIC}}).}

  \item{tolerance}{The tolerance, where \code{tolerance[n]} is the
    tolerance employed for subsets of size \code{n}.}

  \item{pradius}{Preordering radius, by default \code{round(n/3)},
    where \code{n} is the number of regressors.}

  \item{nbest}{Number of best subsets to report for each subset size.}

  \item{\dots}{Ignored.}

}


\details{

  The function \code{xsubset} computes all variable-subsets for linear
  regression models.  The function is generic and provides various
  methods to conveniently specify the regressor and response variables.
  The standard \code{formula} interface (see \code{\link[stats]{lm}})
  can be used, or the information can be extracted from an already
  fitted \code{lm} object.  The regressor matrix and response variable
  can also be passed in directly.  In all cases, the default method is
  eventually called.

  By default, the method computes the \code{nbest} best subset models
  for every subset size, where the "best" models are the models with the
  lowest RSS.  The scope of the search can be limited to certain subset
  sizes by setting \code{size}.  Alternatively, the overall (over all
  sizes) \code{nbest} subset models may be computed according to an
  information criterion of the AIC family specified by \code{criterion}.
  Note that the \code{size} argument is ignored when a criterion other
  than the RSS is employed.

  \bold{Preordering:} The function will preorder the variables to reduce
  execution time if \code{pradius}\eqn{>0} (\code{pradius}\eqn{\le n}),
  where \eqn{n} is the total number of regressors.  Best execution times
  are usually attained for \code{pradius}\eqn{=n/3} approximately, which
  is the default value.

  \bold{Heuristic:} If the exact solution is not needed, a tolerance may
  be specified for each subset size (expanded if necessary) to speed up
  the procedure.  The solution satisfies
  RSS\eqn{(n)\le(1+\tau_n)}RSS\eqn{^*(n)}, where RSS\eqn{(n)} is the
  residual sum of squares (RSS) of the solution subset of size \eqn{n},
  RSS\eqn{^*} is the RSS of the exact solution, and \eqn{\tau_n} the
  tolerance for subset size \eqn{n}.  By default, \eqn{\tau_n=0} for all
  \eqn{n}.  Setting \eqn{\tau_n=+\infty} effectively excludes size
  \eqn{n} from the subset search.  If a selection criterion other than
  RSS is used, a single tolerance value for all subset sizes is
  employed.

  After the input arguments have been processed, the call is forwarded
  to \code{\link{xselect}}.
  
  A set of standard extractor functions for fitted model objects is
  available for objects of class \code{"xsubset"}.  See
  \code{\link{Methods}} for more details.

  Note: The number of observations (rows) must be larger than the number
  of variables (columns).

}

\value{

  An object of class \code{"xsubset"}, i.e. a list containing the
  unprocessed (raw) and processed (as forwarded to \code{xselect})
  arguments, as well as

  \item{value}{If the selection criterion is the RSS, then \code{value}
  is a two-dimensional array with the RSS of the \code{nbest} best
  subset models for all subset sizes.  Otherwise, a one-dimensional
  array with the AIC of the overall best models.}

  \item{which}{If the selection criterion is the RSS, then \code{which}
  is three-dimensional logical array indication the composition of each
  subset in terms of parameters.  Otherwise, a two-dimensional array.}

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

\seealso{\code{\link{xselect}}, \code{\link{summary}}, \link{methods}.}

\examples{

## load data (with logs for relative potentials)
data("AirPollution", package = "xsubset")

#################
## basic usage ##
#################

## simplest form
xs <- xsubset(AirPollution, y.which = "mortality")

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
lm8 <- refit(xs, summary(xs)$best)
summary(lm8)
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
plot(summary(xs, size = 2:15))

## weights (artificial example)
set.seed(123)
w <- runif(nrow(AirPollution), max = 3)
w[w < 0.3] <- 0
xs_wls <- xsubset(mortality ~ precipitation + temperature1 + age,
                  weights = w, data = AirPollution)

}

\keyword{regression}
