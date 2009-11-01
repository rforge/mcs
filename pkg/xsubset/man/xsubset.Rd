\name{xsubset}

\alias{xsubset}
\alias{xsubset.default}
\alias{xsubset.formula}
\alias{xsubset.lm}
\alias{print.xsubset}

\title{Variable Subset Selection}

\description{
  Variable subset selection via best subset regression in linear regression models.
}

\usage{
xsubset(object, \dots)

\method{xsubset}{formula}(formula, data, subset, na.action, weights, offset, \dots,
  model = TRUE, y = TRUE, x = FALSE)
  
\method{xsubset}{lm}(object, \dots, model = TRUE, y = TRUE, x = FALSE)

\method{xsubset}{default}(object, y, weights = NULL, offset = NULL,
  size = NULL, pradius = NULL, tol = 0,
  criterion = NULL, include = NULL, exclude = NULL,
  na.action = na.omit, \dots)
}

\arguments{
  \item{object}{an object. For the \code{lm} method, this needs to be of class \code{lm},
    for the default method it needs to be the regressor matrix containing all
    variables to be included in the subset selection. If an intercept should
    be included it has to be in the first column of \code{x}.}
  \item{y}{numeric (in the default method). Dependent variable in the linear
    regression.}
  \item{formula}{symbolic description of the variables included in the
    subset selection.}
  \item{data, subset, na.action}{arguments controlling formula processing
    via \code{\link[stats]{model.frame}}.}
  \item{weights}{optional numeric vector of weights,}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor.}
  \item{model, x}{logicals. If \code{TRUE} (in the \code{formula} and \code{lm}
    method) the model frame or model matrix, respectively, is returned (similarly
    for the dependent variable \code{y}).}
  \item{size}{numeric. Vector of subset sizes (number of variables, including
    the intercept if any) to be evaluated in the subset selection. The default is to
    include all possible sizes (as implied by \code{include} and \code{exclude}).}
  \item{pradius}{numeric. Preordering radius, by default \code{round(NROW(x)/3)}.}
  \item{tol}{numeric. Tolerance when minimizing the residual sum of squares (RSS).
    When \code{tol = 0} the RSS is minimized, otherwise a solution may be found
    which is up to \code{tol} times the minimal RSS worse.}
  \item{criterion}{numeric or function. A specification of the model selection
    criterion: It can either be a a numeric penalty for calling
    \code{AIC(loglik, k = criterion)}, where \code{logLik} is an object of class
    \code{"logLik"}, see \code{\link[stats]{logLik}} and \code{\link[stats]{AIC}}.
    Or it can be a function for calling \code{criterion(loglik)}.
    The default is to use the BIC, i.e., \code{AIC(logLik, k = log(nobs))} where
    \code{nobs} is the number of observations (with \code{weights > 0}).
    In both cases the best model is chosen by the minimum criterion.}
  \item{include, exclude}{index vectors. Index of variables that are forced to
    be included or excluded from the model, respectively. The indexes pertain
    to the columns of the regressor matrix and can be (as usual) integer indexes,
    character names, or a logical vector selecting the desired columns. By default
    no variables are excluded and the intercept (if any) is forced to be included.}
  \item{\dots}{further arguments passed to methods}
}

\details{
  The function \code{xsubset} computes the best variable-subset models
  for linear regression models. The function is generic and provides
  various methods for conveniently specifying the regressor matrix and
  response. Either a standard \code{formula} interface as in
  \code{\link[stats]{lm}} can be used, or the information can be extracted
  from an already fitted \code{lm} object. In both cases, the default
  method is called as the workhorse function.

  The number of observations (rows) has to be larger than the number of
  variables (columns). The best subset model is derived for each number of
  variables by default, the search scope can be limited by setting \code{size}
  to a selection of potential subset sizes.

  The function will preorder the variables to reduce execution time if
  \code{pradius}\eqn{>0} (\code{pradius}\eqn{\le n}).  Best execution times
  are usually attained for approximately \code{pradius}\eqn{=n/3} which
  is hence chosen as the default.
  If the exact solutions are not needed, and to prevent the function
  from performing an exhaustive search, a tolerance \code{tol}\eqn{>0}
  can be specified for each subset \code{size} (expanded if necessary).
  
  A set of standard extractor functions for fitted model objects is available for
  objects of class \code{"xsubset"}, including methods to the generic functions
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{fitted}}, \code{\link[stats]{model.frame}}, \code{\link[stats]{model.matrix}},
  \code{\link[stats]{AIC}}, \code{\link[stats]{logLik}}, \code{\link[stats]{deviance}}.  
  See \code{\link{summary.xsubset}} for more details on all methods.
}

\value{
  An object of class \code{"xsubset"}, i.e., a list with components including
  \item{call}{the original function call,}
  \item{size}{a vector with evaluated variable subset sizes,}
  \item{rss}{a vector with residual sum of squares of best subset corresponding to \code{size},}
  \item{which}{a list with variable indices of best subsets corresponding to \code{size},}
  \item{include}{a vector with integer indexes of the variables (columns of \code{x})
    forced to be included in the model,}
  \item{exclude}{a vector with integer indexes of the variables (columns of \code{x})
    forced to be excluded from the model,}
  \item{criterion}{a matrix with one row of criterion values (by default BIC) of the
    best subsets corresponding to \code{size},}
  \item{best}{best subset \code{size} (for which \code{criterion} is minimal),}
  \item{x}{model matrix with all variables (always in default method, only if
    \code{x = TRUE} in the other methods),}
  \item{y}{the response count vector (always in default method, only if
    \code{y = TRUE} in the other methods),}
  \item{weights}{the case weights used (if any),}
  \item{offset}{the offset vector used (if any),}
  \item{nobs}{number of observations,}
  \item{nreg}{number of variables/regressors,}
  \item{intercept}{logical indicating whether an intercept is included in the regression,}
  \item{pradius}{preordering radius used,}
  \item{tol}{tolerance parameter vector used,}
  \item{nvisited}{number of visited nodes in the search algorithm,}
  \item{na.action}{NA action used (if any),}
  \item{formula}{the original formula (in \code{formula} and \code{lm} methods),}
  \item{terms}{the associated terms (in \code{formula} and \code{lm} methods),}
  \item{model}{the full model frame (only if \code{model = TRUE} in the
    \code{formula} and \code{lm} methods).}
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

## fit best subsets
xs <- xsubset(mortality ~ ., data = AirPollution)

## visualize/summarize RSS and BIC
plot(xs)
summary(xs)

## extract information (by default for best BIC model)
logLik(xs)
coef(xs)

## extract information for all subsets fitted
deviance(xs)
AIC(xs)
AIC(xs, k = log(nrow(AirPollution)))

## refit model (best BIC by default)
lm6 <- refit(xs)
summary(lm6)
## (Note that the p-values are not valid due to model selection.)

#####################
## further options ##
#####################

## run full analysis with AIC model selection
xs_aic <- xsubset(mortality ~ ., data = AirPollution, criterion = AIC)
## or equivalently
xs_aic <- xsubset(mortality ~ ., data = AirPollution, criterion = 2)

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
