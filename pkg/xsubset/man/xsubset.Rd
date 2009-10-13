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
  size = NULL, pradius = NULL, tol = 0, na.action = na.omit, \dots)
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
  \item{weights}{optional numeric vector of weights (currently not
    implemented yet).}
  \item{offset}{optional numeric vector with an a priori known component to be
    included in the linear predictor.}
  \item{model, x}{logicals. If \code{TRUE} (in the \code{formula} and \code{lm}
    method) the model frame or model matrix, respectively, is returned (similarly
    for the dependent variable \code{y}).}
  \item{size}{numeric. Vector of subset sizes (number of variables, excluding
    the intercept) to be evaluated in the subset selection. The default is to
    include all sizes from zero the number of available variables.}
  \item{pradius}{numeric. Preordering radius, by default \code{round(NROW(x)/3)}.}
  \item{tol}{numeric. Tolerance parameter.}
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
  are usually attained for approximately \code{pradius}\eqn{=n/3}.
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
  \item{bic}{a vector with BIC values of best subsets corresponding to \code{size},}
  \item{best}{best BIC subset \code{size},}
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
for(i in 12:14) AirPollution[[i]] <- log(AirPollution[[i]])

## fit best subsets
xs <- xsubset(mortality ~ ., data = AirPollution)

## visualize/summarize RSS and BIC
plot(xs)
summary(xs)

## extract information (by default for best BIC model)
deviance(xs)
logLik(xs)
AIC(xs)
AIC(xs, k = log(nrow(AirPollution)))

## extract information for all subsets fitted
deviance(xs, size = 0:15)
AIC(xs, size = 0:15)
AIC(xs, size = 0:15, k = log(nrow(AirPollution)))

## refit model (best BIC by default)
lm5 <- refit(xs)
summary(lm5)
## (Note that the p-values are not valid due to model selection.)
}

\keyword{regression}
