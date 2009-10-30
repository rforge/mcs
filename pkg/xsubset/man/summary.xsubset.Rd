\name{summary.xsubset}

\alias{summary.xsubset}
\alias{print.summary.xsubset}
\alias{plot.xsubset}

\alias{AIC.xsubset}
\alias{deviance.xsubset}
\alias{logLik.xsubset}
\alias{refit.xsubset}
\alias{coef.xsubset}
\alias{vcov.xsubset}
\alias{residuals.xsubset}
\alias{fitted.xsubset}
\alias{model.frame.xsubset}
\alias{model.matrix.xsubset}

\title{Methods for xsubset objects}

\description{
  Methods for best subset regression objects fitted by \code{\link{xsubset}}.
}

\usage{
\method{summary}{xsubset}(object, size = NULL, \dots)

\method{refit}{xsubset}(object, size = NULL, \dots)

\method{plot}{xsubset}(x, size = NULL, legend = TRUE,
  xlab = "Number of regressors in model", ylab = "",
  main = "BIC and residual sum of squares",
  col = c(1, 4), lty = 1, type = "b", \dots)
}

\arguments{
  \item{object, x}{an object of class \code{xsubset}.}
  \item{size}{numeric. Vector of subset sizes, defaulting to all subset
    sizes evaluated in \code{object} (or \code{x}, respectively).}
  \item{legend}{logical. Should a legend be drawn?}
  \item{xlab, ylab, main}{character with plot annotation.}
  \item{col, lty}{graphical parameters, expanded to length 2 for BIC
    and RSS, respectively.}
  \item{type}{plot type.}
  \item{\dots}{further arguments passed to methods}
}

\details{
  Many methods that work for \code{lm} objects also work for \code{xsubset}
  objects. The idea is always that the linear model is refitted for the
  specified \code{size} and then the information is extracted. If some method
  is not available, \code{refit} can always be called explicitly which returns
  a \code{lm} object. For some
  methods, \code{size} can also be a vector of sizes that should be assessed.
  The default for \code{size} is typically the best BIC model, sometimes (as
  in the \code{summary} and \code{plot} method) all sizes available in the 
  \code{xsubset} object are assessed. 

  Currently, the methods to the generic functions include
  \code{\link[base]{print}}, \code{\link[base]{summary}}, \code{\link[stats]{coef}}, 
  \code{\link[stats]{vcov}}, \code{\link[stats]{logLik}}, \code{\link[stats]{residuals}}, 
  \code{\link[stats]{fitted}}, \code{\link[stats]{model.frame}}, \code{\link[stats]{model.matrix}},
  \code{\link[stats]{AIC}}, \code{\link[stats]{logLik}}, \code{\link[stats]{deviance}}.  
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
deviance(xs, size = 1:16)
AIC(xs, size = 1:16)
AIC(xs, size = 1:16, k = log(nrow(AirPollution)))

## refit model (best BIC by default)
lm5 <- refit(xs)
summary(lm5)
## (Note that the p-values are not valid due to model selection.)

## extract model information
coef(xs, size = 3)
vcov(xs, size = 3)
residuals(xs, size = 3)
fitted(xs, size = 3)
model.matrix(xs, size = 3)
}

\keyword{regression}
