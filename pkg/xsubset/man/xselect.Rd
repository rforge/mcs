\name{xselect}

\alias{xselect}

\title{Best-Subset Selection}

\description{

  Best-subset selection via all-subsets regression in linear
  regression models.

}

\usage{

  xselect(nobs, nvar, x, y, mark, criterion, tolerance, pradius, nbest)

}

\arguments{

  \item{nobs, nvar}{Number of observations (rows) and regressors
  (columns), respectively.}

  \item{x, y}{Model matrix and dependent variable, respectively.}

  \item{mark}{The first \code{mark} columns of \code{x} are forced into
    the model.}

  \item{criterion}{Information criterion.  Available criteria are
    \code{"RSS"} (default), \code{"AIC"}, \code{"BIC"} and \code{"Cp"}.}

  \item{tolerance}{The tolerance, where \code{tolerance[j]} is the
    tolerance employed for subsets of size \code{j}.}

  \item{pradius}{Preordering radius, by default \code{round(nreg/3)},
    where \code{nreg} is the number of regressors.}

  \item{nbest}{Number of subsets to report for every size.}

}

\details{

  This is a low-level routine.  For a more convenient calling interface
  see \code{\link{xsubset}}.

}


\value{

  A list with the following components:

  \item{nobs, nvar}{TODO}

  \item{mark}{TODO}
  
  \item{criterion}{TODO}

  \item{tolerance}{TODO}

  \item{nbest}{TODO}

  \item{value}{TODO}

  \item{which}{TODO}

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

\seealso{\code{\link{xsubset}}, \code{\link{summary}}, \link{methods}.}

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
