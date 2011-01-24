\name{AirPollution}
\alias{AirPollution}

\title{Air Pollution and Mortality}

\description{
  
  Data relating air pollution and mortality, frequently used for illustrations
  in ridge regression and related tasks.
  
}

\usage{

  data("AirPollution")

}

\format{

  A data frame containing 60 observations on 16 variables.

  \describe{

    \item{precipitation}{Average annual precipitation in inches.}

    \item{temperature1}{Average January temperature in degrees
      Fahrenheit.}

    \item{temperature7}{Average July temperature in degrees Fahrenheit.}

    \item{age}{Percentage of 1960 SMSA population aged 65 or older.}

    \item{household}{Average household size.}

    \item{education}{Median school years completed by those over 22.}

    \item{housing}{Percentage of housing units which are sound and with
      all facilities.}

    \item{population}{Population per square mile in urbanized areas,
      1960.}

    \item{noncauc}{Percentage of non-Caucasian population in urbanized
      areas, 1960.}
    
    \item{whitecollar}{Percentage employed in white collar occupations.}

    \item{income}{Percentage of families with income < USD 3000.}

    \item{hydrocarbon}{Relative hydrocarbon pollution potential.}

    \item{nox}{Relative nitric oxides potential.}

    \item{so2}{Relative sulphur dioxide potential.}

    \item{humidity}{Annual average percentage of relative humidity at
      13:00.}

    \item{mortality}{Total age-adjusted mortality rate per 100,000.}

  }

}

\source{\url{http://lib.stat.cmu.edu/datasets/pollution}}

\references{

  McDonald, G.C. and Schwing, R.C. (1973). Instabilities of
  Regression Estimates Relating Air Pollution to Mortality.
  \emph{Technometrics}, \bold{15}, 463--482.

  Miller, A.J. (2002).  \emph{Subset Selection in Regression}.
  New York: Chapman and Hall. 
  Related software can be found online at
  \url{http://users.bigpond.net.au/amiller/}.

}


\examples{

## load data (with logs for relative potentials)
data("AirPollution", package = "xsubset")
for(i in 12:14) AirPollution[[i]] <- log(AirPollution[[i]])

## fit subsets
xs <- xsubset(mortality ~ ., data = AirPollution)
plot(xs)

## summary with BIC
sx <- summary(xs, penalty = log(nrow(AirPollution)))
print(sx)

## refit best model
lm6 <- refit(xs, size = sx$best)
summary(lm6)

}

\keyword{datasets}
