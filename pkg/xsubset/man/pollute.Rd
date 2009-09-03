\name{pollute}

\docType{data}

\alias{pollute}

\title{Pollute dataset}

\description{The Pollute dataset.}

\usage{data(pollute)}

\format{Dataframe with 60 observations and 16 variables.}

\details{

  This is the pollution data so loved by writers of papers on ridge
  regression.

  Variables in order:
  \describe{
    \item{PREC}{Average annual precipitation in inches.}
    \item{JANT}{Average January temperature in degrees F.}
    \item{JULT}{Same for July}
    \item{OVR65}{\% of 1960 SMSA population aged 65 or older.}
    \item{POPN}{Average household size.}
    \item{EDUC}{Median school years completed by those over 22.}
    \item{HOUS}{\% of housing units which are sound & with all facilities.}
    \item{DENS}{Population per sq. mile in urbanized areas, 1960.}
    \item{NONW}{\% non-white population in urbanized areas, 1960.}
    \item{WWDRK}{\% employed in white collar occupations.}
    \item{POOR}{\% of families with income < \$3000.}
    \item{HC}{Relative hydrocarbon pollution potential.}
    \item{NOX}{Same for nitric oxides.}
    \item{SO2}{Same for sulphur dioxide.}
    \item{HUMID}{Annual average \% relative humidity at 1pm.}
    \item{MORT}{Total age-adjusted mortality rate per 100,000.}
  }

}

\source{\url{http://lib.stat.cmu.edu/datasets/pollution}}

\references{

  McDonald, G. C. and Schwing, R. C. (1973).  Instabilities of
  regression estimates relating air pollution to mortality.
  Technometrics, 15, 463-482.

  Miller, A. J. (2002).  Subset Selection in Regression.  New York:
  Chapman and Hall.  Related software can be found online at
  \url{http://users.bigpond.net.au/amiller/}.

}

\keyword{datasets}
