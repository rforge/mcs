\name{ChicagoHouses}
\alias{ChicagoHouses}

\title{Chicago House Price Data}

\description{Data on house prices in Chicago,
originally provided by Ms. Terry Tasch of Long-Kogan Realty, Chicago.}

\usage{
  data("ChicagoHouses")
}

\format{
A data frame containing 26 observations on 14 variables.
\describe{
  \item{price}{Selling price of house in thousands of dollars.}
  \item{bedrooms}{Number of bedrooms.}
  \item{floor}{Floor space in square feet (computed from dimensions
    of each room and then augmented by 10 percent).}
  \item{fireplaces}{Number of fireplaces.}
  \item{rooms}{Number of rooms.}
  \item{stormwindows}{Factor. Does the house have storm windows?}
  \item{front}{Front footage of lot in feet.}
  \item{tax}{Annual taxes.}
  \item{bathrooms}{Number of bathrooms.}
  \item{construction}{Factor indicating type of construction (frame or brick).}
  \item{garage}{Garage size (none, one car, etc.).}
  \item{condition}{Factor indicating condition (needs work or otherwise).}
  \item{location}{Factor indication location (zones A, B, or other).}
}
}

\source{Sen and Srivastava (1990), Exhibit 2.2, p. 32.}

\references{
  Sen, A. and Srivastava, M. (1990). \emph{Regression Analysis}.
  Springer-Verlag, New York.
}

\examples{
## load data and transform some variables to logs
data("ChicagoHouses", package = "xsubset")
for(i in c(1, 3, 7)) ChicagoHouses[[i]] <- log(ChicagoHouses[[i]])

## fit subsets
xs <- xsubset(price ~ ., data = ChicagoHouses)
plot(xs)
summary(xs)

## refitting best BIC model
summary(refit(xs, size = summary(xs)$best))
}

\keyword{datasets}
