\name{house}
\alias{house}

\title{The House Data}

\description{The House dataset. FIXME.}

\usage{
  data("house")
}

\format{
A data frame containing 26 observations on 14 variables.
\describe{
  \item{price}{FIXME}
  \item{bdr}{bedrooms, FIXME}
  \item{flr}{FIXME}
  \item{fp}{FIXME}
  \item{rms}{rooms, FIXME}
  \item{st}{FIXME}
  \item{lot}{lotsize, FIXME}
  \item{tax}{FIXME}
  \item{bth}{FIXME}
  \item{con}{FIXME}
  \item{gar}{garage, FIXME}
  \item{cdn}{FIXME}
  \item{l1}{FIXME}
  \item{l2}{FIXME}
}
}

\source{\url{http://utstat.toronto.edu/~srivasta/sasreg/CH3E1.SAS}}

\references{
  Sen, A. and Srivastava, M. (1990). \emph{Regression Analysis}.
  Springer-Verlag, New York.
}

\examples{
## load data and transform some variables to logs
data("house", package = "xsubset")
for(i in c(1, 3, 7)) house[[i]] <- log(house[[i]])

## fit subsets
xs <- xsubset(price ~ ., data = house)
plot(xs)
summary(xs)

## refitting best BIC model
summary(refit(xs))
}

\keyword{datasets}
