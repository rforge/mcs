\name{IbkTemperature}
\alias{IbkTemperature}

\title{Temperature Observations and Numerical Weather Predictions for Innsbruck}

\description{
00UTC temperature observations and corresponding 24-hour
reforecast ensemble means from the Global Ensemble Forecast System
(GEFS, Hamill et al. 2013) for SYNOP station Innsbruck
Airport (11120; 47.260, 11.357) from 2011-01-01 to 2015-12-31.
}

\usage{data("IbkTemperature")}

\format{
A data frame containing 1824 daily observations/forecasts for 42 variables.
The first column (\code{temp}) contains temperature observations at
00UTC (coordinated universal time) in degree Celsius, columns 2--37 are
24-hour lead time GEFS reforecast ensemble means for different variables (see below).
Columns 38--42 are deterministic time trend/season patterns.
\describe{
  \item{temp}{Observed temperature at Innsbruck Airport.} 
  \item{tp}{Total accumulated precipitation.} 
  \item{t2m}{Temperature at 2 meters.}
  \item{u10m}{U-component of wind at 10 meters.} 
  \item{v10m}{V-component of wind at 10 meters.} 
  \item{u80m}{U-component of wind at 80 meters.} 
  \item{v80m}{U-component of wind at 80 meters.} 
  \item{cape}{Convective available potential energy.} 
  \item{ci}{Convective inhibition.} 
  \item{sdlwrf}{Surface downward long-wave radiation flux.} 
  \item{sdswrf}{Surface downward short-wave radiation flux.} 
  \item{sulwrf}{Surface upward long-wave radiation flux.} 
  \item{suswrf}{Surface upward short-wave radiation flux.}
  \item{ghf}{Ground heat flux.} 
  \item{slhnf}{Surface latent heat net flux.} 
  \item{sshnf}{Surface sensible heat net flux.}
  \item{mslp}{Mean sea level pressure.} 
  \item{psfc}{Surface pressure.} 
  \item{pw}{Precipitable water.}
  \item{vsmc}{Volumetric soil moisture content.} 
  \item{sh2m}{Specific humidity at 2 meters.} 
  \item{tcc}{Total cloud cover.} 
  \item{tcic}{Total column-integrated condensate.} 
  \item{tsfc}{Skin temperature.} 
  \item{tmax2m}{Maximum temperature.} 
  \item{tmin2m}{Minimum temperature.} 
  \item{st}{Soil temperature (0-10 cm below surface).} 
  \item{ulwrf}{Upward long-wave radiation flux.}  
  \item{wr}{Water runoff.} 
  \item{we}{Water equivalent of accumulated snow depth.} 
  \item{wp}{Wind mixing energy.} 
  \item{w850}{Vertical velocity at 850 hPa surface.} 
  \item{t2pvu}{Temperature on 2 PVU surface.} 
  \item{p2pvu}{Pressure on 2 PVU surface.} 
  \item{u2pvu}{U-component of wind on 2 PVU surface.} 
  \item{v2pvu}{U-component of wind on 2 PVU surface.} 
  \item{pv}{Potential vorticity on 320 K isentrope.} 
  \item{time}{Time in years.}
  \item{sin, cos}{Sine and cosine component of annual harmonic pattern.}
  \item{sin2, cos2}{Sine and cosine component of bi-annual harmonic pattern.}
}
}

\source{
Observations: \url{http://www.ogimet.com/synops.phtml.en}

Reforecasts: \url{http://www.esrl.noaa.gov/psd/forecasts/reforecast2/}
}

\references{
  Hamill TM, Bates GT, Whitaker JS, Murray DR, Fiorino M, Galarneau Jr. TJ, Zhu Y, Lapenta W (2013).
    NOAA's Second-Generation Global Medium-Range Ensemble Reforecast Data Set.
    \emph{Bulletin of the American Meteorological Society}, \bold{94}(10), 1553--1565.
    \doi{10.1175/BAMS-D-12-00014.1}
}

\examples{
## load data and omit only a couple of observations with some missing values
data("IbkTemperature", package = "lmSubsets")
IbkTemperature <- na.omit(IbkTemperature)

## fit a simple climatological model for the temperature
## with a linear trend and annual/bi-annual harmonic seasonal pattern
m0 <- lm(temp ~ time + sin + cos + sin2 + cos2, data = IbkTemperature)

## fit a simple MOS with 2-meter temperature forecast in addition
## to the climatological model
m1 <- lm(temp ~ t2m + time + sin + cos + sin2 + cos2, data = IbkTemperature)

## graphical comparison and MOS summary
plot(temp ~ time, data = IbkTemperature, type = "l", col = "darkgray")
lines(fitted(m1) ~ time, data = IbkTemperature, col = "darkred")
lines(fitted(m0) ~ time, data = IbkTemperature, lwd = 2)
summary(m1)

## best subset selection of remaining variables for the MOS
## (i.e., forcing the regressors of m1 into the model)
ms2 <- lmSubsets(temp ~ ., data = IbkTemperature,
  include = c("t2m", "time", "sin", "cos", "sin2", "cos2"))
plot(summary(ms2))
image(ms2, size = 8:20)
## -> Note that soil temperature and maximum temperature are selected
## in addition to the 2-meter temperature

## best subset selection of all variables
ms3 <- lmSubsets(temp ~ ., data = IbkTemperature)
plot(summary(ms3))
image(ms3, size = 2:20)
## -> Note that 2-meter temperature is not selected into the best
## BIC model but soil-temperature (and maximum temperature) are used instead

## refit the best BIC subset selections
m2 <- refit(lmSelect(ms2, penalty = "BIC"))
m3 <- refit(lmSelect(ms3, penalty = "BIC"))

## compare BIC
BIC(m0, m1, m2, m3)

## compare RMSE
sqrt(sapply(list(m0, m1, m2, m3), deviance)/nrow(IbkTemperature))

## compare coefficients
cf0 <- coef(m0)
cf1 <- coef(m1)
cf2 <- coef(m2)
cf3 <- coef(m3)
names(cf2) <- gsub("^x", "", names(coef(m2)))
names(cf3) <- gsub("^x", "", names(coef(m3)))
nam <- unique(c(names(cf0), names(cf1), names(cf2), names(cf3)))
cf <- matrix(NA, nrow = length(nam), ncol = 4,
  dimnames = list(nam, c("m0", "m1", "m2", "m3")))
cf[names(cf0), 1] <- cf0
cf[names(cf1), 2] <- cf1
cf[names(cf2), 3] <- cf2
cf[names(cf3), 4] <- cf3
print(round(cf, digits = 3), na.print = "")
}

\keyword{datasets}
