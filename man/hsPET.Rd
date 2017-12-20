\name{hsPET}
\alias{hsPET}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of the potential evapotranspiration according to the Hargreaves-Samani method
}
\description{
The function uses the Hargreaves-Samani equation to estimate the potential evapotranspiration from minimum and maximum temperature records, provided the latitude of the target location. Optionally, the extraterrestrial solar radiation is also returned.
}
\usage{
hsPET(date, lat = 42, t.min, t.max, ret.rad = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{date}{
A vector of dates. This includes any format directly coercible to the \code{POSIXlt} time class (see details of function \code{\link{fwi}}). 
}
  \item{lat}{
Latitude of the location, in decimal degrees
}
  \item{t.min}{
Daily series of minimum temperature
}
  \item{t.max}{
Daily series of maximum temperature
}
  \item{ret.rad}{
Logical flag. Should the estimated extraterrestrial radiation values be returned. Default to \code{FALSE}
}
}

\value{
A list containing the following items:
\item{date }{A \code{POSIXlt} vector of dates}
\item{HSevap }{A numeric vector of the estimated potential evapotranspiration, in mm per unit time}

Optionally, if argument \code{ret.rad = TRUE}, a third element is included in the output list: 

\item{Rad }{A numeric vector of the estimated extraterrestrial radiation, in MJ m-2 per unit time}
%% ...
}
\references{
Hargreaves, G.H. and Allen, R.G. (2003) History and evaluation of Hargreaves evapotranspiration equation. \emph{Journal of Irrigation and Drainage Engineering-ASCE}, 129, 53-63.
}
\author{
J. Bedia and M. Iturbide
}
\note{
The difference between the maximum and minimum air temperature (Tmax - Tmin) can be used as an indicator of the fraction of extraterrestrial radiation that reaches the Earth's surface. This principle has been utilized by Hargreaves and Samani to develop estimates of PET using only air temperature data. In general, temperature-based methods remain empirical and require local calibration in order to achieve satisfactory results. A possible exception is the 1985 Hargreaves' method which has shown reasonable ETo results with a global validity. 

Daily estimates by the Hargreaves equations are subject to error caused by the influence of the temperature range caused by the movement of weather fronts and by large variations in wind speed or cloud cover. Therefore, the Hargreaves methods are recommended for use with five-day or longer time steps.
}


\examples{
data(meteoSpain)

daily.pet <- hsPET(date=meteoSpain$date, lat=41.9, t.min=meteoSpain$tn,
            t.max=meteoSpain$tx, ret.rad=FALSE)
str(daily.pet)

# Daily PET values
plot(daily.pet$date, daily.pet$HSevap, ty='l')

# Mean annual values of daily PET
x <- tapply(daily.pet$HSevap, daily.pet$date$year, mean)
plot(unique(daily.pet$date$year + 1900), x, ty='o', ylab = "HS-PET (mm day-1)",
            xlab = "year")
abline(h=mean(x), lty=2)
# End
}
