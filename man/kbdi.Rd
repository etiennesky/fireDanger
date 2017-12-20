\name{kbdi}
\alias{kbdi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of the Keetch-Byram and McArthur's Drought Indices, and the McArthur's Forest Fire Danger Index (FFDI)
}
\description{
The function computes three popular indices used for the assessment of forest fire danger potential: the Keetch-Byram (KBDI) and McArthur's Drought Indices, and the McArthur's Forest Fire Danger Index.
}
\usage{
kbdi(date, t, p, h, w, wrs = 150, start.date = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{date}{
A vector of dates. This includes any format directly coercible to the \code{POSIXlt} time class (see details of function \code{\link{fwi}}
}
  \item{t}{
Temperature time series, in degrees C
}
  \item{p}{
Daily rainfall time series, in mm
}
  \item{h}{
Relative humidity time series, in \%
}
  \item{w}{
Wind speed time series, in km/h
}
  \item{wrs}{
Optional. Minimum weekly rainfall (mm accumulated in a week period), assumed to lead soil moisture to field capacity (see details)
}
  \item{start.date}{
Optional. Start date of KBDI computation (see details). Same format as argument \code{date}.
}
}
\details{
The Keetch-Byram index (KBDI) number indicates the amount of net rainfall (in mm) that is required to reduce the index to zero. The initialization of KBDI usually involves setting it to zero after a period of substantial precipitation. Thus, it is necessary to go back in time until a day is reached on which it is reasonably certain that the upper soil layers were saturated to initialize the index. In snowfree areas, the authors suggest 6 to 8 inches (152.4 to 203.2 mm) in a period of a week. In areas of heavy snowfall, the best choice is probably starting the index just after the snow melts in spring. For this reason, the function has two alternative arguments: \code{wrs}, for the precipitation-based approach, or \code{start.date}, in order to set a fixed date to start the index. 

KBDI is an input for the subsequent computation of the McArthur's drought code and FFDI.
}
\value{
The function returns a dataframe with the following columns:
\item{date}{A \code{POSIXlt} vector of dates}
\item{net.precip}{A numeric vector with the net precipitation (in mm)}
\item{KBDI}{A numeric vector with the Keetch-Byram drought index values}
\item{McArthurDF}{A numeric vector with the MacArthur's drought factor values}
\item{FFDI}{A numeric vector with the MacArthur's Forest Fire Danger Index values}
}
\references{
Dennison, P.E., Roberts, D.A., Thorgusen, S.R., Regelbrugge, J.C., Weise, D. and Christopher, L. (2003) Modeling seasonal changes in live fuel moisture and equivalent water thickness using a cumulative water balance index. \emph{Remote Sensing of Environment}, 88, 442-452.

Dowdy, A.J., Mills, G.A., Finkele, K. and de Groot, W. (2010). Index sensitivity analysis applied to the Canadian Forest Fire Weather Index and the McArthur Forest Fire Danger Index. \emph{Meteorological Applications}, 17, 298-312.

Keetch, J.J. and Byram, G.M. (1968). A drought index for forest fire control. Tech. Rep., USDA Forest Service Research Paper SE-38, North Carolina, USA.
}
\author{
J. Bedia
}
\note{
The physical theory for the KBDI (Keetch and Byram 1968) is based on a number of assumptions. The first assumption is that soil moisture is at field capacity with a water depth equivalent to about 200mm. The second assumption is that the rate of moisture loss in an area depends on the vegetation cover in the area, and vegetation density is a function of the mean annual rainfall. Hence, daily transpiration is approximated by an inverse exponential function of the mean annual rainfall. Finally, the evaporation rate of soil moisture with time is assumed to be an estimation of relative evapotranspiration from exponential function of the daily maximum air temperature. Sensitivity analyses from earlier researchers have revealed that KBDI decays exponentially with an assumed maximum soil water deficit and is sensitive to the daily maximum air temperature (Dennison \emph{et al.} 2003). Its values range from 0 to 800 (in hundredths of inch, equivalent to 203.2mm), with 800 indicating extreme drought and zero indicating saturated soil.

On the other hand, McArthur's drought factor equation incorporates the KBDI, and it is in turn used to compute the McArthur's forest fire danger index (FFDI), used in Australia for several decades, generally providing a good indication of the difficulty of fire suppression over a wide range of conditions.
}

\seealso{
\code{\link{fwi}}, \code{\link{angstronIndex}}, \code{\link{nesterovIndex}}, \code{\link{fmi}}, for alternative drought/fire danger indices.
}
\examples{

data(fwiSpain)
x <- fwiSpain$observed$Salamanca
# Initialization of the KBDI the first day of the series, just as an example
calc <- kbdi(date=x$fecha, t=x$T, p=x$P, h=x$H, w=x$W, start.date=x$fecha[1])
par(mfrow=c(2,2))
plot(x$fecha, calc$net.precip, ty='l', main = "Net precipitation",
            ylab = "Net Precip", xlab = "date")
plot(x$fecha, calc$KBDI, ty='l', main = "Keetch-Byram Drought Index",
            ylab = "KBDI", xlab = "date")
plot(x$fecha, calc$McArthurDF, ty='l', main ="McArthur's Drought Factor",
            ylab = "MADF", xlab = "date")
plot(x$fecha, calc$FFDI, ty='l', main = "McArthur's Forest Fire Danger Index",
            ylab = "FFDI", xlab = "date")
#End
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
