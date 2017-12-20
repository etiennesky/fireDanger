\name{angstromIndex}
\alias{angstromIndex}
\title{
Computation of the Swedish Angstrom Drought Index
}
\description{
The function computes Angstrom Index from daily time series of temperature and relative humidity. 
}
\usage{
angstromIndex(t, rh)
}
\arguments{
  \item{t}{A numeric vector of temperature records, in degree C
}
  \item{rh}{
A numeric vector of relative humidity records, in \%
}
}
\details{
The function will not compute days with missing values for any of the input variables, and will jump to the next day with all complete records.
}
\value{
A numeric vector with the values of the index
}
\references{
Willis, C., van Wilgem, B., Tolhurst, K., Everson, C., D'Abreton, P., Pero, L. and Fleming, G. (2001). The development of a national fire danger rating system for South Africa. Department of Water affairs and Forestry, Pretoria, South Africa.
}
\author{
J. Bedia
}
\note{
The Angstrom Index is useful only as an indication of the dryness of the air. It does not account for fuel moisture status or the effect of wind.
The resulting values are categorized as follows:

AI > 4.0 -- Fire occurrence unlikely

4.0 < AI < 2.5 -- Fire conditions unfavorable

2.5 < AI < 2.0 -- Fire conditions favorable

AI < 2.0 -- Fire occurrence very likely
}
\examples{
data(fwiSpain)

ai <- angstromIndex(t=fwiSpain$observed$Salamanca$T,
            rh=fwiSpain$observed$Salamanca$H)
par(mfrow=c(3,1))
plot(fwiSpain$observed$Salamanca$fecha, fwiSpain$observed$Salamanca$T,
            ty='l', main = "temperature")
plot(fwiSpain$observed$Salamanca$fecha, fwiSpain$observed$Salamanca$H,
            ty='l', main = "relHumidity")
plot(fwiSpain$observed$Salamanca$fecha, ai, ty='l',
            main = "Swedish Angstrom drought Index")

# Classification into the standard risk classes
abline(h = c(2,2.5,4), lty = 2)
riskClasses <- findInterval(ai, c(2,2.5,4))
riskCols <- c("red","orange","yellow","green")
for (i in 1:length(sort(unique(riskClasses)))) {
      ind <- which(riskClasses == sort(unique(riskClasses))[i])
      points(fwiSpain$observed$Salamanca$fecha[ind], ai[ind], col = riskCols[i])
}
legend("topleft", c("likely", "favourable", "unfavourable", "unlikely"),
            col = riskCols, pch=21, title="Fire occurrence probability")
# End
}
