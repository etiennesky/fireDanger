\name{fmi}
\alias{fmi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of the Fuel Moisture Index
}
\description{
Computes the Fuel Moisture Index
}
\usage{
fmi(t, rh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
A numeric vector of temperature records, in degree C
}
  \item{rh}{
A numeric vector of relative humidity records, in \%
}
}
\value{
A numeric vector with the FMI values
}
\references{
Sharples, J.J., McRae, R.H.D., Weber, R.O., Gill, A.M.(2009) A simple index for assessing fuel moisture content. \emph{Environmental Modelling and Software}, 24, 637-646.

Sharples, J.J., McRae, R.H.D., Weber, R.O., Gill, A.M.(2009a) A simple index for assessing fire danger rating. \emph{Environmental Modelling and Software}. DOI:10.1016/j.envsoft.2008.11.004
}
\author{
J. Bedia
}
\note{
Fire moisture index is a dimensionless index and should not be considered as giving a direct estimate of fuel moisture content, as such. However, the question of principle interest was how FMI compared to the output of more complicated fuel moisture models. This question was answered in Sharples \emph{et al.} (2009a) where it was shown that FMI compared remarkably well to the output of existing models, especially the model for the moisture content of eucalypt litter
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{angstronIndex}}, for a similar code in the sense that both rate air rather than fuel moisture
}
\examples{
data(fwiSpain)
par(mfrow=c(4,2))
for (i in 1:length(fwiSpain$observed)) {
      x <- fwiSpain$observed[[i]]
      fmi.series <- fmi(t=x$T, rh=x$H)
      plot(x$fecha, fmi.series, ty='l', main = names(fwiSpain$observed[i]),
                  ylim = c(0,35))
}
# End
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ts }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
