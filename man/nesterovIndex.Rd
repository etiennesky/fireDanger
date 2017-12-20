\name{nesterovIndex}
\alias{nesterovIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the Nesterov Index and the Modified Nesterov Index 
}
\description{
Computes the Nesterov Index (NI) and the Modified Nesterov Index (MNI). NI was developed in former Soviet Union as an empirical function reflecting the relationship between observed weather conditions and fire occurrence.
}
\usage{
nesterovIndex(t, rh, p, modified = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{A numeric vector of temperature records, in degree C}
  \item{rh}{A numeric vector of relative humidity records, in \%}
  \item{p}{A numeric vector of precipitation records, in mm}
  \item{modified}{Logical flag indicating if the original/modified (\code{FALSE/TRUE}) version of the index is computed. Default to original (\code{FALSE}).
}
}
\details{
NI is a cumulative index, but summation is performed for those days when the daily precipitation (p) does not exceed 3 mm. At p > 3 mm, the NI value is reset zero. Usually, the values from NI are divided into five ranges to provide an estimate of fire danger potential. Conditions with NI < 300 (regime I) are not considered hazardous. Conditions in the ranges 300-1000, 1000-4000, 4000-10000, and > 10000 are considered regimes with low (II), moderate (III), high (IV), and extreme (V) level of fire hazard.

Previous studies reveal that NI may be unstable in some cases, and a modification of this index has been proposed by introducing to its values a K scale coefficient (see Groisman \emph{et al.} 2007) in the range [0,1], accounting for the amount of precipitation and previous dryness in a more detailed way than the original NI. This alternative version is obtained by setting the argument \code{modified = TRUE}.  
}
\value{
A numeric vector with the index values
}
\references{
Groisman, P.Y., Sherstyukov, B.G., Razuvaev, V.N., Knight, R.W., Enloe, J.G., Stroumentova, N.S., Whitfield, P.H., Forland, E., Hannsen-Bauer, I., Tuomenvirta, H., Aleksandersson, H., Mescherskaya, A.V., Karl, T.R., 2007. Potential forest fire danger over Northern Eurasia: Changes during the 20th century. \emph{Global and Planetary Change} 56, 371-386
}
\author{
J. Bedia
}

\seealso{
\code{\link{kbdi}} and some components of \code{\link{fwi}} (DC) for other cumulative indices.
}
\examples{
data(fwiSpain)

# Difference between NI and MNI
par(mfrow = c(4,2))
for (i in 1:length(fwiSpain$observed)) {
      x <- fwiSpain$observed[[i]]
      ni <- nesterovIndex(t=x$T, rh=x$H, p=x$P)
      mni <- nesterovIndex(t=x$T, rh=x$H, p=x$P, modified = TRUE)
      plot(x$fecha, ni - mni, ty = "l", main = names(fwiSpain$observed)[i])
}
#End
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
