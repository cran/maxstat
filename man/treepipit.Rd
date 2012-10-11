\name{treepipit}
\alias{treepipit}
\docType{data}
\title{ Tree Pipit Data }
\description{
  Counts of tree pipits at 86 raster points in oak forests.
}
\usage{data("treepipit")}
\format{
  A data frame with 86 observations on the following 2 variables.
  \describe{
    \item{counts}{number of tree pipits counted.}
    \item{coverstorey}{canopy overstorey in percent.}
  }
}
\details{
  The influence of canopy overstorey on the number of bird 
  individuals is of special interest.
}
\source{
  Data collected and kindly provided by Joerg Mueller
  <mue@lwf.uni-muenchen.de>.
}
\references{
  Mueller J. and Hothorn T. (2004), Maximally Selected Two-sample Statistics as a New 
  Tool for the Identification and Assessment of Habitat Factors with an Application to 
  Breeding-bird Communities in Oak Forests,
  \emph{European Journal of Forest Research}, \bold{123}(3), 219--228.
}
\examples{

mod <- maxstat.test(counts ~ coverstorey, data = treepipit, 
                    smethod = "Data", pmethod = "HL", minprop = 0.2,
                    maxprop = 0.8)
print(mod)
plot(mod)
}
\keyword{datasets}
