\name{treepipit}
\alias{treepipit}
\docType{data}
\title{ Tree Pipit Data }
\description{
  Counts of tree pipits at 86 raster points in oak forests.
}
\usage{data(treepipit)}
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
  Mueller J. and Hothorn T. (2003), On the Identification and Assessment of 
  Habitat Patterns with Impact in Breeding Bird Communities in Oak Forests.
  \emph{submitted manuscript}.
}
\examples{
data(treepipit)
mod <- maxstat.test(counts ~ coverstorey, data = treepipit, 
                    smethod = "Data", pmethod = "HL", minprop = 0.2,
                    maxprop = 0.8)
print(mod)
plot(mod)
}
\keyword{datasets}
