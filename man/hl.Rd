\name{hl}
\alias{hl}
\title{ Approximation by Hothorn and Lausen (2001) }
\description{
 Technical parts
}
\usage{
hl(H, E, S, msample, N, b)
}
\arguments{
  \item{H}{ SR-matrix }
  \item{E}{ matrix of expections }
  \item{S}{ matrix of variances }
  \item{msample}{ possible cutpoints }
  \item{N}{ sample size }
  \item{b}{ critical value }
}
\details{
  see \code{\link{qmaxstat}}.
}
\value{
A list wih components
  \item{lower}{lower bound}
  \item{upper}{upper bound}
}
\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\keyword{distribution}
