\name{pexactgauss}
\alias{pexactgauss}
\alias{qexactgauss}
\title{ Computing Maximally Selected Gauss Statistics  }
\description{

Computes the exact probability that a maximally selected gauss statistic is
greater or equal to \code{b}.

}
\usage{
pexactgauss(b, N, m, maxpts)
qexactgauss(p, N, m)
}
\arguments{
  \item{b}{ quantile. }
  \item{p}{ probability. } 
  \item{N}{ number of observations. }
  \item{m}{ cutpoints. }
  \item{maxpts}{maximum number of function values as integer}
}
\details{

This is the exact distribution of a maximally selected Gauss statistic and
the asymptotic distribution for maximally selected rank statistics.
Correlations from Schlittgen (1999), normal probabilities by Genz/Bretz (see
\code{\link[mvtnorm]{pmvnorm}} for details).

}
\value{

The probability that, under the hypothesis of independence, a maximally
selected gauss statistic greater equal \code{b} is observed. 

}
\references{ 
	
        Genz, A. (1992). Numerical computation of multivariate normal probabilities.
        \emph{Journal of Computational and Graphical Statistics}, \bold{1}, 141--150

        Genz, A. (1993). Comparison of methods for the computation of multivariate
        normal probabilities. \emph{Computing Science and Statistics}, \bold{25},
        400--405

	Schlittgen, R. (1999),
	Regression Trees for Survival Data - an Approach to Select
	Discontinuous Split Points by Rank Statistics.
	\emph{Biometrical Journal}, \bold{41}(8), 943--954.
}

\examples{

pexact <- pexactgauss(2.5, 20, 2:18)

}
\keyword{distribution}
