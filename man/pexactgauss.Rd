\name{pexactgauss}
\alias{pexactgauss}
\alias{qexactgauss}
\title{ Computing Maximally Selected Gauss Statistics  }
\description{

Computes the exact probability that a maximally selected gauss statistic is
greater or equal to \code{b}.

}
\usage{
pexactgauss(b, N, m, \dots)
qexactgauss(p, N, m, \dots)
}
\arguments{
  \item{b}{ quantile. }
  \item{p}{ probability. } 
  \item{N}{ number of observations. }
  \item{m}{ cutpoints. }
  \item{\dots}{additional parameters to be passed to
\code{\link[mvtnorm]{pmvnorm}}.}
}
\details{

This is the exact distribution of a maximally selected Gauss statistic and
the asymptotic distribution for maximally selected rank statistics.
The correlation matrix can be computed using \code{\link{cmatrix}}, 
normal probabilities are derived from the procedures by Genz/Bretz (see
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

}
\examples{

pexact <- pexactgauss(2.5, 20, 2:18, abseps=0.01)

}
\keyword{distribution}
