\name{pexactgauss}
\alias{pexactgauss}
\alias{qexactgauss}
\title{ Computing Maximally Selected Gauss Statistics  }
\description{

Computes the exact probability that a maximally selected gauss statistic is
greater or equal to \code{b}.

}
\usage{
pexactgauss(b, x, minprop=0.1, maxprop=0.9, \dots)
qexactgauss(p, x, minprop=0.1, maxprop=0.9, \dots)
}
\arguments{
  \item{b}{ quantile. }
  \item{p}{ probability. } 
  \item{x}{ the prognostic factor(s) under test. }
  \item{minprop}{ at least \code{minprop}*100\% of the observations in the
                  first group. }
  \item{maxprop}{ not more than \code{minprop}*100\% of the observations in
                  the first group. }
  \item{\dots}{additional parameters to be passed to
\code{\link[mvtnorm]{pmvnorm}}.}
}
\details{

This is the exact distribution of a maximally selected Gauss statistic and
the asymptotic distribution for maximally selected rank statistics.
Normal probabilities are derived from the procedures by Genz/Bretz (see
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

x <- rnorm(20)

pexact <- pexactgauss(2.5, x, abseps=0.01)

}
\keyword{distribution}
