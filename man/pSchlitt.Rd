\name{pSchlitt}
\alias{pSchlitt}
\alias{qSchlitt}
\title{ Approximating Maximally Selected Statistics  }
\description{

Approximates the probability that a maximally selected rank statistic is
greater or equal to \code{b}.

}
\usage{
pSchlitt(b, N, m, method=c("genz", "drezner"))
qSchlitt(p, N, m)
}
\arguments{
  \item{b}{ quantile. }
  \item{p}{ probability. } 
  \item{N}{ number of observations. }
  \item{m}{ cutpoints. }
  \item{method}{ method of computing bivariate normal probabilities. Results
                 should be equal. }
}
\details{

Method based on a paper by Schlittgen. There are doubts if this works as
advertised. 

}
\value{

The probability that, under the hypothesis of independence, a maximally
selected statistic greater equal \code{b} is observed.

}
\references{ 
	
	Schlittgen, R. (1999).
	Regression Trees for Survival Data - an Approach to Select
	Discontinuous Split Points by Rank Statistics.
	\emph{Biometrical Journal}, \bold{41}(8), 943-954.
}
\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\examples{

# no example given in the paper!

pSchlitt(2.5, N = 200, m = seq(from=50, to=150, by=10))


}
\keyword{distribution}
