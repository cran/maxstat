\name{pLausen94}
\alias{pLausen94}
\alias{qLausen94}
\title{ Approximating Maximally Selected Statistics }
\description{

Approximates the probability that a maximally selected rank statistic is
greater or equal to \code{b}.

}
\usage{
pLausen94(b, N, minprop=0.1, maxprop=0.9)
qLausen94(p, N, minprop=0.1, maxprop=0.9)
}
\arguments{
  \item{b}{ quantile. }
  \item{p}{ probability. }
  \item{N}{ number of cutpoints.}
\item{minprop}{ at least \code{minprop}*100\% of the observations in the
first group. }
  \item{maxprop}{ not more than \code{minprop}*100\% of the observations in
the first group. }
}
\details{

Approximation based on an improved Bonferroni inequality.

}
\value{
The probability that, under the hypothesis of independence, a maximally
selected statistic greater equal \code{b} is observed.
}
\references{ 

	Worsley, K.J. (1982),
	An Improved Bonferroni Inequality and Applications.
	\emph{Biometrika}, \bold{69}, 297--302

	Lausen, B. (1990),
	Maximal Selektierte Rangstatistiken. Dissertation. 
	Universit\"at Dortmund 

}

\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\examples{

p <- pLausen94(2.5, 20, 0.25, 0.75)

if (round(p, 3) != 0.073) stop("error checking pLausen94")

}

\keyword{distribution}
