\name{pLausen94}
\alias{pLausen94}
\alias{qLausen94}
\title{ Approximating Maximally Selected Statistics }
\description{

Approximates the probability that a maximally selected rank statistic is
greater or equal to \code{b}.

}
\usage{
pLausen94(b, N, minprop=0.1, maxprop=0.9, m=NULL)
qLausen94(p, N, minprop=0.1, maxprop=0.9, m=NULL)
}
\arguments{
  \item{b}{ quantile. }
  \item{p}{ probability. }
  \item{N}{ number of observations.}
\item{minprop}{ at least \code{minprop}*100\% of the observations in the
first group. }
  \item{maxprop}{ not more than \code{minprop}*100\% of the observations in
the first group. }
  \item{m}{a integer vector containing the sample sizes in the first groups
for each cutpoint considered. If \code{is.null(m)} a continuous predictor is
assumed. }
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

        Lausen, B., Sauerbrei, W. & Schumacher, M. (1994).
        Classification and Regression Trees (CART) used for the
        exploration of prognostic factors measured on different scales.
        in: P. Dirschedl & R. Ostermann (Eds),
        \emph{Computational Statistics}, Heidelberg, Physica-Verlag,
        483--496

}

\examples{

p <- pLausen94(2.5, 20, 0.25, 0.75)

# Lausen 94, page 489

if (round(p, 3) != 0.073) stop("error checking pLausen94")

# the same

p2 <- pLausen94(2.5, 200, 0.25, 0.75, m=seq(from=50, to=150, by=10))

stopifnot(all.equal(round(p,3), round(p2,3)))

}

\keyword{distribution}
