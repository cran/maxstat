\name{pmaxstat}
\alias{pmaxstat}
\alias{qmaxstat}
\title{ Approximating Maximally Selected Statistics }
\description{

Approximates the probability that a maximally selected rank
statistic is greater or equal to \code{b}.

}
\usage{
pmaxstat(b, scores, msample, quant=FALSE)
qmaxstat(p, scores, msample)
}
\arguments{
 \item{b}{ quantile. }
 \item{p}{propability. }
 \item{scores}{ integer valued scores assigned to the observations. }
 \item{msample}{ all possible splitpoints. }
 \item{quant}{ logical. Returns the results of SR instead of P-values. Only
               to be used in \code{qmaxstat}.}
}
\details{

Small sample approximation by Hothorn and Lausen (2002).

}
\value{

An upper limit for the probability that, under the hypothesis of 
independence, a maximally selected statistic greater equal \code{b} 
is observed. \code{qmaxstat} needs optimization.

}

\references{

        Hothorn, T. and Lausen, B. (2002).
        On the Exact Distribution of Maximally Selected Rank Statistics.
        \emph{submitted}
}

\examples{

pmaxstat(2.5, 1:20, 5:15)

}
\keyword{distribution}
