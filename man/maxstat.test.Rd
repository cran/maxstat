\name{maxstat.test}
\alias{maxstat.test}
\title{ Maximally Selected Rank and Gauss Statistics }
\description{

Performs a maximally selected rank or Gauss test on bivariate observations. 

}
\usage{
maxstat.test(x, y, cens = NULL, smethod=c("Gauss", "Wilcoxon", "Median",
            "NormalQuantil","LogRank"), pmethod=c("none", "Lau92",
            "Lau94", "exactGauss", "HL", "min"),
            minprop = 0.1, maxprop=0.9, plot=F, xlab=NULL,...)
}
\arguments{
  \item{x}{ vector of data values, independent variable. }
  \item{y}{ numeric vector of data values, dependent variable. }
  \item{cens}{censoring: 0 cencored, 1 dead (only meaningful with LogRank)}
  \item{smethod}{ kind of statistic to be computed. }
  \item{pmethod}{ kind of p-value approximation to be used. }
  \item{minprop}{ at least \code{minprop}*100\% of the observations in the
                  first group. }
  \item{maxprop}{ not more than \code{minprop}*100\% of the observations in
                  the first group. }
  \item{plot}{ logical. Should the process be plotted? }
  \item{xlab}{use alternative \code{xlab}}
  \item{...}{additional plotting parameters}
}
\details{

The assessment of the predictive power of a variable \code{x} for a
dependent variable \code{y} can be determined by a maximally selected rank
or Gauss test. 

\code{smethod} determines the kind of statistic to be used. \code{Gauss}
means a maximally selected Gauss statistic, \code{Wilcoxon} and
\code{Median} simply denote Wilcoxon and Median scores.  
\code{NormalQuantile} and \code{LogRank} denote v.d. Waerden and Log-rank
scores. 

\code{pmethod} specifies which kind of approximation of the p-value should
be used. \code{Lau92} means the limiting distribution by a Brownian bridge
(see \code{\link{pLausen92}}), \code{Lau94} the approximation based on an
improved Bonferroni inequality (see \code{\link{pLausen94}}).

\code{exactGauss} returns the exact p-value for a maximally selected Gauss
statistic.

\code{HL} is a small sample approximation based on the Streitberg-R\"ohmel
algorithm (see \code{\link[exactRankTest]{pperm}}) introduced by Hothorn &
Lausen (2001). This requires integer
valued scores. For v. d. Waerden and Log-rank scores we try to find
integer valued scores having the same shape. This results in slightly
different scores (and a different test), the procedure is described in
Hothorn (2001).  

All the approximations are known to be conservative, so \code{min} gives
the minimum p-value of all procedures.

}
\value{
A list of class \code{htest} containing the following components:
  \item{statistic }{the value of the test statistic.}
  \item{p.value }{the p-value for the test.}
  \item{method}{ the type of test and p-value approximation applied.}
  \item{estimate}{ the estimated cutpoint (of \code{x}) which separates
                   \code{y} best.} 
}
\references{ 

        Hothorn, T. and Lausen, B. (2001).
        On the Exact Distribution of Maximally Selected Rank Statistics. 
        \emph{submitted}

	Lausen, B. and Schumacher, M. (1992).
        Maximally Selected Rank Statistics.
        \emph{Biometrics}, \bold{48}, 73--85

	Hothorn, T. (2001). On Exact Rank Tests in R.
	\emph{R News}, \bold{1}, 11--12

}

\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }

\examples{

x <- sort(runif(20))
y <- c(rnorm(10), rnorm(10, 2))

maxstat.test(x, y, smethod="Wilcoxon", pmethod="HL", minprop=0.25, maxprop=0.75)

}
\keyword{htest}
