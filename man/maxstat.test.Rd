\name{maxstat.test}
\alias{maxstat.test}
\alias{maxstat.test.formula}
\alias{maxstat.test.default}
\title{ Maximally Selected Rank and Gauss Statistics }
\description{
Performs a maximally selected rank or Gauss test on bivariate observations. 
}
\usage{
\method{maxstat.test}{default}(x, y, event = NULL, smethod=c("Gauss", "Wilcoxon", "Median",
            "NormalQuantil","LogRank"), pmethod=c("none", "Lau92",
            "Lau94", "exactGauss", "HL", "min"),
            minprop = 0.1, maxprop=0.9, alpha = NULL, \dots)
\method{maxstat.test}{formula}(formula, data, subset, na.action, \dots)
}
\arguments{
  \item{x}{ numeric vector of data values, independent variable. }
  \item{y}{ numeric vector of data values, dependent variable. }
  \item{event}{censoring: \code{0=alive}, \code{1=dead} (only meaningful with LogRank)}
  \item{smethod}{ kind of statistic to be computed. }
  \item{pmethod}{ kind of p-value approximation to be used. }
  \item{minprop}{ at least \code{minprop}*100\% of the observations in the
                  first group. }
  \item{maxprop}{ not more than \code{minprop}*100\% of the observations in
                  the first group. }
  \item{alpha}{significance niveau, the appropriate quantile is computed if
   \code{alpha} is specified. Used for plotting within \code{\link{plot.maxtest}}.}
  \item{formula}{ a formula describing the model to be tested of the form
    \code{y ~ x} where \code{y} is the dependent and \code{x} the independent
    variable. For survival problems, i.e. using the log-rank statistic, the
formula is of the form \code{Surv(time, event) ~ x}, see above.}
  \item{data}{an optional data frame containing the variables in the
    model formula.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to
    \code{getOption("na.action")}.}
  \item{\dots}{additional parameters to methods.}
}
\details{

The assessment of the predictive power of a variable \code{x} for a
dependent variable \code{y} can be determined by a maximally selected rank
or Gauss test. 

\code{smethod} determines the kind of statistic to be used. 
\code{Wilcoxon} and \code{Median} denote maximally selected 
Wilcoxon and Median statistics. \code{NormalQuantile} and 
\code{LogRank} denote v.d. Waerden and log-rank
scores. \code{Gauss} means a maximally selected Gauss statistic.

\code{pmethod} specifies which kind of approximation of the p-value should
be used. \code{Lau92} is the limiting distribution by a Brownian bridge
(see Lausen and Schumacher, 1992, and \code{\link{pLausen92}}), 
\code{Lau94} the approximation based on an improved Bonferroni 
inequality (see Lausen, Sauerbrei and Schumacher, 1994, and \code{\link{pLausen94}}).

\code{exactGauss} returns the exact p-value for a maximally selected Gauss
statistic, see Hothorn and Lausen (2002).

\code{HL} is a small sample approximation based on the Streitberg-R\"ohmel
algorithm (see \code{\link[exactRankTest]{pperm}}) introduced by Hothorn and
Lausen (2002). This requires integer
valued scores. For v. d. Waerden and Log-rank scores we try to find
integer valued scores having the same shape. This results in slightly
different scores (and a different test), the procedure is described in
Hothorn (2001) and Hothorn and Lausen (2002).  

All the approximations are known to be conservative, so \code{min} gives
the minimum p-value of all procedures.

For survival problems, i.e. using a maximally selected log-rank statistic,
the interface is similar to \code{\link[survival]{survfit}}. The depended
variable is a survival object \code{Surv(time, event)} (note: in contrast to
\code{Surv} in package \code{survival} only \code{time} and \code{event}
can be specified). The argument 
\code{event} may be a numeric vector of \code{0} (alive) and \code{1}
(dead) or a vector of logicals with \code{TRUE} indicating death.

}
\value{
A list of class \code{maxtest} containing the following components
is returned:
  \item{statistic }{the value of the test statistic.}
  \item{p.value }{the p-value for the test.}
  \item{method}{ the type of test and p-value approximation applied.}
  \item{estimate}{ the estimated cutpoint (of \code{x}) which separates
                   \code{y} best.} 
\code{\link{plot.maxtest}} and \code{\link{print.maxtest}} can be used for
plotting and printing.
}
\references{ 

        Hothorn, T. and Lausen, B. (2002).
        On the Exact Distribution of Maximally Selected Rank Statistics. 
        \emph{submitted}, preprint available from
        \url{http://www.mathpreprints.com/math/Preprint/hothorn/20020227/2}

	Lausen, B. and Schumacher, M. (1992).
        Maximally Selected Rank Statistics.
        \emph{Biometrics}, \bold{48}, 73--85

        Lausen, B., Sauerbrei, W. and Schumacher, M. (1994).
        Classification and Regression Trees (CART) used for the
        exploration of prognostic factors measured on different scales.
        in: P. Dirschedl and R. Ostermann (Eds), 
        \emph{Computational Statistics}, Heidelberg, Physica-Verlag,
        483--496

	Hothorn, T. (2001). On Exact Rank Tests in R.
	\emph{R News}, \bold{1}, 11--12

}

\examples{

x <- sort(runif(20))
y <- c(rnorm(10), rnorm(10, 2))

mod <- maxstat.test(y ~ x, smethod="Wilcoxon", pmethod="HL",
                    minprop=0.25, maxprop=0.75, alpha=0.05)
print(mod)
plot(mod)

}
\keyword{htest}
