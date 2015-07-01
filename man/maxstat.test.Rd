\name{maxstat.test}
\alias{maxstat.test}
\alias{maxstat.test.data.frame}
\alias{maxstat.test.default}
\alias{maxstat}
\title{ Maximally Selected Rank and Statistics }
\description{
Performs a test of independence of a response and one or more covariables
using maximally selected rank statistics.
}
\usage{
\method{maxstat.test}{data.frame}(formula, data, subset, na.action, \dots)
maxstat(y, x=NULL, weights = NULL, smethod=c("Wilcoxon", "Median",
            "NormalQuantil","LogRank", "Data"), pmethod=c("none", "Lau92",
            "Lau94", "exactGauss", "HL", "condMC", "min"), iscores=(pmethod=="HL"),
            minprop = 0.1, maxprop=0.9, alpha = NULL, keepxy=TRUE, \dots)
}
\arguments{
  \item{y}{ numeric vector of data values, dependent variable. }
  \item{x}{ numeric vector of data values, independent variable. }
  \item{weights}{ an optional numeric vector of non-negative weights, 
                  summing to the number of observations. }
  \item{smethod}{ kind of statistic to be computed, i.e. defines the scores
to be used for computing the statistic. }
  \item{pmethod}{ kind of p-value approximation to be used. }
  \item{iscores}{ logical: should the scores be mapped into integers
\code{1:length(x)}? This is TRUE by default for \code{pmethod=="HL"} and
FALSE otherwise.}
  \item{minprop}{ at least \code{minprop}*100\% of the observations in the
                  first group. }
  \item{maxprop}{ not more than \code{minprop}*100\% of the observations in
                  the first group. }
  \item{alpha}{significance niveau, the appropriate quantile is computed if
   \code{alpha} is specified. Used for plotting within \code{\link{plot.maxtest}}.}
  \item{keepxy}{logical: return \code{y} and \code{x} as elements of the
\code{maxtest} object.}
  \item{formula}{ a formula describing the model to be tested of the form
    \code{lhs ~ rhs} where \code{lhs} is the response variable. 
    For survival problems, i.e. using the log-rank statistic, the
    formula is of the form \code{Surv(time, event) ~ rhs}, see above.}
  \item{data}{an data frame containing the variables in the
    model formula. \code{data} is required.}
  \item{subset}{an optional vector specifying a subset of observations
    to be used.}
  \item{na.action}{a function which indicates what should happen when
    the data contain \code{NA}s.  Defaults to
    \code{getOption("na.action")}.}
  \item{\dots}{additional parameters to be passed to
               \code{\link[mvtnorm]{pmvnorm}} or \code{B}, an integer
               defining the number of Monte-Carlo replications.}
}
\details{

The assessment of the predictive power of a variable \code{x} for a
dependent variable \code{y} can be determined by a maximally selected rank
statistic. 

\code{smethod} determines the kind of statistic to be used. 
\code{Wilcoxon} and \code{Median} denote maximally selected 
Wilcoxon and Median statistics. \code{NormalQuantile} and 
\code{LogRank} denote v.d. Waerden and log-rank
scores. 

\code{pmethod} specifies which kind of approximation of the p-value should
be used. \code{Lau92} is the limiting distribution by a Brownian bridge
(see Lausen and Schumacher, 1992, and \code{\link{pLausen92}}), 
\code{Lau94} the approximation based on an improved Bonferroni 
inequality (see Lausen, Sauerbrei and Schumacher, 1994, and \code{\link{pLausen94}}).

\code{exactGauss} returns the exact p-value for a maximally selected Gauss
statistic, see Hothorn and Lausen (2003).

\code{HL} is a small sample approximation based on the Streitberg-R\"ohmel
algorithm (see \code{\link[exactRankTests]{pperm}}) introduced by Hothorn and
Lausen (2003). This requires integer
valued scores. For v. d. Waerden and Log-rank scores we try to find
integer valued scores having the same shape. This results in slightly
different scores (and a different test), the procedure is described in
Hothorn (2001) and Hothorn and Lausen (2003).  

All the approximations are known to be conservative, so \code{min} gives
the minimum p-value of all procedures.

\code{condMC} simulates the distribution via conditional Monte-Carlo. 

For survival problems, i.e. using a maximally selected log-rank statistic,
the interface is similar to \code{\link[survival]{survfit}}. The depended
variable is a survival object \code{Surv(time, event)}. The argument 
\code{event} may be a numeric vector of \code{0} (alive) and \code{1}
(dead) or a vector of logicals with \code{TRUE} indicating death.

If more than one covariable is specified in the right hand side of
\code{formula} (or if \code{x} is a matrix or data frame), the variable with
smallest p-value is selected and the p-value for the global test problem of
independence of \code{y} and every variable on the right hand side is
returned (see Lausen et al., 2002). 

}
\value{
An object of class \code{maxtest} or \code{mmaxtest} (if more than one
covariable was specified) containing the following components
is returned:
  \item{statistic }{the value of the test statistic.}
  \item{p.value }{the p-value for the test.}
  \item{smethod}{ the type of test applied.}
  \item{pmethod}{ the type of p-value approximation applied.}
  \item{estimate}{ the estimated cutpoint (of \code{x}) which separates
                   \code{y} best. For numeric data, the groups are defined
                   by \code{x} less or equal to \code{estimate} and \code{x} greater \code{estimate}.}
  \item{maxstats}{a list of \code{maxtest} objects, one for each
                  covariable.}
  \item{whichmin}{an integer specifying the element of \code{maxstats} with
                  smallest p-value.}
  \item{p.value}{the p-value of the global test.}
  \item{univp.values}{the p-values for each of the variables under test.}
  \item{cm}{the correlation matrix the p-value is based on.}
\code{\link{plot.maxtest}} and \code{\link{print.maxtest}} can be used for
plotting and printing. If \code{keepxy = TRUE}, there are elements \code{y}
and \code{x} giving the response and independent variable.
}
\references{ 

        Hothorn, T. and Lausen, B. (2003).
        On the Exact Distribution of Maximally Selected Rank Statistics. 
        \emph{Computational Statistics & Data Analysis}, \bold{43}, 
        121--137.

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

        Lausen, B., Hothorn, T., Bretz, F. and Schmacher, M. (2004).
        Assessment of Optimally Selected Prognostic Factors. \emph{Biometrical Journal},
        \bold{46}(3), 364--374.

}

\examples{

set.seed(29)

x <- sort(runif(20))
y <- c(rnorm(10), rnorm(10, 2))
mydata <- data.frame(cbind(x,y))

mod <- maxstat.test(y ~ x, data=mydata, smethod="Wilcoxon", pmethod="HL",
                    minprop=0.25, maxprop=0.75, alpha=0.05)
print(mod)
plot(mod)

# adjusted for more than one prognostic factor.
library("survival")
mstat <- maxstat.test(Surv(time, cens) ~ IPI + MGE, data=DLBCL, 
                      smethod="LogRank", pmethod="exactGauss", 
                      abseps=0.01)
plot(mstat)

### sphase
set.seed(29)
data("sphase", package = "TH.data")

maxstat.test(Surv(RFS, event) ~ SPF, data=sphase, smethod="LogRank",
             pmethod="Lau94")

maxstat.test(Surv(RFS, event) ~ SPF, data=sphase, smethod="LogRank",
             pmethod="Lau94", iscores=TRUE)

maxstat.test(Surv(RFS, event) ~ SPF, data=sphase, smethod="LogRank",
             pmethod="HL")

maxstat.test(Surv(RFS, event) ~ SPF, data=sphase, smethod="LogRank",
             pmethod="condMC", B = 9999)

plot(maxstat.test(Surv(RFS, event) ~ SPF, data=sphase, smethod="LogRank"))


}
\keyword{htest}
