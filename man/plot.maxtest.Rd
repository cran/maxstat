\name{plot.maxtest}
\alias{plot.maxtest}
\alias{print.maxtest}
\title{ Print and Plot Standardized Statistics }
\description{
Printing and ploting method of objects of class \code{maxtest}}
\usage{
\method{plot}{maxtest}(x, xlab=NULL, ylab=NULL, \dots)
\method{print}{maxtest}{x, \dots}
}
\arguments{
  \item{x}{ an object of class \code{maxtest}.}
  \item{xlab}{ label of x-axis. }
  \item{ylab}{ label of y-axis. }
  \item{\dots}{ additional arguments to \code{plot} or \code{print.htest}.}
}
\details{

The standardized statistics are plotted. If \code{alpha} was given in 
\code{\link{maxstat.test}} the appropriate significance bound is plotted as
a red line. \code{print.maxtest} is just a wrapper to \code{print.htest}.

}

\examples{

x <- sort(runif(20))
y <- rbinom(20, 1, 0.5)

mod <- maxstat.test(y ~ x, smethod="Median", pmethod="HL", alpha=0.05)
print(mod)
plot(mod)

}
\keyword{htest}
