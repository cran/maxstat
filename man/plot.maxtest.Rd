\name{plot.maxtest}
\alias{plot.maxtest}
\alias{print.maxtest}
\alias{plot.mmaxtest}
\alias{print.mmaxtest}
\title{ Print and Plot Standardized Statistics }
\description{
Printing and ploting method of objects of class \code{maxtest}}
\usage{
\method{plot}{maxtest}(x, xlab=NULL, ylab=NULL, \dots)
\method{print}{maxtest}(x, digits = getOption("digits"), \dots)
\method{plot}{mmaxtest}(x, xlab=NULL, ylab=NULL, nrow=2, \dots)
\method{print}{mmaxtest}(x, digits = getOption("digits"), \dots)
}
\arguments{
  \item{x}{ an object of class \code{maxtest} or \code{mmaxtest}.}
  \item{xlab}{ label of x-axis. }
  \item{ylab}{ label of y-axis. }
  \item{nrow}{ number of rows for multiple plots at one page. }
  \item{digits}{ number of significant digits to be printed. }
  \item{\dots}{ additional arguments to \code{plot} or \code{print.htest}.}
}
\details{

The standardized statistics are plotted. If \code{alpha} was given in 
\code{\link{maxstat.test}} the appropriate significance bound is plotted as
a red line. \code{print.maxtest} is just a wrapper to \code{print.htest}.

}

\examples{

set.seed(29)

x <- sort(runif(20))
y <- rbinom(20, 1, 0.5)
mydata <- data.frame(c(x,y))

mod <- maxstat.test(y ~ x, data=mydata, smethod="Median", 
                    pmethod="HL", alpha=0.05)
print(mod)
plot(mod)

}
\keyword{htest}
