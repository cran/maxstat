\name{cmatrix}
\alias{cmatrix}
\title{ Correlation Matrix }
\description{
 Correlation matrix of maximally selected rank statistics. 
}
\usage{
cmatrix(X)
}
\arguments{
  \item{X}{matrix of cutpoints}
}
\details{
  The correlations between all two-sample rank statistics induced by the
cutpoints in \code{X} are computed. 

}
\value{
The correlation matrix with dimension depending on ties in \code{X} is
returned.
}
\references{

        Hothorn, T. and Lausen, B. (2002).
        On the Exact Distribution of Maximally Selected Rank Statistics.
        \emph{Computational Statistics \& Data Analysis}, (in press).
        Preprint available from
        \url{http://www.mathpreprints.com/math/Preprint/hothorn/20020227/2/}

        Lausen, B., Hothorn, T., Bretz, F. and Schmacher, M. (2002).
        Assessment of Optimally Selected Prognostic Factors.
	\emph{submitted}. Preprint available from
	\url{http://www.mathpreprints.com/math/Preprint/blausen/20021007/1/}
}
 
\author{ Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> }
\examples{

X <- matrix(rnorm(30), ncol=3) 
print(system.time(a <- cmatrix(X)))	# this function

# coded by just typing the definition of the correlation

testcorr <- function(X) {
  a <- unlist(test <- apply(X, 2, index), recursive=FALSE)
  cnull <- rep(0, nrow(X))
  mycorr <- diag(length(a))
  for (i in 1:(length(a)-1)) {
    for (j in (i+1):length(a)) {
      cone <- cnull
      cone[a[[i]]] <- 1
      ctwo <- cnull
      ctwo[a[[j]]] <- 1
      sone <- sqrt(sum((cone - mean(cone))^2))
      stwo <- sqrt(sum((ctwo - mean(ctwo))^2))
      tcorr <- sum((cone - mean(cone))*(ctwo - mean(ctwo)))
      tcorr <- tcorr/(sone * stwo)
      mycorr[i,j] <- tcorr
    }
  }
  mycorr
}

print(system.time(tc <- testcorr(X)))
stopifnot(all.equal(tc, a))

}
\keyword{misc}
