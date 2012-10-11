\name{corrmsrs}
\alias{corrmsrs}
\title{ Correlation Matrix }
\description{
 Correlation matrix of maximally selected rank statistics. 
}
\usage{
corrmsrs(X, minprop=0.1, maxprop=0.9)
}
\arguments{
  \item{X}{the vector, matrix or data.frame of prognostic factors under
test.}
  \item{minprop}{ at least \code{minprop}*100\% of the observations in the
                  first group. }
  \item{maxprop}{ not more than \code{minprop}*100\% of the observations in
                  the first group. }

}
\details{
  The correlations between all two-sample rank statistics induced by all
possible cutpoints in \code{X} are computed. 

}
\value{
The correlation matrix with dimension depending on ties in \code{X} is
returned.
}
\references{

        Hothorn, T. and Lausen, B. (2003).
        On the Exact Distribution of Maximally Selected Rank Statistics. 
        \emph{Computational Statistics & Data Analysis}, \bold{43}, 
        121--137.

        Lausen, B., Hothorn, T., Bretz, F. and Schmacher, M. (2004).
        Assessment of Optimally Selected Prognostic Factors. \emph{Biometrical Journal},
        \bold{46}(3), 364--374.
}
 
\examples{

set.seed(29)

# matrix of hypothetical prognostic factors

X <- matrix(rnorm(30), ncol=3) 

# this function

print(system.time(a <- corrmsrs(X, minprop=0, maxprop=0.999)))

# coded by just typing the definition of the correlation

testcorr <- function(X) {
  wh <- function(cut, x)
    which(x <= cut)
  index <- function(x) {
    ux <- unique(x)
    ux <- ux[ux < max(ux)]
    lapply(ux, wh, x = x)
  }
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
tc <- tc + t(tc)
diag(tc) <- 1
stopifnot(all.equal(tc, a))

}
\keyword{misc}
