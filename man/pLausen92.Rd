\name{pLausen92}
\alias{pLausen92}
\alias{qLausen92}
\title{ Approximating Maximally Selected Statistics }
\description{

Approximates the probability that a maximally selected rank statistic is
greater or equal to \code{b}. 

}
\usage{
pLausen92(b, minprop=0.1, maxprop=0.9)
qLausen92(p, minprop=0.1, maxprop=0.9)
}
\arguments{
  \item{b}{ quantile. }
  \item{p}{ probability. }
  \item{minprop}{ at least \code{minprop}*100\% of the observations in the
first group. }
  \item{maxprop}{ not more than \code{minprop}*100\% of the observations in
the first group. }
}
\details{

Large sample approximation by Miller and Siegmund (1982) based on a Brownian
bridge, cf. Lausen and Schumacher (1992). 

}
\value{
The probability that, under the hypothesis of independence, a maximally
selected statistic greater equal \code{b} is observed.
}
\references{ 

	Miller, R. and Siegmund, D. (1982),
	Maximally Selected Chi Square Statistics.
	\emph{Biometrics}, \bold{38}, 1011--1016

	Lausen, B. and Schumacher, M. (1992),
	Maximally Selected Rank Statistics.
	\emph{Biometrics}, \bold{48}, 73--85

}

\examples{

# Compute quantiles. Should be equal to Table 2 in Lausen and Schumacher

load(file.path(system.file(package = "maxstat"), "results", "LausenTab2.rda"))

a <- rev(c(0.01, 0.025, 0.05, 0.1))
prop <- rbind(c(0.25, 0.75), c(0.4, 0.6), c(0.1, 0.9), c(0.4, 0.9))
Quant <- matrix(rep(0, length(a)*nrow(prop)), nrow=length(a)) 

for (i in 1:length(a)) {                                            
  for (j in 1:nrow(prop)) {                            
    Quant[i,j] <- qLausen92(a[i], minprop=prop[j,1], maxprop=prop[j,2]) 
  }
}

Quant <- round(Quant, 3)
rownames(Quant) <- a
colnames(Quant) <- c("A2575", "A46", "A19", "A49")
Quant <- as.data.frame(Quant)
rownames(LausenTab2) <- a

Quant

LausenTab2

if(!all.equal(LausenTab2, Quant)) stop("error checking pLausen92")

}

\keyword{distribution}
