\name{maxstat-internal}
\alias{cmaxstat}
\alias{hl}
\alias{index}
\alias{wh}
\title{Internal maxstat functions}
\description{
 Internal maxstat functions.
}
\usage{
cmaxstat(y, x=NULL, pmethod=c("none", "Lau92", "Lau94",
         "exactGauss", "HL", "min"), minprop = 0.1,
         maxprop=0.9, alpha = NULL, \dots)
hl(scores, H, E, S, msample, N, b)
index(x)
wh(cut, x)
}
\details{
  These functions are not to be called by the user.
}
\keyword{internal}
