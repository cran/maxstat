\name{maxstat-internal}
\alias{intrank}
\alias{p2norm}
\alias{p2normG}
\alias{pbvmax}
\alias{msurv}
\alias{hl}
\alias{pSchlitt}
\alias{qSchlitt}
\alias{corlm}
\alias{wh}
\alias{index}
\alias{corrmatrix}
\alias{forsel}

\title{Internal maxstat functions}
\description{
 Internal maxstat functions.
}
\usage{
intrank(x)
p2norm(h, k, rho, maxp=3000)
p2normG(h, k, rho, maxp=3000)
pbvmax(cor, t, method=c("drezner", "genz"))
msurv(time, event)
hl(scores, H, E, S, msample, N, b)
pSchlitt(b, N, m, method=c("genz", "drezner"))
qSchlitt(p, N, m)
corlm(Xone, Xtwo)
wh(cut, x)
index(x)
corrmatrix(X)
forsel(y, Xin, Xout, ret=FALSE)

}
\details{
  These functions are not to be called by the user.
}
\keyword{internal}
