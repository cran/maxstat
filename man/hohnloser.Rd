\name{hohnloser}
\alias{hohnloser}
\title{
Left ventricular ejection fraction of patients with malignant ventricular
tachyarrhythmias
}
\description{
A data frame with the left ventricular ejection fraction of patients with
malignant ventricular including recurrence-free month and censoring.
}
\usage{
data(hohnloser)
}
\format{
\describe{
\item{\code{EF}}{left ventricular ejection in percent}
\item{\code{month}}{recurrence-free month}
\item{\code{cens}}{censoring: 0 cencored, 1 not censored}
}
The data used here is published in Table 1 of Lausen and Schumacher (1992).
}
\source{
Hohnloser, S.H., Raeder, E.A., Podrid, P.J., Graboys, T.B. and Lown, B.
(1987), Predictors of antiarrhythmic drug efficacy in patients with
malignant ventricular tachyarrhythmias. \emph{American Heart Journal}
\bold{114}, 1--7
}
\references{
Lausen, B. and Schumacher, M. (1992), Maximally Selected Rank Statistics.
\emph{Biometrics} \bold{48}, 73--85
}

\examples{

data(hohnloser)

# The approximation used in Lausen & Schumacher

maxstat.test(hohnloser$EF, hohnloser$month, cens=hohnloser$cens,
smethod="LogRank", pmethod="Lau92")

# improved Bonferroni inequality

maxstat.test(hohnloser$EF, hohnloser$month, cens=hohnloser$cens,
smethod="LogRank", pmethod="Lau94")

# The approximation by Hothorn & Lausen

maxstat.test(hohnloser$EF, hohnloser$month, cens=hohnloser$cens,
smethod="LogRank", pmethod="HL", plot=T)

maxstat.test(hohnloser$EF, hohnloser$month, cens=hohnloser$cens,
smethod="LogRank", pmethod="exactGauss")

}

\keyword{datasets}
