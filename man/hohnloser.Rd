\name{hohnloser}
\alias{hohnloser}
\title{
Left ventricular ejection fraction of patients with malignant ventricular
tachyarrhythmias.
}
\description{
A data frame with the left ventricular ejection fraction of patients with
malignant ventricular tachyarrhythmias including recurrence-free month and censoring.
}
\usage{
data("hohnloser")
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

The data was first published by Hohnloser et al. (1987), 
the data used here is published in Table 1 of Lausen and Schumacher (1992).

}
\references{

Hohnloser, S.H., Raeder, E.A., Podrid, P.J., Graboys, T.B. and Lown, B.
(1987), Predictors of antiarrhythmic drug efficacy in patients with
malignant ventricular tachyarrhythmias. \emph{American Heart Journal}
\bold{114}, 1--7


Lausen, B. and Schumacher, M. (1992), Maximally Selected Rank Statistics.
\emph{Biometrics} \bold{48}, 73--85
}

\examples{

set.seed(29)

library("survival")

# limiting distribution

maxstat.test(Surv(month, cens) ~ EF, data=hohnloser, 
smethod="LogRank", pmethod="Lau92")

# with integer valued scores for comparison

maxstat.test(Surv(month, cens) ~ EF, data=hohnloser, 
smethod="LogRank", pmethod="Lau92", iscores=TRUE)

# improved Bonferroni inequality

maxstat.test(Surv(month, cens) ~ EF, data=hohnloser,
smethod="LogRank", pmethod="Lau94")

maxstat.test(Surv(month, cens) ~ EF, data=hohnloser,
smethod="LogRank", pmethod="Lau94", iscores=TRUE)


# small sample solution by Hothorn & Lausen

maxstat.test(Surv(month, cens) ~ EF, data=hohnloser,
smethod="LogRank", pmethod="HL")

# normal approximation

maxstat.test(Surv(month, cens) ~ EF, data=hohnloser,
smethod="LogRank", pmethod="exactGauss")

maxstat.test(Surv(month, cens) ~ EF, data=hohnloser,
smethod="LogRank", pmethod="exactGauss", iscores=TRUE)

# conditional Monte-Carlo
maxstat.test(Surv(month, cens) ~ EF, data=hohnloser,
smethod="LogRank", pmethod="condMC", B = 9999)

}

\keyword{datasets}
