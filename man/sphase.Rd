\name{sphase}
\alias{sphase}
\non_function{}
\title{ S-phase fraction of tumor cells }
\usage{data(sphase)}
\description{
  S-phase fraction of tumor cells in breast cancer patients. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{SPF}{S-phase fraction}
    \item{RFS}{recurrence free survival}
    \item{cens}{censoring indicator (\code{1} event)}
  }
}
\details{

}
\source{
  J. Pfisterer, F. Kommoss, W. Sauerbrei, D. Menzel, M. Kiechle, E. Giese,
M. Hilgarth and A. Pfleiderer (1995). DNA flow cytometry in node positive
breast cancer: Prognostic value and correlation to morphological and
clinical factors. \emph{Analytical and Quantitative Cytology and Histology}
bold{7}(6), 406--412.
}
\references{

}
\examples{
data(sphase)
maxstat.test(Surv(RFS, cens) ~ SPF, data=sphase, smethod="LogRank",
pmethod="Lau94")
maxstat.test(Surv(RFS, cens) ~ SPF, data=sphase, smethod="LogRank",
pmethod="Lau94", iscores=TRUE)
maxstat.test(Surv(RFS, cens) ~ SPF, data=sphase, smethod="LogRank",
pmethod="HL")
plot(maxstat.test(Surv(RFS, cens) ~ SPF, data=sphase, smethod="LogRank"))

}
\keyword{datasets}
