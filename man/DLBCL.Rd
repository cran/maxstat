\name{DLBCL}
\alias{DLBCL}
\title{
Diffuse large B-cell lymphoma
}
\description{
A data frame with gene expression data from DLBCL (diffuse large B-cell
lymphoma) patients.
}
\usage{
data(DLBCL)
}
\format{
\describe{
\item{\code{DLCLid}}{DLCL identifier}
\item{\code{GEG}}{Gene Expression Group}
\item{\code{time}}{survival time in month}
\item{\code{cens}}{censoring: 0 cencored, 1 dead}
\item{\code{IPI}}{International Prognostic Index}
\item{\code{MGE}}{Mean Gene Expression}
}
Except of \code{MGE}, the data is published at
\url{http://llmpp.nih.gov/lymphoma/data.shtml}. 
}

\source{
Ash A. Alizadeh et. al (2000), Distinct types of diffuse large
B-cell lymphoma identified by gene
expression profiling. \emph{Nature}, \bold{403}, 504--509
}
\references{
Lausen, B. and Schumacher, M. (1992), Maximally Selected Rank Statistics.
\emph{Biometrics} \bold{48}, 73--85
}

\examples{

data(DLBCL)

# remove NA's

DLBCL <- DLBCL[!is.na(DLBCL$time),]

# compute the cutpoint


postscript("statDLBCL.ps",horizontal=F, width=8, height=8)
par(mai=c(1.0196235, 1.0196235, 0.8196973, 0.4198450))

mod <- maxstat.test(DLBCL$MGE, DLBCL$time,
             cens=DLBCL$cens, smethod="LogRank", plot=T, cex.lab=1.6,
             cex.axis=1.6, xlab="Mean gene expression")

dev.off()

# significance of the cutpoint
# Limiting distribution

maxstat.test(DLBCL$MGE, DLBCL$time,
             cens=DLBCL$cens, smethod="LogRank", pmethod="Lau92")

# improved Bonferroni inequality

maxstat.test(DLBCL$MGE, DLBCL$time,
             cens=DLBCL$cens, smethod="LogRank", pmethod="Lau94")

# small sample solution Hothorn & Lausen (2001)

maxstat.test(DLBCL$MGE, DLBCL$time,
             cens=DLBCL$cens, smethod="LogRank", pmethod="HL")

maxstat.test(DLBCL$MGE, DLBCL$time,
             cens=DLBCL$cens, smethod="LogRank", pmethod="exactGauss")

# Nature article survival analysis

splitGEG <- rep(1, nrow(DLBCL))
DLBCL <- cbind(DLBCL, splitGEG)
DLBCL$splitGEG[DLBCL$GEG == "Activated B-like"] <- 0

plot(survfit(Surv(time, cens) ~ splitGEG, data=DLBCL),
     xlab="Survival time in month", ylab="Probability")

text(90, 0.7, "GC B-like")
text(60, 0.3, "Activated B-like")

splitIPI <- rep(1, nrow(DLBCL))
DLBCL <- cbind(DLBCL, splitIPI)
DLBCL$splitIPI[DLBCL$IPI <= 2] <- 0

plot(survfit(Surv(time, cens) ~ splitIPI, data=DLBCL),
     xlab="Survival time in month", ylab="Probability")

text(90, 0.7, "Low clinical risk")
text(60, 0.25, "High clinical risk")

# survival analysis using the cutpoint 

splitMGE <- rep(1, nrow(DLBCL))
DLBCL <- cbind(DLBCL, splitMGE)
DLBCL$splitMGE[DLBCL$MGE <= mod$estimate] <- 0

postscript("survDLBCL.ps",horizontal=F, width=8, height=8)
par(mai=c(1.0196235, 1.0196235, 0.8196973, 0.4198450))

plot(survfit(Surv(time, cens) ~ splitMGE, data=DLBCL),
xlab = "Survival time in month",
ylab="Probability", cex.lab=1.6, cex.axis=1.6)

text(90, 0.9, expression("Mean gene expression" > 0.186), cex=1.6)   
text(90, 0.45, expression("Mean gene expression" <= 0.186 ), cex=1.6)   

dev.off()
}
\keyword{datasets}
