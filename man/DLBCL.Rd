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
\item{\code{DLCLid}}{DLBCL identifier}
\item{\code{GEG}}{Gene Expression Group}
\item{\code{time}}{survival time in month}
\item{\code{cens}}{censoring: 0 cencored, 1 dead}
\item{\code{IPI}}{International Prognostic Index}
\item{\code{MGE}}{Mean Gene Expression}
}
}
\source{

Except of \code{MGE}, the data is published at
\url{http://llmpp.nih.gov/lymphoma/data.shtml}. \code{MGE} is the mean of
the gene expression.

}
\references{

Ash A. Alizadeh et. al (2000), Distinct types of diffuse large
B-cell lymphoma identified by gene
expression profiling. \emph{Nature}, \bold{403}, 504--509

}

\examples{

data(DLBCL)

# compute the cutpoint and plot the empirical process 

mod <- maxstat.test(Surv(time, cens) ~ MGE, data=DLBCL, smethod="LogRank")

print(mod)

\dontrun{
  # postscript("statDLBCL.ps", horizontal=F, width=8, height=8)
  pdf("statDLBCL.pdf", width=8, height=8)
}
par(mai=c(1.0196235, 1.0196235, 0.8196973, 0.4198450))
plot(mod, cex.lab=1.6, cex.axis=1.6, xlab="Mean gene expression",lwd=2)
\dontrun{
  dev.off()
}

# significance of the cutpoint
# limiting distribution

maxstat.test(Surv(time, cens) ~ MGE, data=DLBCL,
             smethod="LogRank", pmethod="Lau92", iscores=TRUE)

# improved Bonferroni inequality, plot with significance bound

maxstat.test(Surv(time, cens) ~ MGE, data=DLBCL,
             smethod="LogRank", pmethod="Lau94", iscores=TRUE)

mod <- maxstat.test(Surv(time, cens) ~ MGE, data=DLBCL, smethod="LogRank",
                    pmethod="Lau94", alpha=0.05)
plot(mod, xlab="Mean gene expression")

\dontrun{
#  postscript(file="RNewsStat.ps",horizontal=F, width=8, height=8)
   pdf("RNewsStat.pdf", width=8, height=8)

}
par(mai=c(1.0196235, 1.0196235, 0.8196973, 0.4198450))
plot(mod, xlab="Mean gene expression", cex.lab=1.6, cex.axis=1.6)
\dontrun{
  dev.off()
}

# small sample solution Hothorn & Lausen

maxstat.test(Surv(time, cens) ~ MGE, data=DLBCL,
             smethod="LogRank", pmethod="HL")

# normal approximation

maxstat.test(Surv(time, cens) ~ MGE, data=DLBCL,
             smethod="LogRank", pmethod="exactGauss", iscores=TRUE,
             abseps=0.01)

# conditional Monte-Carlo

maxstat.test(Surv(time, cens) ~ MGE, data=DLBCL,
             smethod="LogRank", pmethod="condMC", B = 9999) 

# survival analysis and plotting like in Alizadeh et al. (2000)

if(require(survival, quietly = TRUE)) {

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

  \dontrun{
   # postscript("survDLBCL.ps",horizontal=F, width=8, height=8)
    pdf("survDLBCL.pdf", width=8, height=8)

  }
  par(mai=c(1.0196235, 1.0196235, 0.8196973, 0.4198450))

  plot(survfit(Surv(time, cens) ~ splitMGE, data=DLBCL),
       xlab = "Survival time in month",
       ylab="Probability", cex.lab=1.6, cex.axis=1.6, lwd=2)

  text(90, 0.9, expression("Mean gene expression" > 0.186), cex=1.6)   
  text(90, 0.45, expression("Mean gene expression" <= 0.186 ), cex=1.6)   

  \dontrun{
    dev.off()
  }
}

}
\keyword{datasets}
