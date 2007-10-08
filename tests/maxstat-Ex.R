### * <HEADER>
###
attach(NULL, name = "CheckExEnv")
assign(".CheckExEnv", as.environment(2), pos = length(search())) # base
## This plot.new() patch has no effect yet for persp();
## layout() & filled.contour() are now ok
assign("plot.new",
       function() {
	   .Internal(plot.new())
	   pp <- par(c("mfg","mfcol","oma","mar"))
	   if(all(pp$mfg[1:2] == c(1, pp$mfcol[2]))) {
               outer <- (oma4 <- pp$oma[4]) > 0; mar4 <- pp$mar[4]
               mtext(paste("help(", ..nameEx, ")"), side = 4,
                     line = if(outer)max(1, oma4 - 1) else min(1, mar4 - 1),
                     outer = outer, adj = 1, cex = .8, col = "orchid")
	   }
       },
       env = .CheckExEnv)
assign("cleanEx",
       function(env = .GlobalEnv) {
	   rm(list = ls(envir = env, all.names = TRUE), envir = env)
           assign(".Random.seed", as.integer(c(0, rep(7654, 3))), envir=.GlobalEnv)
       },
       env = .CheckExEnv)
assign("..nameEx", "__{must remake R-ex/*.R}__", env = .CheckExEnv) # for now
assign("ptime", proc.time(), env = .CheckExEnv)
postscript("maxstat-Examples.ps")
assign("par.postscript", par(no.readonly = TRUE), env = .CheckExEnv)
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"))
library('maxstat')

cleanEx(); ..nameEx <- "DLBCL"

### * DLBCL

### Name: DLBCL
### Title: Diffuse large B-cell lymphoma
### Aliases: DLBCL
### Keywords: datasets

### ** Examples


data(DLBCL)

# compute the cutpoint and plot the empirical process 

mod <- maxstat.test(Surv(time, cens) ~ MGE, data=DLBCL, smethod="LogRank")

print(mod)

##Don't run: 
##D   # postscript("statDLBCL.ps", horizontal=F, width=8, height=8)
##D   pdf("statDLBCL.pdf", width=8, height=8)

par(mai=c(1.0196235, 1.0196235, 0.8196973, 0.4198450))
plot(mod, cex.lab=1.6, cex.axis=1.6, xlab="Mean gene expression",lwd=2)
##Don't run: 
##D   dev.off()


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

##Don't run: 
##D #  postscript(file="RNewsStat.ps",horizontal=F, width=8, height=8)
##D    pdf("RNewsStat.pdf", width=8, height=8)
##D 

par(mai=c(1.0196235, 1.0196235, 0.8196973, 0.4198450))
plot(mod, xlab="Mean gene expression", cex.lab=1.6, cex.axis=1.6)
##Don't run: 
##D   dev.off()


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

  ##Don't run: 
##D    # postscript("survDLBCL.ps",horizontal=F, width=8, height=8)
##D     pdf("survDLBCL.pdf", width=8, height=8)
##D 
##D   
  par(mai=c(1.0196235, 1.0196235, 0.8196973, 0.4198450))

  plot(survfit(Surv(time, cens) ~ splitMGE, data=DLBCL),
       xlab = "Survival time in month",
       ylab="Probability", cex.lab=1.6, cex.axis=1.6, lwd=2)

  text(90, 0.9, expression("Mean gene expression" > 0.186), cex=1.6)   
  text(90, 0.45, expression("Mean gene expression" <= 0.186 ), cex=1.6)   

  ##Don't run: 
##D     dev.off()
##D   
}




par(get("par.postscript", env = .CheckExEnv))
cleanEx(); ..nameEx <- "corrmsrs"

### * corrmsrs

### Name: corrmsrs
### Title: Correlation Matrix
### Aliases: corrmsrs
### Keywords: misc

### ** Examples


# matrix of hypothetical prognostic factors

X <- matrix(rnorm(30), ncol=3) 

# this function

a <- corrmsrs(X, minprop=0, maxprop=0.9999)

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

tc <- testcorr(X)
tc <- tc + t(tc)
diag(tc) <- 1
stopifnot(all.equal(tc, a))




cleanEx(); ..nameEx <- "hohnloser"

### * hohnloser

### Name: hohnloser
### Title: Left ventricular ejection fraction of patients with malignant
###   ventricular tachyarrhythmias.
### Aliases: hohnloser
### Keywords: datasets

### ** Examples


data(hohnloser)

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




cleanEx(); ..nameEx <- "maxstat.test"

### * maxstat.test

### Name: maxstat.test
### Title: Maximally Selected Rank and Statistics
### Aliases: maxstat.test maxstat.test.data.frame maxstat.test.default
###   maxstat
### Keywords: htest

### ** Examples


x <- sort(runif(20))
y <- c(rnorm(10), rnorm(10, 2))
mydata <- data.frame(cbind(x,y))

mod <- maxstat.test(y ~ x, data=mydata, smethod="Wilcoxon", pmethod="HL",
                    minprop=0.25, maxprop=0.75, alpha=0.05)
print(mod)
plot(mod)

# adjusted for more than one prognostic factor.

data(DLBCL)

mstat <- maxstat.test(Surv(time, cens) ~ IPI + MGE, data=DLBCL, 
                      smethod="LogRank", pmethod="exactGauss", 
                      abseps=0.01)
plot(mstat)




cleanEx(); ..nameEx <- "pLausen92"

### * pLausen92

### Name: pLausen92
### Title: Approximating Maximally Selected Statistics
### Aliases: pLausen92 qLausen92
### Keywords: distribution

### ** Examples


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




cleanEx(); ..nameEx <- "pLausen94"

### * pLausen94

### Name: pLausen94
### Title: Approximating Maximally Selected Statistics
### Aliases: pLausen94 qLausen94
### Keywords: distribution

### ** Examples


p <- pLausen94(2.5, 20, 0.25, 0.75)

# Lausen 94, page 489

if (round(p, 3) != 0.073) stop("error checking pLausen94")

# the same

p2 <- pLausen94(2.5, 200, 0.25, 0.75, m=seq(from=50, to=150, by=10))

stopifnot(all.equal(round(p,3), round(p2,3)))




cleanEx(); ..nameEx <- "pexactgauss"

### * pexactgauss

### Name: pexactgauss
### Title: Computing Maximally Selected Gauss Statistics
### Aliases: pexactgauss qexactgauss
### Keywords: distribution

### ** Examples


x <- rnorm(20)

pexact <- pexactgauss(2.5, x, abseps=0.01)




cleanEx(); ..nameEx <- "plot.maxtest"

### * plot.maxtest

### Name: plot.maxtest
### Title: Print and Plot Standardized Statistics
### Aliases: plot.maxtest print.maxtest plot.mmaxtest print.mmaxtest
### Keywords: htest

### ** Examples


x <- sort(runif(20))
y <- rbinom(20, 1, 0.5)
mydata <- data.frame(c(x,y))

mod <- maxstat.test(y ~ x, data=mydata, smethod="Median", 
                    pmethod="HL", alpha=0.05)
print(mod)
plot(mod)




cleanEx(); ..nameEx <- "pmaxstat"

### * pmaxstat

### Name: pmaxstat
### Title: Approximating Maximally Selected Statistics
### Aliases: pmaxstat qmaxstat
### Keywords: distribution

### ** Examples


pmaxstat(2.5, 1:20, 5:15)




cleanEx(); ..nameEx <- "sphase"

### * sphase

### Name: sphase
### Title: S-phase fraction of tumor cells
### Aliases: sphase
### Keywords: datasets

### ** Examples

data(sphase)
maxstat.test(Surv(RFS, cens) ~ SPF, data=sphase, smethod="LogRank",
pmethod="Lau94")
maxstat.test(Surv(RFS, cens) ~ SPF, data=sphase, smethod="LogRank",
pmethod="Lau94", iscores=TRUE)
maxstat.test(Surv(RFS, cens) ~ SPF, data=sphase, smethod="LogRank",
pmethod="HL")
maxstat.test(Surv(RFS, cens) ~ SPF, data=sphase, smethod="LogRank",
pmethod="condMC", B = 9999)
plot(maxstat.test(Surv(RFS, cens) ~ SPF, data=sphase, smethod="LogRank"))

cleanEx(); ..nameEx <- "treepipit"

### * treepipit
 
### Name: treepipit
### Title: Tree Pipit Data
### Aliases: treepipit
### Keywords: datasets
 
### ** Examples

data(treepipit)
mod <- maxstat.test(counts ~ coverstorey, data = treepipit, 
                    smethod = "Data", pmethod = "HL", minprop = 0.2,
                    maxprop = 0.8)
print(mod)
plot(mod)

### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", env = .CheckExEnv),"\n")
dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
