# $Id: maxstat.R,v 1.31 2001/12/14 13:11:54 hothorn Exp $

maxstat.test <- function(x, ...) UseMethod("maxstat.test")

maxstat.test.default <-
function(x, y=NULL, event = NULL, smethod=c("Gauss", "Wilcoxon", "Median",
         "NormalQuantil","LogRank"), pmethod=c("none", "Lau92", "Lau94",
          "exactGauss", "HL", "min"), minprop = 0.1, maxprop=0.9, alpha =
          NULL, ...)
{
  smethod <- match.arg(smethod)
  pmethod <- match.arg(pmethod)

  if (is.null(y) || is.null(x)) stop("no data given")
  if (length(y) != length(x)) stop("unequal sample sizes")

  xname <- deparse(substitute(x))
  yname <- deparse(substitute(y))  
  if (length(xname) == 1 & length(yname) == 1)
    DNAME <- paste(xname, "and", yname)
  else
    DNAME <- "y by x" 

  N <- length(y)

  if (!is.null(event)) {
    if (smethod != "LogRank") warning("eventoring only meaningful with LogRank")
    if (!is.logical(event)) {
      one <- event[event == 1]
      zero <- event[event == 0]
      if (length(c(one, zero)) != N)
        stop("only 0/TRUE (alive) and 1/FALSE (dead) allowed for eventoring")
    } else {
      event <- as.numeric(event)
    }
  }

  y <- y[order(x)]
  if (!is.null(event)) event <- event[order(x)]
          else event <- rep(1, length(y))
  x <- sort(x)
  ties <- duplicated(x)

  m <- which(!ties)
  if (minprop == 0 & maxprop==1) m <- m[2:(length(m)-1)] else {
    if (all(m < floor(N*minprop))) stop("minprop too large")
    if (all(m > floor(N*maxprop))) stop("maxprop too small")
    m <- m[m >= floor(N*minprop)]
    m <- m[m <= floor(N*maxprop)]
  }

  if(length(m) < 1) stop("no data between minprop, maxprop")

  if (smethod == "Gauss") {
    cu <- cumsum(y)
    G <- cu[m]/m - (sum(y) - cu[m])/(N-m)

    Test <- abs(sqrt(m*(N-m)/N)*G)
    STATISTIC <- max(G)
          ESTIMATOR <- min(which(G == STATISTIC))
    names(ESTIMATOR) <- c("estimated cutpoint")

    if (pmethod == "none")
      PVAL <- NA
    if (pmethod == "Lau92")
      PVAL <- pLausen92(STATISTIC, minprop, maxprop)
    if (pmethod == "Lau94")
      PVAL <- pLausen94(STATISTIC, N, minprop, maxprop, m=m)
    if (pmethod == "exactGauss")
      PVAL <- pexactgauss(STATISTIC, N, m)
    if (pmethod == "HL") {
      PVAL <- NA
      warning("cannot compute HL p-value for Gauss statistic")
    }
    if (pmethod == "min") {
      warning("minimum p-value meaningless, reporting the exact one.")
      PVAL <- pexactgauss(STATISTIC, N, m)
    }
  } else {
    if (smethod == "Wilcoxon") {
      scores <- rank(y)
      if (!all(round(scores) == scores))
        scores <- 2*scores
      scores <- scores - min(scores)
    }
    if (smethod == "Median") {
      scores <- rank(y)
      scores[scores <= (N+1)/2] <- 0
      scores[scores > 0] <- 1
    }
    if (smethod == "NormalQuantil") {
      scores <- qnorm(rank(y)/(N+1))
      scores <- scores - min(scores)
      scores <- round(scores*N/max(scores))
    }
    if (smethod == "LogRank") {
      sc <- rep(0, N)
      intr <- intrank(y)
      for (i in 1:N) {
        indx <- which(y <= y[i])
        sc[i] <- event[i] - sum(event[indx]/(N - intr[indx] + 1))
      }
      scores <- sc
      scores <- scores - min(scores)
      scores <- round(scores*N/max(scores))
    }
    
    E <- m/N*sum(scores)
    V <- m*(N-m)/(N^2*(N-1))*(N*sum(scores^2) - sum(scores)^2)

    Test <- abs((cumsum(scores)[m] - E)/sqrt(V))


    STATISTIC <- max(Test)
    ESTIMATOR <- x[m[min(which(Test == STATISTIC))]]
    names(STATISTIC) <- "M"
    names(ESTIMATOR) <- c("estimated cutpoint")

    if (is.null(alpha)) QUANT <- NA

    if (pmethod == "none") {
      PVAL <- NA
      QUANT <- NA
    }
    if (pmethod == "Lau92") {
      PVAL <- pLausen92(STATISTIC, minprop, maxprop)
      if (!is.null(alpha))
        QUANT <- qLausen92(alpha, minprop, maxprop)
    }
    if (pmethod == "Lau94") {
      PVAL <- pLausen94(STATISTIC, N, minprop, maxprop, m=m)
      if (!is.null(alpha))
         QUANT <- qLausen94(alpha, N, minprop, maxprop, m=m)
    }
    if (pmethod == "exactGauss") {
      PVAL <- pexactgauss(STATISTIC, N, m)
      if (!is.null(alpha))
         QUANT <- qexactgauss(alpha, N, m)
    }
    if (pmethod == "HL") {
      PVAL <- pmaxstat(STATISTIC, scores, m)
      if (!is.null(alpha))
         QUANT <- qmaxstat(alpha, scores, m)
    }
    if (pmethod == "min") {
      PVAL <- min(pLausen92(STATISTIC, minprop, maxprop), 
                  pLausen94(STATISTIC, N, minprop, maxprop, m=m),  
                  pexactgauss(STATISTIC, N, m),
                  pmaxstat(STATISTIC, scores, m))
      if (!is.null(alpha))
         QUANT <- NA
    } 
  }

  RVAL <- list(statistic = STATISTIC, p.value = PVAL,
               method = paste(smethod, "using", pmethod),
               estimate = ESTIMATOR, data.name = DNAME,
               stats = Test, cuts = x[m], quant = QUANT)
  class(RVAL) <- "maxtest"
  RVAL
}

maxstat.test.formula <-
function(formula, data, subset, na.action, ...)
{
    if(missing(formula)
       || (length(formula) != 3)
       || (length(attr(terms(formula[-2]), "term.labels")) != 1)
       || (length(attr(terms(formula[-3]), "term.labels")) != 1))
        stop("formula missing or incorrect")
    if (length(formula[[2]]) > 1) {
        if (formula[[2]][[1]] == as.name("Surv")) 
            formula[[2]][[1]] <- as.name("msurv")
    }
    if(missing(na.action))
        na.action <- getOption("na.action")
    m <- match.call(expand.dots = FALSE)
    if (length(m$formula[[2]]) > 1) {
        if (m$formula[[2]][[1]] == as.name("Surv"))
            m$formula[[2]][[1]] <- as.name("msurv")
    }
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    DNAME <- gsub("msurv", "Surv", DNAME)
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    nc <- ncol(mf[[response]])
    if (!is.null(nc)) {
      if (nc == 2) {
        if (!class(mf[[response]]) == "msurv") stop("response is not of class msurv")
        event <- mf[[response]][,2]
        y <- mf[[response]][,1]
        arg <- list(...)
        if (!is.null(arg$event)) warning("survival formula and argument event given. Using survival formula.")
        arg <- c(arg, event = list(event))
      }
    } else {
      y <- mf[[response]]
      arg <- list(...)
    }
    DATA <- list(x = mf[[-response]], y = y)
    names(DATA) <- c("x", "y")
    y <- do.call("maxstat.test", c(DATA, arg))
    y$data.name <- DNAME
    return(y)
}



print.maxtest <- function(x, digits = 4, quote = TRUE, prefix = "", ...) {
  x$stats <- NULL
  x$cuts <- NULL
  x$quant <- NULL
  class(x) <- "htest"
  print(x, digits = digits, quote = quote, prefix = prefix, ...)
} 

plot.maxtest <- function(x, xlab=NULL, ylab=NULL, ...) {
  xname <- unlist(strsplit(x$data.name, "by"))[2]
  if (is.na(x$quant)) {
    smethod <- gsub(" ", "", unlist(strsplit(x$method, "using"))[1])
    if (smethod == "LogRank") smethod <- "log-rank"
    if (is.null(ylab)) ylab <- paste("Standardized", smethod, "statistic")
    if (is.null(xlab)) xlab <- xname
    plot(x$cuts, x$stats, type="b", xlab=xlab, ylab=ylab, ...)
    lines(c(x$estimate, x$estimate), c(0, x$statistic), lty=2)
  } else {
    smethod <- gsub("LogRank", "log-rank", x$method)
    smethod <- gsub(" ", "",unlist(strsplit(smethod, "using")))
    ylim <- c(min(x$quant, min(x$stats)), max(x$quant, max(x$stats)))
    ylim <- c(ylim[1]*0.95, ylim[2]*1.05)
    xlength <- range(x$cuts)
    if (is.null(ylab)) ylab <- paste("Standardized", smethod[1],
                                     "statistic using", smethod[2])
    if (is.null(xlab)) xlab <- xname
    plot(x$cuts, x$stats, type="b", xlab=xlab, ylab=ylab, ylim=ylim, ...)
    lines(c(x$estimate, x$estimate), c(0, x$statistic), lty=2)
    lines(xlength, c(x$quant, x$quant), col="red")
  }
}

msurv <- function(time, event) {
  RET <- cbind(time, event)
  class(RET) <- "msurv"
  RET
}

pLausen92 <- function(b, minprop=0.1, maxprop=0.9)
{
  db <- dnorm(b)
  p <- 4*db/b + db*(b - 1/b)*log((maxprop*(1 - minprop))/((1-maxprop)*minprop))
  max(p,0)
}

qLausen92 <- function(p, minprop=0.1, maxprop=0.9)
{
  test <- function(x)
    abs(pLausen92(x, minprop, maxprop) - p)

  return(optim(2, test)$par)
}

pLausen94 <- function(b, N, minprop=0.1, maxprop=0.9, m=NULL)
{
  if(is.null(m))
    m <- floor(N*minprop):floor(N*maxprop)
  m1 <- m[1:(length(m)-1)]
  m2 <- m[2:length(m)]
  t <- sqrt(1 - m1*(N-m2)/((N-m1)*m2))
  D <- sum(1/pi*exp(-b^2/2)*(t - (b^2/4 -1)*(t^3)/6))
  1 - (pnorm(b) - pnorm(-b)) + D
}

qLausen94 <- function(p, N, minprop=0.1, maxprop=0.9, m=NULL)
{
  test <- function(x)
    abs(pLausen94(x, N, minprop, maxprop, m) - p)

  return(optim(2, test)$par)
}
  
p2normG <- function(h,k,rho, maxp=3000)
{
  # multivariate normal according to Genz/Bretz
  sigma <- diag(2)
  sigma[1,2] <- rho
  sigma[2,1] <- rho
  prob <- pmvnorm(lower=c(-Inf, -Inf), upper=c(h, k), mean=c(0,0), corr=sigma)
  prob
}

p2norm <- function(h,k,rho, maxp=3000)
{
  # multivariate normal according to Drezner (cited in Schlittgen)
  fct <- function(r) 
    1/sqrt(1 - r^2)*exp(-(h^2 -2*sqrt(1-r^2)*h*k + k^2)/(2*r^2))
     
  if (rho < 1) {
    integral <- integrate(fct, 0, sqrt(1-rho^2), maxpts=maxp)
    return((pnorm(-max(h, k)) - 1/(2*pi)*integral$value))
  } else return(pnorm(-max(h, k)))
}

pbvmax <- function(cor, t, method=c("drezner", "genz"))
{
  if (method == "drezner") {
    prob <- p2norm(t, t, cor) - p2norm(-t, t, cor)
    prob <- prob - p2norm(t, -t, cor) + p2norm(-t,-t,cor)
  } else {
    prob <- p2normG(t, t, cor) - p2normG(-t, t, cor)
    prob <- prob - p2normG(t, -t, cor) + p2normG(-t,-t,cor)
  }
  prob
}

pSchlitt <- function(b, N, m, method=c("genz", "drezner"))
{
  n <- m[2:length(m)]
  m <- m[1:(length(m)-1)]
  mcorr <- sqrt((m/N)*(1 - n/N))/sqrt(n/N*(1 - m/N))

  prob <- sapply(mcorr, pbvmax, t=b, method=method)
  # note: length(m) = k - 1
  prob <- sum(log(prob)) - (length(m) + 1 - 2)*log(pnorm(b) - pnorm(-b))
  1-exp(prob)
}

qSchlitt <- function(p, N, m)
{
  test <- function(x)
    abs(pSchlitt(x, N, m) - p)

  return(optim(2, test)$par)
}

pexactgauss <- function(b, N, m, maxpts=25000)
{
  if (!require(mvtnorm)) stop("package mvtnorm not loaded")
  n <- m[2:length(m)]
  mm <- m[1:(length(m)-1)]
  mcorr <- sqrt((mm/N)*(1 - n/N))/sqrt(n/N*(1 - mm/N))

  corrmatrix <- diag(length(m))

  for (i in 1:(length(m)-1))
    corrmatrix[i,(i+1):length(m)] <- cumprod(mcorr[i:(length(m)-1)])

  p <- pmvnorm(mean=rep(0, length(m)),
               corr=t(corrmatrix), lower=rep(-b, length(m)),
               upper=rep(b, length(m)), maxpts=maxpts)
  attributes(p) <- NULL
  1 - p
}

qexactgauss <- function(p, N, m)
{
  test <- function(x)
    abs(pexactgauss(x, N, m) - p)

  return(optim(2, test)$par)
}


pmaxstat <- function(b, scores,  msample, quant=FALSE)
{
  if (!require(exactRankTests)) stop("package exactRankTests not loaded")

  # for intergers only

  if (!all(round(scores) == scores))
    stop("scores are not integers in pmaxstat")
  if (length(scores) < length(msample))
    stop("incorrect number of cutpoints in pmaxstat")
  if (length(b) != 1)
    stop("b is not a single number in pmaxstat")

  N <- length(scores)

  scores <- equiscores(scores, N)$scores

  if(sum(scores) > sum(1:200)) { 
    warning("Cannot compute SR p-value. Sum of scores > 20100")
    p <- list(1, 1)
    names(p) <- c("upper", "lower")
    return(p)
  }

  H <- rep(0, sum(scores)*N)

  totsum <- sum(scores)
  sc <- rep(1, N)

  # Streitberg / Roehmel in C, package "exactRankTest"

  if (!is.loaded("cpermdist2")) stop("Function cpermdist2 from package exactRankTests not found!")

  H <- .C("cpermdist2", H = as.double(H), as.integer(N),
                as.integer(totsum), as.integer(sc),
    as.integer(scores), as.integer(N),
                as.integer(length(H)))$H

  # add last row, column for compatibility

  H <- matrix(H, nrow=N, ncol=totsum, byrow=TRUE)
  H <- rbind(H, rep(0, ncol(H)))
  H <- cbind(H, c(rep(0, nrow(H)-1), 1))

  # sample sizes

  m <- 1:(N-1)

  # Expectation and Variance of a Linear Rank Test

  E <- m/N*sum(scores)
  V <- m*(N-m)/(N^2*(N-1))*(N*sum(scores^2) - sum(scores)^2)
  S <- rep(1:(ncol(H)-1), nrow(H) -2)
  S <- matrix(S, nrow(H) -2, ncol(H)-1, byrow=TRUE)
  EM <- matrix(rep(E, ncol(H) -1), nrow(H) -2, ncol(H) - 1)
  VM <- matrix(rep(V, ncol(H) -1), nrow(H) -2, ncol(H) - 1 )
  S <- (S- E)/sqrt(V)

  # remove technical parts

  H <- H[2:(nrow(H)-1), ]
  H <- H[, 2:(ncol(H))]

  # S is the matrix of the standardized values

  S <- abs(S)
  S[H == 0] <- 0

  # extend to number of permutations

  H <- H*gamma(m+1)*gamma(N -m +1)

  if (quant)
    return(list(scores=scores, H=H, E=E, S=S, msample=msample, N=N))

  # those are in general not needed

  H[S <= b] <- 0


  # delete those, which are in m+1 and + max(scores) still > b
  # well, that's the trick

  sH <- apply(H, 1, sum)

  for (i in min(msample):(nrow(H)-1)) {
    indx <- which(H[i,] > 0)
    if (length(indx) > 0) {
      indxmax <- indx[indx < E[i]]
      indxmax <- indxmax[S[i+1, indxmax + max(scores)] > b]
      if (length(indxmax) > 0 & all(!is.na(indxmax)))
        sH[i+1] <- sH[i+1] - sum(H[i, indxmax]) 
      indxmin <- indx[indx > E[i]]
      indxmin <- indxmin[S[i+1, indxmin + min(scores)] > b]
      if (length(indxmin) > 0 & all(!is.na(indxmin)))
        sH[i+1] <- sH[i+1] - sum(H[i, indxmin])
    }
  }

  # only meaningful sample sizes

  sH <- sH[msample]

  gaN <- gamma(N+1)   
  lower <- min(sum(sH)/gaN, 1)  # <- this is a better approx.
  #  upper <- max(apply(H, 1, sum))/gaN
  #  p <- list(upper, lower)
  #  names(p) <- c("upper", "lower")
  #  cat("hl working: ", all.equal(hl(scores, H, E, S, msample, N, b), p), "\n")
  lower
}

qmaxstat <- function(p, scores, msample)
{
  if (p > 1 | p < 0) stop("p not in [0,1]")
  sr <- pmaxstat(3, scores, msample, quant=TRUE)
  tr <- rev(sort(unique(round(sr$S,5))))
  i <- 1
  pw <- 0
  while (pw < p) {
    pw <- hl(sr$scores, sr$H, sr$E, sr$S, sr$msample, sr$N, tr[i])$lower
    i <- i+1
  }
  return(tr[i-1])
}

intrank <- function(x)
{
  r <- rank(x)
  for (i in sort(unique(r[duplicated(r)]))) {
    n <- length(r[r == i])
    r[r == i] <- i + sum(1:(n-1))/n
  }
  r
}

hl <- function(scores, H, E, S, msample, N, b)
{
  H[S <= b] <- 0

  # delete those, which are in m+1 and + max(scores) still > b
  # well, that's the trick

  sH <- apply(H, 1, sum)

  for (i in min(msample):(nrow(H)-1)) {
    indx <- which(H[i,] > 0)
    if (length(indx) > 0) {
      indxmax <- indx[indx < E[i]]
      indxmax <- indxmax[S[i+1, indxmax + max(scores)] > b]
      if (length(indxmax) > 0 & all(!is.na(indxmax)))
        sH[i+1] <- sH[i+1] - sum(H[i, indxmax]) 
      indxmin <- indx[indx > E[i]]
      indxmin <- indxmin[S[i+1, indxmin + min(scores)] > b]
      if (length(indxmin) > 0 & all(!is.na(indxmin)))
        sH[i+1] <- sH[i+1] - sum(H[i, indxmin])
    }
  }

  # only meaningful sample sizes

  sH <- sH[msample]

  gaN <- gamma(N+1)   
  lower <- min(sum(sH)/gaN, 1)  # <- this is a better approx.
  upper <- max(apply(H, 1, sum))/gaN
  p <- list(upper, lower)
  names(p) <- c("upper", "lower")
  p
}
