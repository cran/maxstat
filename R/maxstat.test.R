# $Id: maxstat.test.R,v 1.10 2002/09/17 09:10:27 hothorn Exp $

maxstat.test <- function(formula, data, ...) 
  UseMethod("maxstat.test", data)

maxstat.test.default <- function(formula, data, ...)
 stop(paste("Do not know how to handle objects of class", class(data)))

maxstat.test.data.frame <-
function(formula, data, subset, na.action, ...)
{
    cl <- match.call()
    if(missing(formula)
       || (length(formula) != 3)
       || (length(attr(terms(formula[-3]), "term.labels")) != 1))
        stop("formula missing or incorrect")
    if(missing(na.action))
        na.action <- getOption("na.action")
    m <- match.call(expand.dots = FALSE)
    mt <- terms(formula, data=data)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")
    progfact <- attr(attr(mf, "terms"), "term.labels")
    MULTIMAX <- (length(progfact) > 1)
    RNAME <- names(mf)[response]
    PNAMES <- names(mf)[-response]
    X <- mf[,progfact] 
    y <- mf[[response]]
    arg <- list(...)
    DATA <- list(y=y, x=X)
    names(DATA) <- c("y", "x")
    mod <- do.call("maxstat", c(DATA, arg))
    if (MULTIMAX) {
      mod$maxstats <- sapply(mod$maxstats, function(x) { 
                             x$data.name <- paste(c(RNAME, x$data.name), 
                               collapse=" by "); list(x); } )
      mod$data.name <- paste(RNAME, paste(PNAMES, collapse=" + "), collapse
                             = " by ")
    } else {
      mod$data.name <- paste(c(RNAME, PNAMES), collapse=" by ")
    }
    mod$call <- cl
    return(mod)
}


maxstat <- function(y, x=NULL, smethod=c("Wilcoxon", 
         "Median", "NormalQuantil","LogRank"), 
         pmethod=c("none", "Lau92", "Lau94", "exactGauss", "HL", "min"), 
         iscores=(pmethod == "HL"), minprop=0.1, maxprop=0.9,  
         alpha=NULL, ...) 
{

  if (is.null(x)) stop("no data given")
  MULTIMAX <- is.matrix(x) || is.data.frame(x)
  smethod <- match.arg(smethod)
  pmethod <- match.arg(pmethod)
  
  scores <- cscores(y, type=smethod, int=FALSE)
  if (iscores & sum(scores - floor(scores)) != 0) {
    # check for midranks (Wilcoxon z.B.)
    fscore <- scores - floor(scores)
    if (all(fscore[fscore != 0] == 0.5))
      scores <- 2*scores
    else {
      # and handle real scores the way Hothorn & Lausen 2002 suggest
      scores <- scores - min(scores)
      scores <- round(scores*length(scores)/max(scores))
    }
  }

  if (MULTIMAX) {
    if (pmethod=="none") 
      stop("pmethod not specified.")
    if (!is.null(alpha) & MULTIMAX) 
      warning("cannot compute quantiles for more than one variable")
    mmax <- vector(mode="list", length=ncol(x))
    pvalues <- rep(0, ncol(x))
    statistics <- rep(0, ncol(x))
    cx <- x
    factors <- sapply(x, is.factor)
    for (i in 1:ncol(x)) {
      mmax[[i]] <- cmaxstat(scores, x[,i], pmethod, minprop, 
                           maxprop, alpha, ...)
      mmax[[i]]$data.name <- colnames(x)[i]
      mmax[[i]]$smethod <- smethod
      mmax[[i]]$pmethod <- pmethod 
      if (factors[i]) {
        cx[,i] <- as.integer(cx[,i])
        cr <- range(as.integer(mmax[[i]]$cuts))
      } else { 
        cr <- range(mmax[[i]]$cuts)
      }
      cx[cx[,i] < cr[1],i] <- cr[1]
      cx[cx[,i] > cr[2],i] <- cr[2]
      pvalues[i] <- mmax[[i]]$p.value
      statistics[i] <- mmax[[i]]$statistic
    }
    STATISTIC <- max(statistics)
    cm <- cmatrix(cx)
    cm <- cm + t(cm)
    diag(cm) <- diag(cm)/2
    p <- pmvnorm(lower=-STATISTIC, upper=STATISTIC, mean=rep(0,ncol(cm)),
                 corr=cm, ...)
    if (attr(p, "msg") != "Normal Completion") {
      msg <- paste("pvmnorm: ", attr(p, "msg"), collapse=" ")
      warning(msg)
    }
    RET <- list(maxstats=mmax, whichmin=which.min(pvalues), p.value=1 - p)
    class(RET) <- "mmaxtest"
  } else {
    RET <- cmaxstat(scores, x, pmethod, minprop, maxprop, alpha, ...)
  }
  RET$smethod <- smethod
  RET$pmethod <- pmethod 
  RET
}


cmaxstat <- function(y, x=NULL, pmethod=c("none", "Lau92", "Lau94",
          "exactGauss", "HL", "min"), minprop = 0.1, 
          maxprop=0.9, alpha = NULL, ...)
{
  pmethod <- match.arg(pmethod)

  if (is.null(y) || is.null(x)) stop("no data given")

  xname <- deparse(substitute(x))
  yname <- deparse(substitute(y))  
  if (length(xname) == 1 & length(yname) == 1)
    DNAME <- paste(xname, "and", yname)
  else
    DNAME <- "y by x" 

  N <- length(y)

  y <- y[order(x)]
  x <- sort(x)
  ties <- duplicated(x)

  m <- which(!ties) - 1 
  if (minprop == 0 & maxprop==1) m <- m[2:(length(m)-1)] else {
    if (all(m < floor(N*minprop))) stop("minprop too large")
    if (all(m > floor(N*maxprop))) stop("maxprop too small")
    m <- m[m >= floor(N*minprop)]
    m <- m[m <= floor(N*maxprop)]
  }

  if(length(m) < 1) stop("no data between minprop, maxprop")

  ss <- sum(y)
  E <- m/N*ss
  V <- m*(N-m)/(N^2*(N-1))*(N*sum(y^2) - ss^2)

  Test <- abs((cumsum(y)[m] - E)/sqrt(V))

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
    PVAL <- pexactgauss(STATISTIC, N, m, ...)
    if (!is.null(alpha))
       QUANT <- qexactgauss(alpha, N, m, ...)
  }
  if (pmethod == "HL") {
    PVAL <- pmaxstat(STATISTIC, y, m)
    if (!is.null(alpha))
       QUANT <- qmaxstat(alpha, y, m)
  }
  if (pmethod == "min") {
    PVAL <- min(pLausen92(STATISTIC, minprop, maxprop), 
                pLausen94(STATISTIC, N, minprop, maxprop, m=m),  
                pexactgauss(STATISTIC, N, m, ...),
                pmaxstat(STATISTIC, y, m))
    if (!is.null(alpha))
       QUANT <- NA
  } 

  RVAL <- list(statistic = STATISTIC, p.value = PVAL,
               method = pmethod,
               estimate = ESTIMATOR, data.name = DNAME,
               stats = Test, cuts = x[m], quant = QUANT)
  class(RVAL) <- "maxtest"
  RVAL
}
