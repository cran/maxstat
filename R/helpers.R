


corrmatrix <- function(X) {
  N <- nrow(X)
  if (is.null(N)) N <- length(X)
  k <- ncol(X)
  if (is.null(k)) k <- 1
  cut <- c()
  for (j in 1:k) {
    if (k > 1)
      points <- X[!duplicated(X[,j]),j]
    else
      points <- X[!duplicated(X)]
    points <- points[points != max(points)]
    cut <- rbind(cut, cbind(points, j))
  }
  corrmatrix <- matrix(rep(0, nrow(cut)^2), ncol=nrow(cut))
  for (i in 1:nrow(cut)) {
    for (j in i:nrow(cut)) {
      if (k > 1) {
        indxone <- which(X[,cut[i,2]] <= cut[i,1])
        indxtwo <- which(X[,cut[j,2]] <= cut[j,1])
      } else {
        indxone <- which(X <= cut[i,1])
        indxtwo <- which(X <= cut[j,1])
      }
      corrmatrix[i,j] <- corrgauss(indxone, indxtwo, N)
    }
  }
  t(corrmatrix)
}



wh <- function(cut, x)
  which(x <= cut)

index <- function(x) {
  ux <- unique(x)
  ux <- ux[ux < max(ux)]
  lapply(ux, wh, x = x)
}


corlm <- function(Xone, Xtwo)
{
  i <- ncol(Xone)
  if (i > 1) {
    if (i != ncol(Xtwo) || !all.equal(Xone[,1:(i-1)], Xtwo[,1:(i-1)]))
       stop("error")
    sXone <- solve(t(Xone)%*%Xone)
    sXtwo <- solve(t(Xtwo)%*%Xtwo)
    num <- (sXone%*%t(Xone))[i,]%*%(Xtwo%*%sXtwo)[,i]
    den <- sqrt(sXone[i,i] * sXtwo[i,i])
    return(num/den)
  } else {
    sXone <- solve(t(Xone)%*%Xone)
    sXtwo <- solve(t(Xtwo)%*%Xtwo)
    num <- (sXone*(t(Xone)%*%Xtwo)*sXtwo)
    den <- sqrt(sXone * sXtwo)
    return(num/den)
  }
}

forsel <- function(y, Xin, Xout, ret=FALSE)
{
  k <- ncol(Xout)
  N <- nrow(Xin)
  if (is.null(N)) N <- length(Xin)
  if (!is.null(Xin)) {
    if (length(y) != N || length(y) != nrow(Xout)) stop("error")
  }
  goodvar <- c()
  sigma <- diag(k)
  gauss <- c()
  for (i in 1:k) {
    Xnewi <- cbind(Xin, Xout[,i])
    sXnewi <- solve(t(Xnewi)%*%Xnewi)
    l <- ncol(Xnewi)
    stat <- ((sXnewi%*%t(Xnewi)%*%y)[l])/sqrt(sXnewi[l,l])
    gauss <- c(gauss, stat)
    for (j in min((i+1),k):k)
      sigma[i,j] <- corlm(Xnewi, cbind(Xin, Xout[,j]))
  }
  sigma <- sigma + t(sigma) + diag(k)
  mg <- max(abs(gauss))
  indx <- which(abs(gauss) == mg)
  p <- 1-pmvnorm(rep(-mg, k), rep(mg, k), rep(0, k), t(sigma))$value
  if (p < 0.05) {
    goodvar <- c(goodvar, indx)
    newXin <- cbind(Xin, Xout[,indx])
#    cat("variable selected: ", colnames(Xout)[indx], " at p = ", p, "\n")
    goodvar <- c(goodvar, forsel(y, newXin, Xout[,!(1:k %in% indx)]))
    return(goodvar)
  } else {
    return(NULL)
  }  
}
