#' @title Simulate classes
#'
#' @description Simulate multivariate normal data.
#'
#' @param dims a 10 element vector of group sizes.
#' @param p integer number of variables.
#' @param n1 integer number of samples in each of two classes in training/calibration data.
#' @param n2 integer number of samples in each of two classes in test/validation data.
#' 
#' @details The class simulation is a straigh forward simulation of mulitvariate normal
#' data into two classes for training and test data, respectively.
#' The data simulation uses a strictly structured multivariate normal simulation for 
#' with continuous response data.
#'  
#' @return Returns a list of predictor and response data for training and testing.
#'
#' @author Tahir Mehmood, Kristian Hovde Liland, Solve Sæbø.
#'
#' @references T. Mehmood, K.H. Liland, L. Snipen, S. Sæbø, A review of variable selection 
#' methods in Partial Least Squares Regression, Chemometrics and Intelligent Laboratory Systems
#' 118 (2012) 62-69.
#' T. Mehmood, S. Sæbø, K.H. Liland, Comparison of variable selection methods in partial least
#' squares regression, Journal of Chemometrics 34 (2020) e3226.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{shaving}}, 
#' \code{\link{stpls}}, \code{\link{truncation}},
#' \code{\link{bve_pls}}, \code{\link{ga_pls}}, \code{\link{ipw_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{rep_pls}}, \code{\link{spa_pls}},
#' \code{\link{lda_from_pls}}, \code{\link{lda_from_pls_cv}}, \code{\link{setDA}}.
#' 
#' @examples
#' str(simulate_classes(5,4,4))
#' 
#' @importFrom mvtnorm rmvnorm
#' @export
simulate_classes <- function (p,n1,n2){
  sigma <- diag(1,p)
  sigma[1,2] <- 0.7
  sigma[2,1] <- 0.7
  means1 <- c(0, 2.9,matrix(0,1,p-2))
  means2 <- c(0,-2.9,matrix(0,1,p-2))
  Xtrain <- rbind(rmvnorm(n=n1, mean=means1, sigma), rmvnorm(n=n1, mean=means2, sigma))
  Ytrain <- c(matrix(1,n1,1), -matrix(1,n1,1))
  
  Xtest <- rbind(rmvnorm(n=n2, mean=means1, sigma),rmvnorm(n=n2, mean=means2, sigma)) 
  Ytest <- c(matrix(1,n2,1), -matrix(1,n2,1))
  return(list(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest))
}

#' @rdname simulate_classes
#' @export
simulate_data <- function (dims, n1=150, n2=50){
  n <- n1+n2
  Sigma <- createSmatrix(dims=dims,c(0.5,0.5,0.5,0.2,0,0,0,0,0,0),c(0.35,-0.3,0.25,0.2,0,0,0,0,0,0),1,1)
  Data  <- rmvnorm(n, mean = rep(0, nrow(Sigma)), sigma = Sigma, method=c("chol"))
  # Beta  <- solve(Sigma[-dim(Sigma)[1],-dim(Sigma)[1]])%*%Sigma[dim(Sigma)[1],-dim(Sigma)[1]]
  X <- Data[,-dim(Sigma)[1]]
  y <- Data[,dim(Sigma)[1]]
  Xtrain <- X[1:n1,];    Ytrain <- y[1:n1]
  Xtest  <- X[(n1+1):n,]; Ytest <- y[(n1+1):n]
  
  return(list(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest))
}

#' @title Matrix plotting
#'
#' @description Plot a heatmap with colorbar.
#'
#' @param x a \code{matrix} to be plotted.
#' @param main header text for the plot.
#' @param ... additional arguments (not implemented).
#'  
#' @author Tahir Mehmood, Kristian Hovde Liland, Solve Sæbø.
#'
#' @references T. Mehmood, K.H. Liland, L. Snipen, S. Sæbø, A review of variable selection 
#' methods in Partial Least Squares Regression, Chemometrics and Intelligent Laboratory Systems
#' 118 (2012) 62-69.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{shaving}}, 
#' \code{\link{stpls}}, \code{\link{truncation}},
#' \code{\link{bve_pls}}, \code{\link{ga_pls}}, \code{\link{ipw_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{rep_pls}}, \code{\link{spa_pls}},
#' \code{\link{lda_from_pls}}, \code{\link{lda_from_pls_cv}}, \code{\link{setDA}}.
#' 
#' @examples
#' myImagePlot(matrix(1:12,3,4), 'A header')
#' 
#' @import grDevices graphics
#' @export
myImagePlot <- function(x,main, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- seq_len(ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- seq_len(nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- rev(seq_len(nrow(x)))
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,7,2.5,2))
  image(seq_along(xLabels), seq_along(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max), main=main)
  if( !is.null(title) ){
    title(main=title)
  }
  axis(1, at=seq_along(xLabels), labels=xLabels, cex.axis=0.7)
  axis(2, at=seq_along(yLabels), labels=yLabels, las=1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}

##########################
# Non-exported functions #
##########################

# Environment to handle functions with reduced argument set
.PLSVarSelEnv <- new.env(parent=emptyenv())
PLSVarSelEnv <- function() .PLSVarSelEnv
putPLSVarSelEnv <- function(x, value) assign(x, value, envir=PLSVarSelEnv())
getPLSVarSelEnv <- function(x, mode="any") get(x, envir=PLSVarSelEnv(), mode=mode, inherits=FALSE)

putPLSVarSelEnv("X", 0)
putPLSVarSelEnv("y", 0)

# Class performance
ClassPerf <- function(Predicted, Orignal, index){
  Predicted[index] <- 1
  tb  <- table(factor(Predicted, levels=c(1,0)), factor(Orignal, levels=c(1,0)))
  Sen <- tb[1,1]/ sum(tb[,1])
  Spe <- tb[2,2]/ sum(tb[,2])
  Acc <- (tb[1,1]+tb[2,2])/ sum(tb)
  return(list(Sen=Sen, Spe=Spe, Acc=Acc))
}	

# RMSEP
rmsep <- function(ypred,y){
  sqrt(mean((y-ypred)^2))
}

# Prediction performance
PredPerf <- function(index, Xcal, Xtest, Ycal, Ytest, ncomp, R2 ){
  XcalN   <- Xcal[,index]
  XtestN  <- Xtest[,index]
  mydata  <- data.frame( yy=c(Ycal, Ytest) ,XX=I(rbind(XcalN, XtestN)) , train= c(rep(TRUE, length(Ycal)), rep(FALSE, length(Ytest))))
  fitt    <- plsr(yy~  XX,  ncomp=min(ncomp, (ncol(XcalN)-1)),validation = "CV", data=mydata[mydata$train,])
  usecomp <- which(fitt$valid$PRESS[1,] == min(fitt$valid$PRESS[1,]))[1]
  pred    <- predict(fitt, ncomp=usecomp, newdata=mydata[!mydata$train,])
  RMSEP   <- rmsep(pred[,1,1],Ytest)
  RMSEP_minA <- RMSEP/sqrt(1-R2)   #R2 from the input parameter of the design 
  return(list( RMSEP=RMSEP, RMSEP_minA=RMSEP_minA))
} 


# Create covariance matrix
#' @import bdsmatrix bdsmatrix
createSmatrix <- function(dims, rhos, xycors, xvar=1, yvar=1){
  # Assume there are G groups of genes
  # dims = A G-vector of gene group sizes
  # rhos = A G-vector of correlations within gene groups
  # xycors = A G-vector of correlations between gene groups and y
  # xvar = Common variance for all genes, default = 1
  # yvar = Variance of phenotype variable y, default = 1
  
  p <-sum(dims)
  
  if(length(xvar) > 1){ 
    xvar    <- rep(xvar,times=dims)        
    xvarmat <- outer(sqrt(xvar),sqrt(xvar),"*")
  } else { 
    xvarmat <- matrix(xvar,p,p)}
  covvek <- rep(xycors,times=dims)*sqrt(xvar)*sqrt(yvar)
  
  blocklist <- list()
  for(i in seq_along(dims)){
    blocklist[[i]] <- .makeSmatrix(dims[i],rhos[i])
  }
  
  covmat <- as.matrix(bdsmatrix(dims, unlist(blocklist),rmat=matrix(c(covvek,yvar),ncol=1)))
  covmat[1:p,1:p] <- covmat[1:p,1:p]*xvarmat
  
  #Positive-semi-definiteness check
  eig <- eigen(covmat)
  pd  <- all(eig$values>0)
  if(!pd) print("Warning! Matrix is not positive definite!")
  covmat
}
.makeSmatrix <- function(p,rho){
  # Function to create block of uniformly correlated variables
  # Used in function createSmatrix()
  J <- matrix(1,p,p)
  (1-rho)*diag(p)+rho*J    
}

#Example:
#Sigma <- createSmatrix(c(5,5,5,5),c(0.2,0.2,0.2,0.2),c(0.4,-0.3,0,0),1,1)
#image(Sigma)
#library(mvtnorm)
#N=3
#Data<- rmvnorm(n, mean = rep(0, nrow(Sigma)), sigma = Sigma, method=c("chol"))
#Beta<- solve(Sigma[-11,-11])%*%Sigma[11,-11]



# Prepare package internal namespace
.plsEnv <- new.env(parent=emptyenv())
plsEnv <- function() .plsEnv
putplsEnv <- function(x, value) assign(x, value, envir=plsEnv())
getplsEnv <- function(x, mode="any") get0(x, envir=plsEnv(), mode=mode, inherits=FALSE)
putplsEnv("LQ","lda")

#' Set chosen Discriminant Analysis
#' 
#' The default methods is LDA, but QDA and column of maximum prediction can be chosen.
#'
#' @param LQ character argument 'lda', 'qda', 'max' or NULL
#'
#' @return Returns the default set method.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{shaving}}, 
#' \code{\link{stpls}}, \code{\link{truncation}},
#' \code{\link{bve_pls}}, \code{\link{ga_pls}}, \code{\link{ipw_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{rep_pls}}, \code{\link{spa_pls}},
#' \code{\link{lda_from_pls}}, \code{\link{lda_from_pls_cv}}, \code{\link{setDA}}.
#' 
#' @examples
#' \dontrun{
#' setDA() # Query 'lda', 'qda' or 'max'
#' setDA('qda') # Set default method to QDA
#' }
#' @export
setDA <- function(LQ = NULL){
  if(is.null(LQ)){
    return( getplsEnv("LQ") )
  } else {
    if(!LQ %in% c("lda","qda","max")){
      stop("Argument must be 'lda', 'qda' or 'max'")
    }
    putplsEnv("LQ", LQ)
    return(LQ)
  }
}


t2_calc <- function(x, type, alpha, main){
#  type <- match.arg(type)
  phase = 1
  method = "sw"
  p <- ncol(x)
  m <- nrow(x)
  if (inherits(x, "matrix") || inherits(x, "data.frame")) 
    (x <- array(data.matrix(x), c(m, p, 1)))
  n <- dim(x)[3]
  x.jk <- matrix(0, m, p)
  t2 <- matrix(0, m, 1)
  x.jk <- apply(x, 1:2, mean)
  Xmv <- colMeans(x.jk)
  S <- covar(x, method = method)
  colm <- nrow(x)
  name <- paste("Hotelling Control Chart")
  for (ii in 1:m) {
    t2[ii, 1] <- n * t(x.jk[ii, ] - Xmv) %*% solve(S) %*% 
      (x.jk[ii, ] - Xmv)
  }
  ifelse(n == 1, ifelse(phase == 1, ucl <- ((colm - 1)^2)/colm * qbeta(1 - alpha, p/2, ((colm - p - 1)/2)), 
                        ucl <- ((p * (colm + 1) * (colm - 1))/((colm^2) - colm * p)) * qf(1 - alpha, p, colm - p)), 
         ifelse(phase == 1, ucl <- (p *  (colm - 1) * (n - 1))/(colm * n - colm - p + 1) * qf(1 - alpha, p, colm * n - colm - p + 1), 
                ucl <- (p * (colm + 1) * (n - 1))/(colm * n - colm - p + 1) * qf(1 - alpha, p, colm * n - colm - p + 1)))
  if (any(t2 > ucl)) {
    cat("The following(s) point(s) fall outside of the control limits")
    t3 <- which(t2 > ucl)
    print(t3)
    for (ii in seq_along(t3)) {
      v = 1
      k = 0
      for (i in seq_len(p)) {
        k <- k + factorial(p)/(factorial(i) * factorial(p - 
                                                          i))
      }
      q <- matrix(0, k, p + 3)
      for (i in seq_len(p)) {
        a <- t(combn(p, i))
        for (l in seq_len(nrow(a))) {
          for (j in seq_len(ncol(a))) {
            q[v, j + 3] <- a[l, j]
          }
          v = v + 1
        }
      }
      for (i in seq_len(nrow(q))) {
        b <- subset(q[i, 4:ncol(q)], q[i, 4:ncol(q)] > 
                      0)
        di <- length(b)
        if (length(b) > 1) {
          q[i, 1] <- n * t(Xmv[b] - x.jk[t3[ii], ][b]) %*% 
            solve(S[b, b]) %*% (Xmv[b] - x.jk[t3[ii], 
            ][b])
        }
        else (q[i, 1] <- n * (x.jk[t3[ii], ][b] - Xmv[b])^2/S[b, b])
        ifelse(n == 1, ifelse(phase == 1, q[i, 2] <- ((colm - 1)^2)/colm * qbeta(1 - alpha, di/2, (((2 * (colm - 1)^2)/(3 * colm - 4) - di - 1)/2)), 
                              q[i, 2] <- ((di * (colm + 1) * (colm - 1))/((colm^2) - colm * di)) * qf(1 - alpha, di, colm - di)), 
               ifelse(phase == 1, q[i, 2] <- (di * (colm - 1) * (n - 1))/(colm * n - colm - di + 1) * qf(1 - alpha, di, colm * n - colm - di + 1), 
                      q[i, 2] <- (di * (colm + 1) * (n - 1))/(colm * n - colm - di + 1) * qf(1 - alpha, di, colm * n - colm - di + 1)))
        q[i, 3] <- 1 - pf(q[i, 1], di, colm - 1)
      }
      colnames(q) <- c("t2 decomp", "ucl", "p-value", 
                       1:p)
      print(list(`Decomposition of` = t3[ii]))
      print(round(q, 4))
    }
  }
  t3 <- which(t2 > ucl)
  par(mar = c(4, 5, 3, 5))
  plot(t2, ylim = c(0, 1.1 * max(max(t2), ucl)), main = name, 
       xlab = "Sample", ylab = expression(T^2), type = "o", 
       las = 1)
  points(t3, t2[t3], col = 2)
  segments(0, ucl, m, ucl, col = 2)
  mtext(paste(" UCL=", round(ucl, 2)), side = 4, at = ucl, las = 2)
  outList = list(name, ucl = round(ucl, 2), t2 = round(t2, 2), Xmv = round(Xmv, 2), covariance = signif(S, 2))
  return(outList)
}

covar <- function (x, stat, method, ...) 
{
  p <- ncol(x)
  m <- nrow(x)
  if (inherits(x, "matrix") || inherits(x, "data.frame")) 
    (x <- array(data.matrix(x), c(m, p, 1)))
  n <- dim(x)[3]
  s.jk <- matrix(0, m, p^2)
  SS <- matrix(0, m, 1)
  if (n > 1) {
    arrays <- expand.grid(1:p, 1:p)
    for (i in 1:m) {
      for (j in 1:p^2) {
        s.jk[i, j] <- cov(x[i, arrays[j, 1], ], x[i, 
                                                  arrays[j, 2], ])
      }
    }
    S <- matrix(colMeans(s.jk), p, p)
    for (ii in 1:m) {
      SS[ii] <- det(matrix(s.jk[ii, ], p, p))
    }
    if (missing(stat)) 
      (return(S))
    else (return(SS))
  }
  if (n == 1) {
    if (missing(method)) 
      (method = "sw")
    if (method == "sw") {
      B <- matrix(0, p, p)
      w <- sweep(x, 2, (apply(x, 2, mean)))
      for (i in seq_len(m)) {
        B <- B + w[i, , ] %*% t(w[i, , ])
      }
      S <- B/(m - 1)
    }
    if (method == "hm") {
      V <- matrix(0, m - 1, p)
      for (i in seq_len(m - 1)) {
        V[i, ] <- x[i + 1, , ] - x[i, , ]
      }
      S <- 0.5 * t(V) %*% V/(m - 1)
    }
    return(S)
  }
}