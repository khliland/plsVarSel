#' Hotelling's T^2 based variable selection in PLS -- T^2-PLS)
#'
#' @description Variable selection based on the T^2 statistic. A side effect
#' of running the selection is printing of tables and production of plots.
#' 
#' @param ytr Vector of responses for model training.
#' @param Xtr Matrix of predictors for model training.
#' @param yts Vector of responses for model testing.
#' @param Xts Matrix of predictors for model testing.
#' @param ncomp Number of PLS components.
#' @param alpha Hotelling's T^2 significance levels.
#'
#' @return Parameters and variables corresponding to variable selections
#' of minimum error and minimum variable set.
#' 
#' @references Tahir Mehmood, Hotelling T^2 based variable selection in partial 
#' least squares regression, Chemometrics and Intelligent Laboratory Systems 154 (2016), pp 23-28
#'
#' @examples
#' data(gasoline, package = "pls")
#' library(pls)
#' if(interactive()){
#'   t2 <- T2_pls(gasoline$octane[1:40], gasoline$NIR[1:40,], 
#'              gasoline$octane[-(1:40)], gasoline$NIR[-(1:40),], 
#'              ncomp = 10, alpha = c(0.2, 0.15, 0.1, 0.05, 0.01))
#'   matplot(t(gasoline$NIR), type = 'l', col=1, ylab='intensity')
#'   points(t2$mv[[1]], colMeans(gasoline$NIR)[t2$mv[[1]]], col=2, pch='x')
#'   points(t2$mv[[2]], colMeans(gasoline$NIR)[t2$mv[[2]]], col=3, pch='o')
#' }
#' @importFrom stats cov qbeta
#' @importFrom utils combn
#' @export
T2_pls <- function(ytr, Xtr, yts, Xts, ncomp = 10, alpha = c(0.2, 0.15, 0.1, 0.05, 0.01)){
  pls.cv <- plsr(ytr ~ Xtr, ncomp=ncomp, validation = "CV")
  opt.comp <- which.min(pls.cv$validation$PRESS)
  opt.fit  <- plsr(ytr ~ Xtr, opt.comp)
  pls.lwd  <- data.frame(opt.fit$loading.weights[,1:opt.comp])
  R <- matrix(NA, length(alpha), 5)
  V <- vector('list', length(alpha))
  for(k in 1:length(alpha)){
    T2.chart <- t2_calc(pls.lwd, type = "t2", alpha = alpha[k], main = "T2 chart for PLS loadings")
    # Marked for deletion from CRAN
    #T2.chart <- mult.chart(pls.lwd, type = "t2", alpha = alpha[k], main = "T2 chart for PLS loadings")
    ind.T2 <- which(T2.chart$t2 > T2.chart$ucl)
    if(length(ind.T2)<5){
      ind.T2 <- sort(T2.chart$t2, decreasing=TRUE, index.return=TRUE)$ix[1:5]
    }
    Xtr.T2 <- Xtr[,ind.T2]
    Xts.T2 <- Xts[,ind.T2]
    pls.cv.T2 <- plsr(ytr ~ Xtr.T2, ncomp=min(ncomp, ncol(Xtr.T2)-1), validation = "CV")
    opt.comp.T2 <- which.min(pls.cv.T2$validation$PRESS)
    opt.comp.T2 <- which.min(pls.cv.T2$validation$PRESS)
    opt.fit.T2  <- plsr(ytr ~ Xtr.T2, opt.comp.T2)
    Yhat.T2 <- opt.fit.T2$fitted.values[,,opt.comp.T2]
    opt.MSE.T2 <- rmsep(ytr, Yhat.T2) 
    Yts.hat.T2 <- predict(opt.fit.T2, Xts.T2)[,,opt.comp.T2]
    opt.MSE.ts.T2 <- rmsep(yts, Yts.hat.T2)
    R[k,]  <- c(opt.MSE.T2, opt.MSE.ts.T2, opt.comp.T2, length(ind.T2), alpha[k])
    V[[k]] <- ind.T2
  }
  ind.p <- which.min(R[,2])
  ind.v <- which.min(R[,4])  
  res <- list(mR=R[c(ind.p, ind.v),], mv=V[c(ind.p, ind.v)])
  res
}
