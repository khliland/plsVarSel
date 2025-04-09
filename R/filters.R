####################
## Filter methods ##
####################

#' @aliases SR sMC LW RC
#' @title Filter methods for variable selection with Partial Least Squares.
#'
#' @description Various filter methods extracting and using information from 
#' \code{mvr} objects to assign importance to all included variables. Available 
#' methods are Significance Multivariate Correlation (sMC), Selectivity Ratio (SR), 
#' Variable Importance in Projections (VIP), Loading Weights (LW), Regression Coefficients (RC).
#'
#' @param pls.object \code{mvr} object from PLS regression.
#' @param opt.comp optimal number of components of PLS model.
#' @param p number of variables in PLS model.
#' @param X data matrix used as predictors in PLS modelling.
#' @param alpha_mc quantile significance for automatic selection of variables in \code{sMC}.
#' @param nsel number of variables to select.
#'
#' @return A vector having the same length as the number of variables in the associated
#' PLS model. High values are associated with high importance, explained variance or
#' relevance to the model.
#' 
#' The sMC has an attribute "quantile", which is the associated quantile of the
#' F-distribution, which can be used as a cut-off for significant variables, similar
#' to the cut-off of 1 associated with the VIP.
#' 
#' @details From plsVarSel 0.9.10, the VIP method handles multiple responses
#' correctly, as does the LW method. All other filter methods implemented in 
#' this package assume a single response and will give its results based on the
#' first response in multi-response cases.
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
#' data(gasoline, package = "pls")
#' library(pls)
#' pls  <- plsr(octane ~ NIR, ncomp = 10, validation = "LOO", data = gasoline)
#' comp <- which.min(pls$validation$PRESS)
#' X    <- unclass(gasoline$NIR)
#' vip <- VIP(pls, comp)
#' sr  <- SR (pls, comp, X)
#' smc <- sMC(pls, comp, X)
#' lw  <- LW (pls, comp)
#' rc  <- RC (pls, comp)
#' urc <- URC(pls, comp)
#' frc <- FRC(pls, comp)
#' mrm <- mRMR(pls, 401, X)$score
#' matplot(scale(cbind(vip, sr, smc, lw, rc, urc, frc, mrm)), type = 'l')
#'
#' @export
VIP <- function (pls.object, opt.comp, p = dim(pls.object$coef)[1]) {
  # Variable importance in prediction
  W <- pls.object$loading.weights
  WW <- W * W/apply(W, 2, function(x) sum(x * x))
  Q <- pls.object$Yloadings
  TT <- pls.object$scores
  if(length(dim(Q)) == 0){
    Q2 <- as.numeric(Q) * as.numeric(Q)
  } else {
    Q2 <- rowSums(t(Q * Q))
  }
  Q2TT <- Q2[1:opt.comp] * diag(crossprod(TT))[1:opt.comp]
  VIP <- sqrt(p * apply(sweep(WW[, 1:opt.comp, drop=FALSE],2,Q2TT,"*"), 1, sum)/sum(Q2TT))
  VIP
}


#' @rdname VIP
#' @export
SR <- function(pls.object, opt.comp, X){
  # Selectivity ratio
  X   <- unclass(as.matrix(X))
  RC  <- pls.object$coefficients[,1,opt.comp]
  Wtp <- RC/c(sqrt(crossprod(RC))[1])
  Ttp <- X %*% Wtp
  Ptp <- crossprod(X, Ttp)/c(crossprod(Ttp)[1])
  
  Xtp <- tcrossprod(Ttp, Ptp)
  Xr  <- X - Xtp
  SR  <- colSums(Xtp*Xtp)/colSums(Xr*Xr)
  #  var.test(Xtp, Xr)
  SR
}

#' @rdname VIP
#' @export
sMC <- function(pls.object, opt.comp, X, alpha_mc = 0.05){
  # Significance Multivariate Correlation
  # [smcF smcFcrit SSCregression SSResidual] = smc(b, X)
  # Output:
  # smcF : SMC F-values for the list of variables
  # smcFcrit: F-critical cutoff threshold value for significant important variables (smcF>smcFcrit)
  #
  # In case of publication of any application of this method,
  # please, cite the original work:
  # T.N. Tran*, N.L. Afanador, L.M.C. Buydens, L. Blanchet, 
  # Interpretation of variable importance in Partial Least Squares with Significance Multivariate Correlation (sMC), 
  # Chemometrics and Intelligent Laboratory Systems, Volume 138, 15 November 2014, Pages 153-160
  # DOI: http://dx.doi.org/10.1016/j.chemolab.2014.08.005
  
  b  <- pls.object$coefficients[,1,opt.comp]
  X   <- unclass(as.matrix(X))
  
  n <- dim(X)[1]
  
  yhat <- X%*%b
  Xhat <- tcrossprod(yhat,b)/crossprod(b)[1]
  Xresidual <- X - Xhat
  
  SSCregression <- colSums(Xhat^2)
  SSResidual    <- colSums(Xresidual^2)
  
  MSCregression <- SSCregression # 1 degrees of freedom
  MSResidual    <- SSResidual/(n-2)
  
  smcF     <- MSCregression/MSResidual;
  smcFcrit <- qf(1-alpha_mc,1,n-2)
#  list(smcF=smcF, smcFcrit=smcFcrit)
  attr(smcF, "quantile") <- smcFcrit
  smcF
}

#' @rdname VIP
#' @export
LW <- function(pls.object, opt.comp)
  # Loading weights
  abs(pls.object$loading.weights[ , opt.comp])

#' @rdname VIP
#' @export
RC <- function(pls.object, opt.comp)
  # Regression coefficients
  abs(pls.object$coef[ , 1, opt.comp])

# Remove names and sort vector
simplify <- function(X){
  names(X) <- NULL
  sort(X)
}

#' @rdname VIP
#' @export
URC <- function(pls.object, opt.comp){
  RC  <- pls.object$coefficients[ , 1, opt.comp]
  return(as.numeric(abs(RC) / max(abs(RC))))
}

#' @rdname VIP
#' @export
FRC <- function(pls.object, opt.comp){
  RC  <- pls.object$coefficients[ , 1, opt.comp]
  URC <- abs(RC) / max(abs(RC))
  return(URC / pls.object$validation$PRESS[opt.comp])
}

#' @rdname VIP
#' @importFrom praznik MRMR
#' @export
mRMR <- function(pls.object, nsel, X){
  MRMR(data.frame(X), pls.object$Yscores[,1], k = nsel, threads = 0)
}



###############
## filterPLS ##
###############

#' @title Optimisation of filters for Partial Least Squares
#'
#' @description Extract the index of influential variables based on threshold defiend for
#' LW (loading weights), RC (regression coef), JT (jackknife testing) and VIP (variable 
#' importance on projection).
#'
#' @param y vector of response values (\code{numeric} or \code{factor}).
#' @param X numeric predictor \code{matrix}.
#' @param ncomp integer number of components (default = 10).
#' @param ncomp.opt use the number of components corresponding to minimum error (minimum)
#'  or \code{ncomp} (same).
#' @param validation type of validation in the PLS modelling (default = "LOO").
#' @param LW.threshold threshold for Loading Weights if applied (default = NULL).
#' @param RC.threshold threshold for Regression Coefficients if applied (default = NULL).
#' @param URC.threshold threshold for Unit normalized Regression Coefficients if applied (default = NULL).
#' @param FRC.threshold threshold for Fitness normalized Regression Coefficients if applied (default = NULL).
#' @param JT.threshold threshold for Jackknife Testing if applied (default = NULL).
#' @param VIP.threshold threshold for Variable Importance on Projections if applied (default = NULL).
#' @param SR.threshold threshold for Selectivity Ration if applied (default = NULL).
#' @param sMC.threshold threshold for Significance Multivariate Correlation if applied (default = NULL).
#' @param mRMR.threshold threshold for minimum Redundancy Maximum Releveance if applied (default = NULL).
#' @param WVC.threshold threshold for Weighted Variable Contribution if applied (default = NULL).
#' @param ... additional paramters for \code{pls}, e.g. segmentation or similar.
#'
#' @details Filter methods are applied for variable selection with PLSR. This function can 
#' return selected variables and Root Mean Squared Error of Cross-Validation for various 
#' filter methods and determine optimum numbers of components.
#' 
#' @return Returns a list of lists containing filters (outer list), their selected variables,
#' optimal numbers of components and prediction accuracies.
#'
#' @author Tahir Mehmood, Kristian Hovde Liland, Solve Sæbø.
#'
#' @references T. Mehmood, K.H. Liland, L. Snipen, S. Sæbø, A review of variable selection 
#' methods in Partial Least Squares Regression, Chemometrics and Intelligent Laboratory Systems
#' 118 (2012) 62-69.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC/URC/FRC/mRMR), \code{\link{filterPLSR}}, \code{\link{spa_pls}}, 
#' \code{\link{stpls}}, \code{\link{truncation}}, \code{\link{bve_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{ipw_pls}}, \code{\link{ga_pls}}, \code{\link{rep_pls}}, \code{\link{WVC_pls}}, \code{\link{T2_pls}}.
#'
#' @examples
#' data(gasoline, package = "pls")
#' \dontrun{
#' with( gasoline, filterPLSR(octane, NIR, ncomp = 10, "minimum", validation = "LOO",
#'  RC.threshold = c(0.1,0.5), SR.threshold = 0.5))
#' }
#' 
#' @import pls
#' @export
filterPLSR <- function(y, X, ncomp = 10, ncomp.opt = c("minimum","same"), validation = "LOO", 
                       LW.threshold = NULL, RC.threshold = NULL, URC.threshold = NULL, FRC.threshold = NULL,
                       JT.threshold = NULL, VIP.threshold = NULL, SR.threshold = NULL, sMC.threshold = NULL,
                       mRMR.threshold = NULL, WVC.threshold = NULL, ...){
  
  # Strip X
  X <- unclass(as.matrix(X))
  
  n <- dim(X)[1]
  p <- dim(X)[2]

  modeltype <- "prediction"
  if (is.factor(y)) {
    modeltype <- "classification"
    y.orig <- as.numeric(y)
    y      <- model.matrix(~ y-1)
    tb     <- as.numeric(names(table(y)))
  }
  
  ncomp.opt <- ncomp.opt[1]
  
  # Primary PLSR fitting
  if(!is.null(JT.threshold)){
    pls.object <- plsr(y ~ X, ncomp = ncomp, validation = validation, jackknife = TRUE, ...)
  } else {
    pls.object <- plsr(y ~ X, ncomp = ncomp, validation = validation, ...)
  }
  if (modeltype == "prediction"){
    opt.comp <- which.min(pls.object$validation$PRESS[1,])
  } else if (modeltype == "classification"){
    classes <- lda_from_pls_cv(pls.object, X, y.orig, ncomp)
    opt.comp <- which.max(colSums(classes==y.orig))
  }

  # Apply filter methods
  selections <- list()
  if(!is.null(LW.threshold)){
    LWvalues <- LW(pls.object,opt.comp)
    if(is.logical(LW.threshold))
      LW.threshold <- 0.05
    if(length(LW.threshold) == 1){ # Single threshold
      selections$LW <- simplify(which(LWvalues > LW.threshold))
      if(length(selections$LW) < ncomp)
        selections$LW <- simplify(order(LWvalues, decreasing = TRUE)[1:ncomp])
    } else { # Optimise over several thresholds
      selections$LW <- list()
      selections$LW$comps  <- numeric(length(LW.threshold))
      selections$LW$RMSECV <- numeric(length(LW.threshold))
      names(selections$LW$comps)  <- LW.threshold
      names(selections$LW$RMSECV) <- LW.threshold
      for(i in 1:length(LW.threshold)){
        selections$LW[[i+2]] <- simplify(which(LWvalues > LW.threshold[i]))
        if(length(selections$LW[[i+2]]) < ncomp)
          selections$LW[[i+2]] <- simplify(order(LWvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$LW[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        if (modeltype == "prediction"){
          comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
          selections$LW$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
          selections$LW$comps[i]  <- comp
        } else if (modeltype == "classification"){
          classes <- lda_from_pls_cv(pls.this, X[,selections$LW[[i+2]], drop=FALSE], y.orig, ncomp)
          err  <- colSums(classes!=y.orig)
          comp <- ifelse(ncomp.opt == "minimum", which.min(err), ncomp)
          selections$LW$RMSECV[i] <- err[comp]/n
          selections$LW$comps[i]  <- comp
        }
      }
    }
  }
  
  if(!is.null(RC.threshold)){
    RCvalues <-  RC(pls.object,opt.comp)
    if(is.logical(RC.threshold))
      RC.threshold <- 0.01
    if(length(RC.threshold) == 1){ # Single threshold
      selections$RC <- simplify(which(RCvalues > RC.threshold))
      if(length(selections$RC) < ncomp)
        selections$RC <- simplify(order(RCvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$RC <- list()
      selections$RC$comps  <- numeric(length(RC.threshold))
      selections$RC$RMSECV <- numeric(length(RC.threshold))
      names(selections$RC$comps)  <- RC.threshold
      names(selections$RC$RMSECV) <- RC.threshold
      for(i in 1:length(RC.threshold)){
        selections$RC[[i+2]] <- simplify(which(RCvalues > RC.threshold[i]))
        if(length(selections$RC[[i+2]]) < ncomp)
          selections$RC[[i+2]] <- simplify(order(RCvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$RC[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        if (modeltype == "prediction"){
          comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
          selections$RC$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
          selections$RC$comps[i]  <- comp
        } else if (modeltype == "classification"){
          classes <- lda_from_pls_cv(pls.this, X[,selections$RC[[i+2]], drop=FALSE], y.orig, ncomp)
          err  <- colSums(classes!=y.orig)
          comp <- ifelse(ncomp.opt == "minimum", which.min(err), ncomp)
          selections$RC$RMSECV[i] <- err[comp]/n
          selections$RC$comps[i]  <- comp
        }
      }
    }
  }
  
  if(!is.null(URC.threshold)){
    URCvalues <-  URC(pls.object,opt.comp)
    if(is.logical(URC.threshold))
      URC.threshold <- 0.01
    if(length(URC.threshold) == 1){ # Single threshold
      selections$URC <- simplify(which(URCvalues > URC.threshold))
      if(length(selections$URC) < ncomp)
        selections$URC <- simplify(order(URCvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$URC <- list()
      selections$URC$comps  <- numeric(length(URC.threshold))
      selections$URC$RMSECV <- numeric(length(URC.threshold))
      names(selections$URC$comps)  <- URC.threshold
      names(selections$URC$RMSECV) <- URC.threshold
      for(i in 1:length(URC.threshold)){
        selections$URC[[i+2]] <- simplify(which(URCvalues > URC.threshold[i]))
        if(length(selections$URC[[i+2]]) < ncomp)
          selections$URC[[i+2]] <- simplify(order(URCvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$URC[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        if (modeltype == "prediction"){
          comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
          selections$URC$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
          selections$URC$comps[i]  <- comp
        } else if (modeltype == "classification"){
          classes <- lda_from_pls_cv(pls.this, X[,selections$URC[[i+2]], drop=FALSE], y.orig, ncomp)
          err  <- colSums(classes!=y.orig)
          comp <- ifelse(ncomp.opt == "minimum", which.min(err), ncomp)
          selections$URC$RMSECV[i] <- err[comp]/n
          selections$URC$comps[i]  <- comp
        }
      }
    }
  }

  if(!is.null(FRC.threshold)){
    FRCvalues <-  FRC(pls.object,opt.comp)
    if(is.logical(FRC.threshold))
      FRC.threshold <- 0.01
    if(length(FRC.threshold) == 1){ # Single threshold
      selections$FRC <- simplify(which(FRCvalues > FRC.threshold))
      if(length(selections$FRC) < ncomp)
        selections$FRC <- simplify(order(FRCvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$FRC <- list()
      selections$FRC$comps  <- numeric(length(FRC.threshold))
      selections$FRC$RMSECV <- numeric(length(FRC.threshold))
      names(selections$FRC$comps)  <- FRC.threshold
      names(selections$FRC$RMSECV) <- FRC.threshold
      for(i in 1:length(FRC.threshold)){
        selections$FRC[[i+2]] <- simplify(which(FRCvalues > FRC.threshold[i]))
        if(length(selections$FRC[[i+2]]) < ncomp)
          selections$FRC[[i+2]] <- simplify(order(FRCvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$FRC[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        if (modeltype == "prediction"){
          comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
          selections$FRC$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
          selections$FRC$comps[i]  <- comp
        } else if (modeltype == "classification"){
          classes <- lda_from_pls_cv(pls.this, X[,selections$FRC[[i+2]], drop=FALSE], y.orig, ncomp)
          err  <- colSums(classes!=y.orig)
          comp <- ifelse(ncomp.opt == "minimum", which.min(err), ncomp)
          selections$FRC$RMSECV[i] <- err[comp]/n
          selections$FRC$comps[i]  <- comp
        }
      }
    }
  }

  if(!is.null(VIP.threshold)){
    VIPvalues <- VIP(pls.object, opt.comp)
    if(is.logical(VIP.threshold))
      VIP.threshold <- 1
    if(length(VIP.threshold) == 1){ # Single threshold
      selections$VIP <- simplify(which(VIPvalues > VIP.threshold))
      if(length(selections$VIP) < ncomp)
        selections$VIP <- simplify(order(VIPvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$VIP <- list()
      selections$VIP$comps  <- numeric(length(VIP.threshold))
      selections$VIP$RMSECV <- numeric(length(VIP.threshold))
      names(selections$VIP$comps)  <- VIP.threshold
      names(selections$VIP$RMSECV) <- VIP.threshold
      for(i in 1:length(VIP.threshold)){
        selections$VIP[[i+2]] <- simplify(which(VIPvalues > VIP.threshold[i]))
        if(length(selections$VIP[[i+2]]) < ncomp)
          selections$VIP[[i+2]] <- simplify(order(VIPvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$VIP[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        if (modeltype == "prediction"){
          comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
          selections$VIP$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
          selections$VIP$comps[i]  <- comp
        } else if (modeltype == "classification"){
          classes <- lda_from_pls_cv(pls.this, X[,selections$VIP[[i+2]], drop=FALSE], y.orig, ncomp)
          err  <- colSums(classes!=y.orig)
          comp <- ifelse(ncomp.opt == "minimum", which.min(err), ncomp)
          selections$VIP$RMSECV[i] <- err[comp]/n
          selections$VIP$comps[i]  <- comp
        }
      }
    }
  }
  
  if(!is.null(SR.threshold)){
    SRvalues <- SR(pls.object, opt.comp, X)
    if(is.logical(SR.threshold))
      SR.threshold <- pf(0.99, df1 = n-1, df2 = n-2)
    if(length(SR.threshold) == 1){ # Single threshold
      selections$SR <- simplify(which(SRvalues > SR.threshold))
      if(length(selections$SR) < ncomp)
        selections$SR <- simplify(order(SRvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$SR <- list()
      selections$SR$comps  <- numeric(length(SR.threshold))
      selections$SR$RMSECV <- numeric(length(SR.threshold))
      names(selections$SR$comps)  <- SR.threshold
      names(selections$SR$RMSECV) <- SR.threshold
      for(i in 1:length(SR.threshold)){
        selections$SR[[i+2]] <- simplify(which(SRvalues > SR.threshold[i]))
        if(length(selections$SR[[i+2]]) < ncomp)
          selections$SR[[i+2]] <- simplify(order(SRvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$SR[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        if (modeltype == "prediction"){
          comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
          selections$SR$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
          selections$SR$comps[i]  <- comp
        } else if (modeltype == "classification"){
          classes <- lda_from_pls_cv(pls.this, X[,selections$SR[[i+2]], drop=FALSE], y.orig, ncomp)
          err  <- colSums(classes!=y.orig)
          comp <- ifelse(ncomp.opt == "minimum", which.min(err), ncomp)
          selections$SR$RMSECV[i] <- err[comp]/n
          selections$SR$comps[i]  <- comp
        }
      }
    }
  }
  
  if(!is.null(sMC.threshold)){
    sMCvalues <- sMC(pls.object, opt.comp, X)
    if(is.logical(sMC.threshold))
      sMC.threshold <- pf(0.99, df1 = n-1, df2 = n-2)
    if(length(sMC.threshold) == 1){ # Single threshold
      selections$sMC <- simplify(which(sMCvalues > sMC.threshold))
      if(length(selections$sMC) < ncomp)
        selections$sMC <- simplify(order(sMCvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$sMC <- list()
      selections$sMC$comps  <- numeric(length(sMC.threshold))
      selections$sMC$RMSECV <- numeric(length(sMC.threshold))
      names(selections$sMC$comps)  <- sMC.threshold
      names(selections$sMC$RMSECV) <- sMC.threshold
      for(i in 1:length(sMC.threshold)){
        selections$sMC[[i+2]] <- simplify(which(sMCvalues > sMC.threshold[i]))
        if(length(selections$sMC[[i+2]]) < ncomp)
          selections$sMC[[i+2]] <- simplify(order(sMCvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$sMC[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        if (modeltype == "prediction"){
          comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
          selections$sMC$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
          selections$sMC$comps[i]  <- comp
        } else if (modeltype == "classification"){
          classes <- lda_from_pls_cv(pls.this, X[,selections$sMC[[i+2]], drop=FALSE], y.orig, ncomp)
          err  <- colSums(classes!=y.orig)
          comp <- ifelse(ncomp.opt == "minimum", which.min(err), ncomp)
          selections$sMC$RMSECV[i] <- err[comp]/n
          selections$sMC$comps[i]  <- comp
        }
      }
    }
  }
  
  if(!is.null(JT.threshold)){
    JTvalues <- jack.test(pls.object, ncomp = opt.comp)$pvalues[,1,1]
    if(is.logical(JT.threshold))
      JT.threshold <- 0.05
    if(length(JT.threshold) == 1){ # Single threshold
      selections$JT <- simplify(which(JTvalues > JT.threshold))
      if(length(selections$JT) < ncomp)
        selections$JT <- simplify(order(JTvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$JT <- list()
      selections$JT$comps  <- numeric(length(JT.threshold))
      selections$JT$RMSECV <- numeric(length(JT.threshold))
      names(selections$JT$comps)  <- JT.threshold
      names(selections$JT$RMSECV) <- JT.threshold
      for(i in 1:length(JT.threshold)){
        selections$JT[[i+2]] <- simplify(which(JTvalues > JT.threshold[i]))
        if(length(selections$JT[[i+2]]) < ncomp)
          selections$JT[[i+2]] <- simplify(order(JTvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$JT[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        if (modeltype == "prediction"){
          comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
          selections$JT$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
          selections$JT$comps[i]  <- comp
        } else if (modeltype == "classification"){
          classes <- lda_from_pls_cv(pls.this, X[,selections$JT[[i+2]], drop=FALSE], y.orig, ncomp)
          err  <- colSums(classes!=y.orig)
          comp <- ifelse(ncomp.opt == "minimum", which.min(err), ncomp)
          selections$JT$RMSECV[i] <- err[comp]/n
          selections$JT$comps[i]  <- comp
        }
      }
    }
  }
  
  if(!is.null(mRMR.threshold)){
    mRMRvalues <-  mRMR(pls.object,ncol(X),X)$selection
    if(is.logical(mRMR.threshold))
      mRMR.threshold <- 1
    if(length(mRMR.threshold) == 1){ # Single threshold
      selections$mRMR <- simplify(mRMRvalues[1:mRMR.threshold])#which(mRMRvalues > mRMR.threshold))
      if(length(selections$mRMR) < ncomp)
        selections$mRMR <- simplify(order(mRMRvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$mRMR <- list()
      selections$mRMR$comps  <- numeric(length(mRMR.threshold))
      selections$mRMR$RMSECV <- numeric(length(mRMR.threshold))
      names(selections$mRMR$comps)  <- mRMR.threshold
      names(selections$mRMR$RMSECV) <- mRMR.threshold
      for(i in 1:length(mRMR.threshold)){
        selections$mRMR[[i+2]] <- simplify(mRMRvalues[1:mRMR.threshold])#which(mRMRvalues > mRMR.threshold[i]))
        if(length(selections$mRMR[[i+2]]) < ncomp)
          selections$mRMR[[i+2]] <- simplify(mRMRvalues[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$mRMR[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        if (modeltype == "prediction"){
          comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
          selections$mRMR$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
          selections$mRMR$comps[i]  <- comp
        } else if (modeltype == "classification"){
          classes <- lda_from_pls_cv(pls.this, X[,selections$mRMR[[i+2]], drop=FALSE], y.orig, ncomp)
          err  <- colSums(classes!=y.orig)
          comp <- ifelse(ncomp.opt == "minimum", which.min(err), ncomp)
          selections$mRMR$RMSECV[i] <- err[comp]/n
          selections$mRMR$comps[i]  <- comp
        }
      }
    }
  }

  if(!is.null(WVC.threshold)){
    WVCvalues <-  abs(WVC_pls(y,X,ncomp,TRUE)$B[-1,1,opt.comp])
    if(is.logical(WVC.threshold))
      WVC.threshold <- 0.9
    if(length(WVC.threshold) == 1){ # Single threshold
      selections$WVC <- simplify(which(WVCvalues > WVC.threshold))
      if(length(selections$WVC) < ncomp)
        selections$WVC <- simplify(order(WVCvalues, decreasing = TRUE)[1:ncomp])
    } else {
      selections$WVC <- list()
      selections$WVC$comps  <- numeric(length(WVC.threshold))
      selections$WVC$RMSECV <- numeric(length(WVC.threshold))
      names(selections$WVC$comps)  <- WVC.threshold
      names(selections$WVC$RMSECV) <- WVC.threshold
      for(i in 1:length(WVC.threshold)){
        selections$WVC[[i+2]] <- simplify(which(WVCvalues > WVC.threshold[i]))
        if(length(selections$WVC[[i+2]]) < ncomp)
          selections$WVC[[i+2]] <- simplify(order(WVCvalues, decreasing = TRUE)[1:ncomp])
        pls.this <- plsr(y ~ X[,selections$WVC[[i+2]], drop=FALSE], ncomp = ncomp, validation = validation, ...)
        if (modeltype == "prediction"){
          comp <- ifelse(ncomp.opt == "minimum", which.min(pls.this$valid$PRESS[1,]), ncomp)
          selections$WVC$RMSECV[i] <- RMSEP(pls.this, estimate="CV")$val[1,1,comp+1]
          selections$WVC$comps[i]  <- comp
        } else if (modeltype == "classification"){
          classes <- lda_from_pls_cv(pls.this, X[,selections$WVC[[i+2]], drop=FALSE], y.orig, ncomp)
          err  <- colSums(classes!=y.orig)
          comp <- ifelse(ncomp.opt == "minimum", which.min(err), ncomp)
          selections$WVC$RMSECV[i] <- err[comp]/n
          selections$WVC$comps[i]  <- comp
        }
      }
    }
  }
  
  selections
}

# Previous version
# SR <- function(pls.object, opt.comp, X){
#   # Selectivity ratio
#   X   <- as.matrix(X)
#   RC  <- pls.object$coefficients[,1,opt.comp]
#   Wtp <- RC/norm(matrix(RC))
#   Ttp <- X%*%Wtp
#   Ptp <- (t(X)%*%Ttp)/c((t(Ttp)%*%Ttp))
#   
#   Xtp <- Ttp%*%t(Ptp)
#   Xr  <- X-Xtp
#   SR  <- colSums(Xtp*Xtp)/colSums(Xr*Xr)
#   #  var.test(Xtp, Xr)
#   SR
# }
