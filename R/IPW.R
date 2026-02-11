#' @title Iterative predictor weighting PLS (IPW-PLS)
#' @aliases ipw_pls ipw_pls_legacy
#'
#' @description An iterative procedure for variable elimination.
#'
#' @param y vector of response values (\code{numeric} or \code{factor}).
#' @param X numeric predictor \code{matrix}.
#' @param ncomp integer number of components (default = 10).
#' @param no.iter the number of iterations (default = 10).
#' @param IPW.threshold threshold for regression coefficients (default = 0.1).
#' @param filter which filtering method to use (among "RC", "SR", "LW", "VIP", "sMC")
#' @param scale standardize data (default=TRUE, as in reference)
#'
#' @details This is an iterative elimination procedure where a measure of predictor 
#' importance is computed after fitting a PLSR model (with complexity chosen based
#' on predictive performance). The importance measure is used both to re-scale the 
#' original X-variables and to eliminate the least important variables before
#' subsequent model re-fitting
#' 
#' The IPW implementation was corrected in \code{plsVarSel} version 0.9.5. For backward
#' compatibility the old implementation is included as \code{ipw_pls_legacy}.
#'  
#' @return Returns a vector of variable numbers corresponding to the model 
#' having lowest prediction error.
#'
#' @author Kristian Hovde Liland
#'
#' @references M. Forina, C. Casolino, C. Pizarro Millan, Iterative predictor weighting
#'  (IPW) PLS: a technique for the elimination of useless predictors in regression problems,
#'  Journal of Chemometrics 13 (1999) 165-184.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{shaving}}, 
#' \code{\link{stpls}}, \code{\link{truncation}},
#' \code{\link{bve_pls}}, \code{\link{ga_pls}}, \code{\link{ipw_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{rep_pls}}, \code{\link{spa_pls}},
#' \code{\link{lda_from_pls}}, \code{\link{setDA}}.
#' 
#' @examples
#' data(gasoline, package = "pls")
#' with( gasoline, ipw_pls(octane, NIR) )
#'
#' @export
ipw_pls <- function(y, X, ncomp=10, no.iter=10, IPW.threshold=0.01, filter="RC", scale=TRUE){
  
  # Strip X
  X <- unclass(as.matrix(X))
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  if(is.factor(y)) {
    modeltype <- "classification"
    y.orig <- y
    y <- model.matrix(~y-1,data.frame(y=y))
  } else {
    modeltype <- "prediction"
    if(scale){
      y <- scale(y)
    }
  }
  Xscaled <- scale(X)
  s <- attr(Xscaled, 'scaled:scale') # Store standard deviations
  if(scale){ # Scale X by default
    Xorig <- Xscaled
  } else {
    Xorig <- X
  }
  z <- s/sum(s) # First scaling, before PLS
  early_stopping <- FALSE
  for(i in seq_len(no.iter)){
    # Rescaling
    X  <- Xorig*rep(z,each=n)
    
    # PLS modelling
    pls.object <- plsr(y ~ X, ncomp=ncomp, validation = "CV")
    if (modeltype == "prediction"){
      opt.comp <- which.min(pls.object$validation$PRESS[1,])
    } else if (modeltype == "classification"){
      classes <- lda_from_pls_cv(pls.object, X, y.orig, ncomp)
      opt.comp <- which.max(colSums(classes==y.orig))
    }
    
    # Filter calculation
    weights <- switch(filter,
                      RC  = RC(pls.object, opt.comp),
                      SR  = SR(pls.object, opt.comp, X),
                      LW  = LW(pls.object, opt.comp),
                      VIP = VIP(pls.object, opt.comp, p),
                      sMC = sMC(pls.object, opt.comp, X))
    weights[!is.finite(weights)] <- 0 # Correct for non-finite weights
    
    # Calculate importance
    z <- abs(weights)*s
    z <- z0 <- z/sum(z)
    z[z<IPW.threshold] <- 0
    if(sum(z)==0){
      z0[z0<max(z0)] <- 0
      z <- z0
      warning('The combination of parameters removed all variables, defaulting single variable.')
      early_stopping <- TRUE
      break()
    }
  }
  if(early_stopping){
    ipw.selection <- which.max(z)
  } else {
    ipw.selection <- simplify(which(z >= IPW.threshold))
  }
  return(list(ipw.selection=ipw.selection, ipw.importance=z))
}

#' @rdname ipw_pls
#' @export
ipw_pls_legacy <- function(y, X, ncomp=10, no.iter=10, IPW.threshold=0.1){
  
  # Strip X
  X <- unclass(as.matrix(X))
  
  if(is.factor(y)) {
    modeltype <- "classification"
    y.orig <- y
    y <- model.matrix(~y-1,data.frame(y=y))
  } else {
    modeltype <- "prediction"
    y <- scale(y)
  }
  #X<- scale(X)
  for(i in seq_len(no.iter)){
    pls.object <- plsr(y ~ X, ncomp=ncomp, validation = "CV")
    if (modeltype == "prediction"){
      opt.comp <- which.min(pls.object$validation$PRESS[1,])
    } else if (modeltype == "classification"){
      classes <- lda_from_pls_cv(pls.object, X, y.orig, ncomp)
      opt.comp <- which.max(colSums(classes==y.orig))
    }
    # Press    <- pls.object$valid$PRESS[1,]
    # opt.comp <- which.min(Press)
    # pls.fit  <- plsr(y ~ X, ncomp=opt.comp)
    # RC <- pls.fit$coef[,1,opt.comp]	
    RC <- pls.object$coef[,1,opt.comp]	
    X  <- X*RC
  }
  ipw.selection <- which(abs(RC) >= IPW.threshold)
  if(length(ipw.selection)<= (ncomp +1)){
    ipw.selection <- sort(RC,decreasing = TRUE, index.return = T)$ix [1:ncomp]
  }
  return(list(ipw.selection=simplify(ipw.selection)))
}
