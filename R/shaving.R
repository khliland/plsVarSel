#' @title Repeated shaving of variables
#'
#' @description One of five filter methods can be chosen for repeated shaving of
#' a certain percentage of the worst performing variables. Performance of the
#' reduced models are stored and viewable through \code{print} and \code{plot}
#' methods.
#'
#' @param y vector of response values (\code{numeric} or \code{factor}).
#' @param X numeric predictor \code{matrix}.
#' @param ncomp integer number of components (default = 10).
#' @param method filter method, i.e. SR, VIP, sMC, LW or RC given as \code{character}.
#' @param prop proportion of variables to be removed in each iteration (\code{numeric}).
#' @param min.left minimum number of remaining variables.
#' @param comp.type use number of components chosen by cross-validation, \code{"CV"},
#' or fixed, \code{"max"}.
#' @param validation type of validation for \code{plsr}. The default is "CV". If more
#' than one set of CV segments is wanted, use a vector of lenth two, e.g. \code{c("CV",5)}.
#' @param fixed vector of indeces for compulsory/fixed variables that should always be
#' included in the modelling.
#' @param newy validation response for RMSEP/error computations.
#' @param newX validation predictors for RMSEP/error computations.
#' @param segments see \code{mvr} for documentation of segment choices.
#' @param plsType Type of PLS model, "plsr" or "cppls".
#' @param Y.add Additional response for CPPLS, see \code{plsType}.
#' @param ... additional arguments for \code{plsr} or \code{cvsegments}.
#' @param x object of class \code{shaved} for plotting or printing.
#' @param what plot type. Default = "error". Alternative = "spectra".
#' @param index which iteration to plot. Default = "min"; corresponding to minimum RMSEP.
#' @param log logarithmic x (default) or y scale.
#'
#' @details Variables are first sorted with respect to some importancemeasure, 
#' and usually one of the filter measures described above are used. Secondly, a 
#' threshold is used to eliminate a subset of the least informative variables. Then
#' a model is fitted again to the remaining variables and performance is measured. 
#' The procedure is repeated until maximum model performance is achieved.
#'
#' @return Returns a list object of class \code{shaved} containing the method type,
#' the error, number of components, and number of variables per reduced model. It
#' also contains a list of all sets of reduced variable sets plus the original data.
#'
#' @author Kristian Hovde Liland
#'
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{shaving}}, 
#' \code{\link{stpls}}, \code{\link{truncation}},
#' \code{\link{bve_pls}}, \code{\link{ga_pls}}, \code{\link{ipw_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{rep_pls}}, \code{\link{spa_pls}},
#' \code{\link{lda_from_pls}}, \code{\link{lda_from_pls_cv}}, \code{\link{setDA}}.
#'  
#' @examples
#' data(mayonnaise, package = "pls")
#' sh <- shaving(mayonnaise$design[,1], pls::msc(mayonnaise$NIR), type = "interleaved")
#' pars <- par(mfrow = c(2,1), mar = c(4,4,1,1))
#' plot(sh)
#' plot(sh, what = "spectra")
#' par(pars)
#' print(sh)
#'
#' @importFrom progress progress_bar
#' @export
shaving <- function(y, X, ncomp = 10, method = c("SR", "VIP", "sMC", "LW", "RC"), prop = 0.2, min.left = 2,
                    comp.type = c("CV", "max"), validation = c("CV", 1), fixed = integer(0), newy = NULL, newX = NULL,
                    segments = 10, plsType = "plsr", Y.add = NULL, ...){
  
  modeltype <- "prediction"
  if (is.factor(y)) {
    modeltype <- "classification"
    y.orig <- as.numeric(y)
    y      <- model.matrix(~ y-1)
    tb     <- as.numeric(names(table(y)))
    min.left <- max(min.left, dim(y)[2])
  }

  # Initialization
  n <- nrow(X)
  p <- ncol(X)
  if(min.left < 1){ # Get minimum size
    min.left <- max(1, round(min.left * p))
  }
  comp.type <- comp.type[1]
  if(length(method) == 5){
    method <- method[1] # Default to SR
  }
  
  # Find the number of shaves
  nrep <- 1
  nfix <- length(fixed)
  if(nfix >= p)
    stop('Too many variables are fixed.')
  left <- p-nfix
  while(left > min.left){
    nrep <- nrep + 1
    left <- left - max(1, round(prop * left))
  }
  rems <- integer(nrep)
  left <- p-nfix
  for(i in 1:nrep){
    rems[i] <- left
    left <- left - max(1, round(prop * left))
  }
  rems <- -diff(rems)
  left.vec <- 1:p # Vector of remaining variables
  
  # Prepare for expanded/repated segments
  if(!is.null(segments) && !is.na(segments) && length(segments) == 1 && validation[1] == "CV"){
    rep.seg <- 1
    if(length(validation) == 2){
      rep.seg    <- validation[2]
      validation <- validation[1]
    }
    segs   <- list()
    for(i in 1:rep.seg)
      segs <- c(segs, cvsegments(n, segments, ...))
  } else {
    segs <- segments
  }
  
  # Prepare storage
  error       <- numeric(nrep)
  comps       <- integer(nrep)
  nvar        <- integer(nrep)
  variables   <- list()
  # left        <- p
  
  # Shaving loop
  pb <- progress_bar$new(total = nrep, format = "  [:bar] :percent ")
  for(i in 1:nrep){
    nvar[i] <- length(left.vec)
    variables[[i]] <- left.vec

    comp      <- min(n-1, nvar[i], ncomp)
    if(is.null(Y.add)){
      data      <- data.frame(y = y, X = I(X[, left.vec, drop = FALSE]))
    } else {
      data      <- data.frame(y = y, X = I(X[, left.vec, drop = FALSE]), Y.add = Y.add)
    }
    if(validation == "CV"){
      if(plsType == "cppls"){
        if(is.null(Y.add)){
          pls       <- cppls(y ~ X, ncomp=comp, data=data, validation = validation, segments = segs, ...)
          
        } else {
          pls       <- cppls(y ~ X, ncomp=comp, data=data, validation = validation, segments = segs, Y.add = Y.add, ...)
        }
      } else {
        pls       <- plsr(y ~ X, ncomp=comp, data=data, validation = validation, segments = segs, ...)
      }
    } else {
      if(plsType == "cppls"){
        if(is.null(Y.add)){
          pls       <- cppls(y ~ X, ncomp=comp, data=data, validation = validation, ...)
        } else {
          pls       <- cppls(y ~ X, ncomp=comp, data=data, validation = validation, Y.add = Y.add, ...)
        }
      } else {
        pls       <- pls(y ~ X, ncomp=comp, data=data, validation = validation, ...)
      }
    }
    # factor y
    if(comp.type == "CV"){
      if (modeltype == "prediction"){
        opt.comp <- which.min(pls$validation$PRESS[1,])
      } else if (modeltype == "classification"){
        if(is.null(Y.add)){
          classes <- lda_from_pls_cv(pls, data$X, y.orig, comp)
        } else {
          classes <- lda_from_pls_cv(pls, data$X, y.orig, comp, Y.add)
        }
        opt.comp <- which.max(colSums(classes==y.orig))
      }
    } else {
      opt.comp <- comp
    }
    comps[i]  <- opt.comp
    if(is.null(newX) || is.null(newy)){
      if (modeltype == "prediction"){
        error[i]  <- RMSEP(pls, estimate = "CV")$val[1,1, opt.comp + 1]
      } else if (modeltype == "classification"){
        error[i] <- sum(classes[,opt.comp] != y.orig)/n
      }
    } else {
      if (modeltype == "prediction"){
        val <- data.frame(X=I(newX[,left.vec, drop = FALSE]), y=newy)
        error[i]  <- RMSEP(pls, estimate = "test", newdata = val)$val[1,1, opt.comp + 1]
      } else if (modeltype == "classification"){
        yv.orig <- as.numeric(newy)
        yv      <- model.matrix(~ newy-1)
        val <- data.frame(X=I(newX[,left.vec, drop = FALSE]), y=newy)
        error[i] <- sum(lda_from_pls(pls, yv.orig, val$X, opt.comp)[,opt.comp] != yv.orig)/length(yv.orig)
      }
    }
    
    if(method == "SR"){
      weights <- SR(pls, opt.comp, data$X)
    }
    if(method == "VIP"){
      weights <- VIP(pls, opt.comp, p = dim(pls$coef)[1])
    }
    if(method == "sMC"){
      weights <- sMC(pls, opt.comp, data$X, alpha_mc = 0.05)
    }
    if(method == "LW"){
      weights <- LW(pls, opt.comp)
    }
    if(method == "RC"){
      weights <- RC(pls, opt.comp)
    }
    if(i < nrep)
      rem <- rems[i]
    # rem       <- max(1, round(prop * (left-nfix)))
    left      <- left - rem
    if(nfix == 0){
      left.vec <- left.vec[sort(order(weights)[-(1:rem)])]
    } else {
      left.vec1 <- left.vec[!(left.vec%in%fixed)]
      weights1  <- weights[!(left.vec%in%fixed)]
      left.vec1 <- left.vec1[sort(order(weights1)[-(1:rem)])]
      left.vec  <- sort(c(left.vec1,fixed))
    }
    pb$tick()
  }
  obj <- list(min.error = min(error), min.red = index <- length(error)-which.min(rev(error)), method = method, error = error, comps = comps, nvar = nvar, variables = variables, X = X)
  class(obj) <- c("shaved", "list")
  obj
}

#' @rdname shaving
#' @export
plot.shaved <- function(x, y, what = c('error', 'spectra'), index = "min", log = "x", ...){
  if(length(what) > 1)
    what <- what[1]
  
  if(what == "error"){
    with(x, 
         plot(error ~ nvar, xlim = rev(range(x$nvar)), type = 'b', log = log, ...))
  } else { # "spectra"
    with(x, {
      # Prepare
      if(is.character(index) && index == "min"){
        index <- length(error)-which.min(rev(error))+1
      }
      X.NA <- X
      X.NA[, -variables[[index]]] <- NA
      x <- 1:ncol(X)
      x.NA <- x
      x.NA[-variables[[index]]] <- NA
      xaxt <- par("xaxt")
      if(!is.null(labels <- colnames(X))){
        xaxt <- "n"
      }
      X  <- t(X)
      xr <- range(X)
      
      # Plot
      matplot(x, X, type = 'l', col = 1, xaxt = xaxt, ylim = c(xr[1]-0.05*diff(xr), xr[2]), ...)
      matplot(x, t(X.NA), type = 'l', col = 2, add = TRUE, ...)
      points(x.NA,rep(xr[1]-0.05*diff(xr), length(x)), col = 'blue', pch = 20)
      if(xaxt == "n"){
        ticks <- axTicks(1)
        ticks <- ticks[ticks >= 1 & ticks <= length(labels)]
        axis(1, ticks, labels[ticks])
      }
    })
  }
}

#' @rdname shaving
#' @export
print.shaved <- function(x, ...){
  cat('Shaving using ', x$method, ':\n', sep = "")
  cat('Minimum error = ', (x$min.err), ', achieved after ', x$min.red,
      ' out of ', length(x$error)-1, ' reductions using ',
      x$nvar[length(x$error)-which.min(rev(x$error))+1], ' variables.', sep = "")
}

