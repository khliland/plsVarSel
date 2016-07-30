## PLS-DA + LDA
#' LDA/QDA classification from PLS model
#' 
#' For each number of components LDA/QDA models are created from the 
#' scores of the supplied PLS model and classifications are performed.
#'
#' @param model \code{pls} model fitted with the \code{pls} package
#' @param grouping vector of grouping labels
#' @param newdata predictors in the same format as in the \code{pls} model
#' @param ncomp maximum number of PLS components
#'
#' @return matrix of classifications
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{shaving}}, 
#' \code{\link{stpls}}, \code{\link{truncation}},
#' \code{\link{bve_pls}}, \code{\link{ga_pls}}, \code{\link{ipw_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{rep_pls}}, \code{\link{spa_pls}},
#' \code{\link{lda_from_pls}}, \code{\link{lda_from_pls_cv}}, \code{\link{setDA}}.
#' 
#' @examples
#' data(mayonnaise, package = "pls")
#' mayonnaise <- within(mayonnaise, {dummy <- model.matrix(~y-1,data.frame(y=factor(oil.type)))})
#' pls <- plsr(dummy ~ NIR, ncomp = 10, data = mayonnaise, subset = train)
#' with(mayonnaise, {
#'  classes <- lda_from_pls(pls, oil.type[train], NIR[!train,], 10)
#'  colSums(oil.type[!train] == classes) # Number of correctly classified out of 42
#' })
#' 
#' @export
lda_from_pls <- function(model, grouping, newdata, ncomp){
  LQ <- getplsEnv("LQ")
  
  if(LQ == "max"){
    labels  <- names(table(grouping))
    predVal <- predict(model, newdata = newdata, ncomp = 1:ncomp)
    class   <- apply(predVal,c(1,3),which.max)
    for(i in 1:ncol(class)){
      class[[i]]   <- labels[class[[i]]]
    }
    colnames(class) <- paste("Comp.", 1:ncomp, sep="")
    return(class)
    
  } else { # LDA or QDA
    # Extract and predict scores
    scoresCal <- scores(model)
    scoresVal <- predict(model, newdata = newdata, type = "scores")
    
    # Prepare for storage
    N <- dim(scoresVal)
    class <- matrix(0, N[1],ncomp)
    
    # Create ncomp lda models and predict classes
    for(i in 1:ncomp){
      if(getplsEnv("LQ") == "lda"){
        ldai <- lda(scoresCal[, 1:i, drop = FALSE], grouping, tol = 1.0e-10)
      }
      if(getplsEnv("LQ") == "qda"){
        ldai <- qda(scoresCal[, 1:i, drop = FALSE], grouping, tol = 1.0e-10)
      }
      class[, i] <- predict(ldai, scoresVal[, 1:i, drop = FALSE])$class
    }
    colnames(class) <- paste("Comp.", 1:ncomp, sep="")
    return(class)
  }
}

# Cross-validate PLS-DA + LDA (dirty code)
#' Cross-validated LDA/QDA classification from PLS model
#'
#' For each number of components LDA/QDA models are created from the 
#' scores of the supplied PLS model and classifications are performed.
#' This use of cross-validation has limitations. Handle with care!
#' 
#' @param model \code{pls} model fitted with the \code{pls} package
#' @param X predictors in the same format as in the \code{pls} model
#' @param y vector of grouping labels
#' @param ncomp maximum number of PLS components
#'
#' @return matrix of classifications
#'
#' @examples
#' data(mayonnaise, package = "pls")
#' mayonnaise <- within(mayonnaise, {dummy <- model.matrix(~y-1,data.frame(y=factor(oil.type)))})
#' pls <- plsr(dummy ~ NIR, ncomp = 15, data = mayonnaise, subset = train, 
#'             validation = "CV", segments = 40, segment.type = "consecutive")
#' with(mayonnaise, {
#'  classes <- lda_from_pls_cv(pls, NIR[train,], oil.type[train], 15)
#'  colSums(oil.type[train] == classes) # Number of correctly classified out of 120
#' })
#' 
#' @importFrom MASS lda qda
#' 
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{shaving}}, 
#' \code{\link{stpls}}, \code{\link{truncation}},
#' \code{\link{bve_pls}}, \code{\link{ga_pls}}, \code{\link{ipw_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{rep_pls}}, \code{\link{spa_pls}},
#' \code{\link{lda_from_pls}}, \code{\link{lda_from_pls_cv}}, \code{\link{setDA}}.
#' 
#' @export
lda_from_pls_cv <- function(model, X, y, ncomp){
  N <- dim(model$scores)
  ncomp <- min(min(min(ncomp, N[2]),dim(X)[2]),dim(y)[2])
  classes  <- matrix(0, N[1], ncomp)
  dummy    <- model.matrix(~ factor(y)-1)
  segments <- model$validation$segments # Extract segments
  data     <- data.frame(X = I(X), y = y, dummy = I(dummy))
  for(i in 1:length(segments)){
    # Update model with new data
    model_i <- update(model, subset = NULL, 
                      formula = dummy~X, 
                      data = data[-segments[[i]],,drop=FALSE], 
                      validation = "none", ncomp = ncomp)
    
    # Predict left out
    comp <- 1
    if(!is.null(compp <- dim(scores(model_i))[2]))
      comp <- min(compp,ncomp)
    classes[segments[[i]],1:comp] <- lda_from_pls(model_i, y[-segments[[i]]], data[segments[[i]],, drop = FALSE], comp)
  }
  colnames(classes) <- paste("Comp.", 1:ncomp, sep="")
  classes
}
