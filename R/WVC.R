#' Weighted Variable Contribution in PLS (WVC-PLS)
#'
#' @param y Vector of responses.
#' @param X Matrix of predictors.
#' @param ncomp Number of components.
#' @param normalize Divide WVC vectors by maximum value.
#' @param threshold Set loading weights smaller than threshold to 0 and recompute component.
#'
#' @return loading weights, loadings, regression coefficients, scores and Y-loadings
#' plus the WVC weights.
#' 
#' @description This implements the PLS-WVC2 component dependent version of WVC from Lin et al., i.e.,
#' using Equations 14, 16 and 19. The 
#' implementation is used in T. Mehmood, S. Sæbø, K.H. Liland, Comparison of variable selection methods in partial least
#' squares regression, Journal of Chemometrics 34 (2020) e3226. However, there is a mistake
#' in the notation in Mehmood et al. exchanging the denominator of Equation 19 (w'X'Xw) with (w'X'Yw).
#' 
#' @references Variable selection in partial least squares with the weighted variable contribution to the first singular value of the covariance matrix,
#' Weilu Lin, Haifeng Hang, Yingping Zhuang, Siliang Zhang, Chemometrics and Intelligent Laboratory Systems 183 (2018) 113–121.
#'
#' @examples
#' library(pls)
#' data(mayonnaise, package = "pls")
#' wvc <- WVC_pls(factor(mayonnaise$oil.type), mayonnaise$NIR, 10)
#' wvcNT <- WVC_pls(factor(mayonnaise$oil.type), mayonnaise$NIR, 10, TRUE, 0.5)
#' old.par <- par(mfrow=c(3,1), mar=c(2,4,1,1))
#' matplot(t(mayonnaise$NIR), type='l', col=1, ylab='intensity')
#' matplot(wvc$W[,1:3], type='l', ylab='W')
#' matplot(wvcNT$W[,1:3], type='l', ylab='W, thr.=0.5')
#' par(old.par)
#' 
#' @export
WVC_pls <- function(y, X, ncomp, normalize=FALSE, threshold=NULL){
  ncomp <- min(ncomp, c(ncol(X)-1))
  
  if (!is.factor(y)) {
    y <- as.matrix((as.numeric(y)))
    X <- as.matrix(scale(X))
  } else {
    y <- model.matrix(~y-1,data.frame(y=factor(y)))
    X <- as.matrix(X)
  }
  
  nc <- ncol(X)
  nr <- nrow(X)
  nresp <- ncol(y)
  
  W <- P <- WVC <- matrix(NA, nc, ncomp)
  B <- array(NA, c(nc+1, nresp, ncomp))
  S <- matrix(NA, nr, ncomp)
  Q <- EE <- SS <- matrix(NA, nresp, ncomp)
  
  for(a in 1:ncomp){
    usv <- svd(crossprod(X,y), 1,0)
    w  <- usv$u
    s  <- X%*%w
    p  <- crossprod(X,s) /c(crossprod(s))
    q  <- crossprod(y,s) /c(crossprod(s))
    
    EE[a] <- t(X%*%w)%*%y%*%q / (t(X%*%w)%*%X%*%w ) # alpha (Eq. 19)
    SS[a] <- t(X%*%w)%*%y%*%q                       # s (Eq. 14)
    W[,a] <- w # loading weights
    for(j in 1:nrow(W)){
      if(a==1){
        WVC[j,a] <- sqrt(ncol(X) *EE[1:a]%*% c(W[j,1:a])^2 %*%t(SS[1:a])/ # m * alpha *w^2 * s
                           c(EE[1:a]%*%t(SS[1:a])))                       # alpha * s (Eq. 16)
      } else{
        WVC[j,a] <- sqrt(ncol(X) *EE[1:a]%*% diag(W[j,1:a])^2 %*%SS[1:a]/c(EE[1:a]%*%SS[1:a]))
      }
    }
    if(normalize){
      WVC[,a] <- WVC[,a]/max(WVC[,a])
    }
    if(!is.null(threshold)){
      w[WVC[,a]<threshold] <- 0
      s  <- X%*%w
      p  <- crossprod(X,s) /c(crossprod(s))
      q  <- crossprod(y,s) /c(crossprod(s))
    }
    
    W[,a] <- w # loading weights
    S[,a] <- s # Score vector
    P[,a] <- p # X-loadings
    Q[,a] <- q # Y-loadings
    X <- X - s%*% t(p)
    y <- y - s%*% t(q)
    
    mat <- crossprod(P[,1:a], W[,1:a])
    id  <- subset(as.data.frame(which(cor(mat) > 0.9999, arr.ind=TRUE)), row < col)
    if(nrow(id) >0) mat <- mat[,-unique(id$row)]
    for(j in 1:nresp){
      BB <- W[,1:a]%*% solve(mat)%*%Q[j, 1:a]
      if(dim(id)[1] >0) {
        for(i in 1:length(id)){
          BB <- append(BB, 0, after = id[i])
        }
      }
      AA <- mean(y)-apply(X, 2, mean)%*%BB
      B[,j,a] <- c(AA, BB)
    }
  }
  
  if(is.null(colnames(X))){
    Xname <- as.character(1:nc)
  } else {
    Xname <- colnames(X)
  }
  rownames(WVC) <- rownames(W) <- rownames(P) <- Xname
  rownames(B)   <- c("constant", Xname)
  
  res <- list(W=W, P=P, B=B, S=S, Q=Q, WVC=WVC)
  return(res)
}
