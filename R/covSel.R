#' Covariance Selection - CovSel
#' 
#' @description Sequential selection of variables based on squared covariance
#' with response and intermediate deflation (as in Partial Least Squares).
#'
#' @param X \code{matrix} of input variables
#' @param Y \code{matrix} of response variable(s)
#' @param nvar maximum number of variables
#'
#' @return
#' \item{selected}{an integer vector of selected variables}
#' \item{scores}{a matrix of score vectors}
#' \item{loadings}{a matrix of loading vectors}
#' \item{Yloadings}{a matrix of Y loadings}
#' @export
#' @references J.M. Roger, B. Palagos, D. Bertrand, E. Fernandez-Ahumada. CovSel: Variable selection for highly multivariate and multi-response calibration: Application to IR spectroscopy. Chemom Intel Lab Syst. 2011;106(2):216-223.
#' P. Mishra, A brief note on a new faster covariate's selection (fCovSel) algorithm, Journal of Chemometrics 36(5) 2022.
#'
#' @examples
#' data(gasoline, package = "pls")
#' sels <- with(gasoline, covSel(NIR, octane, 5))
#' matplot(t(gasoline$NIR), type = "l")
#' abline(v = sels$selected, col = 2)
covSel <- function(X, Y, nvar){
  Y <- as.matrix(Y)
  N <- nrow(X)
  X <- X - rep(colMeans(X), each = N)
  Y <- Y - rep(colMeans(Y), each = N)
  if(missing(nvar)){
    nvar <- min(ncol(X), N-1)}
  T <- matrix(0, N, nvar)
  Q <- matrix(0, nvar, ncol(Y))
  P <- matrix(0, nvar, ncol(X))
  sel <- integer(nvar)
  for(i in 1:nvar){
    covs <- rowSums(crossprod(X, Y)^2)
    sel[i] <- which.max(covs)
    T[, i] <- X[, sel[i]]
    if(i > 1)
      T[, i] <- T[, i] - T[, 1:(i-1), drop=FALSE] %*% crossprod(T[, 1:(i-1), drop=FALSE], T[, i])
    T[, i] <- T[, i] / sqrt(sum(T[, i]^2))
    Q[i, ] <- crossprod(T[, i, drop=FALSE], Y)
    Y <- Y - T[, i, drop=FALSE] %*% Q[i, , drop=FALSE]
  }
  return(list(selected = sel, scores = T, loadings = crossprod(X, T), Yloadings = t(Q)))
}
