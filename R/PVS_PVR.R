#' @title Principal Variable Selection (PVS)
#'
#' @description Greedy algorithm for extracting the most dominant (principal)
#' variables (X-columns) with respect to explained X-variance.
#'
#' @param X numeric predictor \code{matrix}.
#' @param nvar integer, the required number of selected variables.
#' @param ncomp integer, number of principal components included in the voting (default = all PCs).
#'
#' @return A list containing:
#' \item{Q}{Orthonormal scores (associated with the selected variables).}
#' \item{R}{Corresponding loadings. NOTE: R[,vperm] is upper triangular.}
#' \item{ids}{Indices arranged in the order of the nvar selected variables.}
#' \item{vperm}{Indices arranged in the order of the nvar selected and all non-selected variables. NOTE: R[,vperm] is upper triangular.}
#' \item{ssEX}{The variances explained by the selected variables.}
#' \item{ni}{The norms of the (residual) selected variables before the score-normalization (Q).}
#' \item{U}{The normalized PCA-scores.}
#' \item{s}{Singular values of the mean centered X.}
#'
#' @author Ulf Indahl, Kristian Hovde Liland.
#'
#' @references Joakim Skogholt, Kristian Hovde Liland, Tormod Næs, Age K. Smilde, Ulf Geir Indahl,
#' Selection of principal variables through a modified Gram–Schmidt process with and without supervision,
#' Journal of Chemometrics, Volume 37, Issue 10, Pages e3510 (2023), https://doi.org/10.1002/cem.3510
#'
#' @seealso \code{\link{PVR}}, \code{\link{VIP}}, \code{\link{filterPLSR}}, \code{\link{shaving}}, 
#' \code{\link{stpls}}, \code{\link{truncation}},
#' \code{\link{bve_pls}}, \code{\link{ga_pls}}, \code{\link{ipw_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{rep_pls}}, \code{\link{spa_pls}}.
#'
#' @examples
#' library(pls)
#' data(gasoline, package = "pls")
#' 
#' # PVS: Select 10 variables using all PCs in voting
#' pvs_result <- PVS(gasoline$NIR, nvar = 10)
#' 
#' # Compare with PCA using pcr() (octane is unused in PCA)
#' pca_result <- pcr(octane ~ NIR, ncomp = 10, data = gasoline, scale = FALSE)
#' 
#' # Plot cumulative variance explained
#' plot(cumsum(pvs_result$ssEX), type = "b", col = "blue", 
#'      xlab = "Number of Variables/Components", 
#'      ylab = "Cumulative % Variance Explained",
#'      main = "PVS vs PCA", ylim = c(0, 100))
#' pca_var <- 100 * cumsum(pca_result$Xvar) / pca_result$Xtotvar
#' lines(seq_along(pca_var), pca_var, type = "b", col = "red")
#' legend("bottomright", legend = c("PVS (10 variables)", "PCA (10 components)"),
#'        col = c("blue", "red"), lty = 1, pch = 1)
#'
#' @export
PVS <- function(X, nvar, ncomp = NULL) {
  # Initializations
  X <- as.matrix(X)
  mX <- colMeans(X)
  X <- sweep(X, 2, mX, "-")  # Centering of the data matrix
  
  svd_result <- svd(X)       # PCA of the datamatrix X
  U <- svd_result$u
  s <- svd_result$d
  
  nvar <- min(nvar, qr(X)$rank)
  eps <- 1e-10               # Tolerance for defining a vanishing norm
  n <- nrow(X)
  p <- ncol(X)
  
  ids <- rep(NA, nvar)          # Indices of the selected variables
  Q <- matrix(0, n, nvar)
  R <- matrix(0, nvar, p)
  ssEX <- numeric(nvar)         # Explained variances
  ni <- numeric(nvar)           # Norms of selected variables
  vAll <- numeric(p)         # Vector for holding the votes for all variables
  
  # PC's for voting - if not correctly specified
  if (is.null(ncomp)) {
    ncomp <- length(s)           # All PC's are included in the voting function
  } else {
    ncomp <- max(ncomp, 2)           # At least 2 PC's must be included in the voting function
  }
  
  T <- U[, seq_len(ncomp), drop = FALSE] * rep(s[seq_len(ncomp)], each = n)  # Non-normalized PCA-scores for implementing the voting
  ssTX <- sum(s^2)           # The total sum-of-squares
  
  a <- 1
  idn <- seq_len(p)          # Indices of candidate variables available for voting/selection
  
  while (a <= nvar) {
    Xidn <- X[, idn, drop = FALSE]
    vAll[idn] <- colSums((crossprod(T, Xidn))^2) / colSums(Xidn^2)  # Update the vAll-votes for the X-variables subject to selection
    si <- which.max(vAll)    # Identify the variable accounting for the maximum variance
    mx <- vAll[si]
    ni[a] <- sqrt(sum(X[, si]^2))  # Norm of selected candidate variable
    
    if (ni[a] > eps) {       # Consider only candidate variables with non-vanishing norms
      ids[a] <- si           # Index of the a-th selected variable
      qa <- X[, si] / ni[a]  # Unit vector in the direction of the selected variable
      R[a, idn] <- crossprod(qa, Xidn)  # Deflation coefficients in the chosen direction
      ssEX[a] <- 100 * mx / ssTX  # Store fraction of explained variance by the chosen variable
      X[, idn] <- Xidn - outer(qa, R[a, idn])  # Deflate X wrt the chosen variable (modified Gram-Schmidt step)
      Q[, a] <- qa           # Store qa into Q
      a <- a + 1             # Update a to prepare for the subsequent selection
    }
    
    idn <- setdiff(idn, si)  # Eliminate selected/vanishing variable from future voting assignments
    vAll[si] <- 0            # Set the future votes of the selected/vanishing variable to 0
  }
  
  vperm <- c(ids, setdiff(seq_len(p), ids))  # Indices of the selected and unselected variables
  
  return(list(Q = Q, R = R, ids = ids, vperm = vperm, ssEX = ssEX, ni = ni, U = U, s = s))
}

# MATLAB original code (commented out):
# function [Q, R, ids, vperm, ssEX, ni, U, s] = PVS(X, nvar, B)
# % Greedy algorithm for extracting the most dominant (principal)
# % variables (X-coulumns) w.r.t. explained X-variance.
# 
# % Inputs:
# % X - the datamatrix.
# % nvar - the required number of selected variables.
# % B - number of PC's included in the voting.
# 
# % Outputs:
# % Q   - orthonormal scores (associated with the selected variables).
# % R   - corresponding loadings. NOTE: R(:,vperm) is upper triangular.
# % ids - indices arranged in the order of the nvar selected variables.
# % vperm - indices arranged in the order of the nvar selected- and all non-selected variables. NOTE: R(:,vperm) is upper triangular.
# % ssEX - the variances explained by the selected variables.
# % ni  - the norms of the (residual) selected variables before the score-normalization (Q).
# % U   - the normalized PCA-scores.
# % s   - singular values of the mean centered X
# 
# % Initializations:
# mX  = mean(X); X = bsxfun(@minus,X,mX);  % Centering of the data matrix.
# [U, S, ~] = svd(X, 'econ'); s = diag(S); % Essentially PCA of the datamatrix X.
# nvar      = min(nvar,rank(X));
# eps    = 1e-10;                          % Tolerance for defining a vanishing norm.
# [n, p] = size(X);                        % Data dimensions.
# ids  = NaN(nvar,1);                         % Indices of the selected variables.
# Q    = zeros(n,nvar); R = zeros(nvar,p);       % See Output-description above.
# ssEX = zeros(nvar,1);                       % See Output-description above.
# ni   = zeros(nvar,1);                       % See Output-description above.
# vAll = zeros(1,p);                       % Vector for holding the votes for all variables.
# % PC's for voting - if not correctly specified:
# if nargin == 3
#     B = max(B,2);  % At least 2 PC's must be included in the voting function.
# elseif nargin < 3
#     B = length(s); % All PC's are included in the voting function.
# end
# T = U(:,1:B).*s(1:B)';                   % Non-normalized PCA-scores for implementing the voting.
# ssTX = sum(s.^2);                        % The total sum-of-squares.
# a = 1;  idn = 1:p;                       % Book-keeping: a-th variable to be selected & idn - indices of candidate variables available for voting/selection.
# while a <= nvar
#     Xidn = X(:,idn);
#     vAll(idn) = sum((T'*Xidn).^2)./sum(Xidn.^2); % Update the vAll-votes for the X-variables subject to selection.
#     [mx, si]  = max(vAll);                     % Identify & store (the index of) the variable accounting for the maximum variance.
#     ni(a)     = norm(X(:,si));                 % norm of selected candiate variable.
#     if ni(a) > eps                             % Consider only candidate variables with non-vanishing norms.
#         ids(a)  = si;                          % Index of the a-th selected variable.
#         qa      = X(:,si)./ni(a);              % Unit vector in the direction of the selected variable.
#         R(a,idn)= (qa'*Xidn);                  % Deflation coefficients in the chosen (v) direction.
#         ssEX(a) = 100*mx/ssTX;                 % Store fraction of explained variance by the chosen variable.
#         X(:,idn)= Xidn - qa.*R(a,idn);         % Deflate X wrt the chosen variable (modified Gram-Schmidt step).
#         Q(:,a)  = qa;                          % Store qa into Q.
#         a       = a+1;                         % Update a to prepare for the subsequent selection.
#     end
#     idn = setdiff(idn,si);                     % Eliminate selected/vanishing variable from future voting assignments.
#     vAll(si)=0;   %X(:,si) = 0;                % Set the future votes of the selected/vanishing variable to 0.                  
# end
# vperm = [ids; setdiff((1:p)', ids)];           % Indices of the selected- and unselected variables.

#' @title Principal Variable Regression (PVR)
#'
#' @description Greedy algorithm for extracting the most dominant variables/columns
#' with respect to simultaneous explained X-variance and squared correlation with Y.
#'
#' @param X numeric predictor \code{matrix}.
#' @param Y numeric response \code{vector} or \code{matrix} (single or multiple responses).
#' @param nvar integer, the required number of selected variables (default = 2).
#' @param ncomp integer, the number of principal components to include in the voting process (default = all PCs).
#'
#' @return A list containing:
#' \item{ids}{The indices of the selected variables.}
#' \item{betas}{The regression coefficients (including the constant term) for prediction of Y from the selected variables.}
#' \item{Q}{Orthonormal scores (associated with the selected variables).}
#' \item{R}{Corresponding loadings. NOTE: R[,vperm] is upper triangular.}
#' \item{vperm}{Indices arranged in the order of the nvar selected and all non-selected variables. NOTE: R[,vperm] is upper triangular.}
#' \item{U}{The normalized PCA-scores.}
#' \item{s}{Singular values of the mean centered X.}
#' \item{ssEX}{The X-variances explained by the selected variables.}
#' \item{ssEY}{The Y-variances explained by the selected variables.}
#' \item{ni}{The norms of the (residual) selected variables before the score-normalization (Q).}
#'
#' @author Ulf Indahl, Kristian Hovde Liland.
#'
#' @references Joakim Skogholt, Kristian Hovde Liland, Tormod Næs, Age K. Smilde, Ulf Geir Indahl,
#' Selection of principal variables through a modified Gram–Schmidt process with and without supervision,
#' Journal of Chemometrics, Volume 37, Issue 10, Pages e3510 (2023), https://doi.org/10.1002/cem.3510
#'
#' @seealso \code{\link{PVS}}, \code{\link{VIP}}, \code{\link{filterPLSR}}, \code{\link{shaving}}, 
#' \code{\link{stpls}}, \code{\link{truncation}},
#' \code{\link{bve_pls}}, \code{\link{ga_pls}}, \code{\link{ipw_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{rep_pls}}, \code{\link{spa_pls}}.
#'
#' @examples
#' library(pls)
#' data(gasoline, package = "pls")
#' 
#' # PVR: Select 10 variables using all PCs in voting
#' pvr_result <- PVR(gasoline$NIR, gasoline$octane, nvar = 10)
#' 
#' # Compare with PCR using all variables
#' pcr_result <- pcr(octane ~ NIR, ncomp = 10, data = gasoline, 
#'                   validation = "CV", scale = FALSE)
#' 
#' # Compare X-variance and Y-variance explained
#' par(mfrow = c(1, 2))
#' plot(cumsum(pvr_result$ssEX), type = "b", col = "blue", 
#'      xlab = "Number of Variables/Components", 
#'      ylab = "Cumulative % X-Variance",
#'      main = "X-Variance: PVR vs PCR",
#'      ylim = c(50, 100))
#' pcr_xvar <- 100 * cumsum(pcr_result$Xvar) / pcr_result$Xtotvar
#' lines(seq_along(pcr_xvar), pcr_xvar, type = "b", col = "red")
#' legend("bottomright", legend = c("PVR (10 vars)", "PCR (10 comps)"),
#'        col = c("blue", "red"), lty = 1, pch = 1)
#' 
#' plot(cumsum(pvr_result$ssEY), type = "b", col = "blue", 
#'      xlab = "Number of Variables/Components", 
#'      ylab = "Cumulative % Y-Variance",
#'      main = "Y-Variance: PVR vs PCR",
#'      ylim = c(0, 100))
#' pcr_yvar <- 100 * R2(pcr_result)$val[1,1,-1]
#' lines(seq_along(pcr_yvar), pcr_yvar, type = "b", col = "red")
#' legend("bottomright", legend = c("PVR (10 vars)", "PCR (10 comps)"),
#'        col = c("blue", "red"), lty = 1, pch = 1)
#' par(mfrow = c(1, 1))
#' 
#' # Predict using selected variables
#' X_selected <- gasoline$NIR[, pvr_result$ids]
#' y_pred_pvr <- cbind(1, X_selected) %*% pvr_result$betas[, ncol(pvr_result$betas)]
#' y_pred_pcr <- predict(pcr_result, ncomp = 10, newdata = gasoline)
#' 
#' # Compare RMSE (training error - same data used for fitting)
#' rmse_pvr <- sqrt(mean((gasoline$octane - y_pred_pvr)^2))
#' rmse_pcr <- sqrt(mean((gasoline$octane - y_pred_pcr)^2))
#' cat("RMSE - PVR:", round(rmse_pvr, 4), "\n")
#' cat("RMSE - PCR:", round(rmse_pcr, 4), "\n")
#'
#' @export
PVR <- function(X, Y, nvar = 2, ncomp = NULL) {
  # Initializations
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  mX <- colMeans(X)
  mY <- colMeans(Y)
  X <- sweep(X, 2, mX, "-")  # Centering of the data matrix
  Y <- sweep(Y, 2, mY, "-")  # Centering of the responses
  
  svd_result <- svd(X)       # PCA of X
  U <- svd_result$u
  s <- svd_result$d
  
  Y0 <- Y
  nvar <- min(nvar, length(s))
  eps <- 1e-10               # Tolerance
  n <- nrow(X)
  p <- ncol(X)
  
  ids <- rep(NA, nvar)          # Indices of the selected variables
  Q <- matrix(0, n, nvar)
  R <- matrix(0, nvar, p)
  ssEX <- numeric(nvar)         # X-variances explained
  ssEY <- numeric(nvar)         # Y-variances explained
  ni <- numeric(nvar)           # Norms of selected variables
  vAll <- numeric(p)         # Vector for holding the votes for all variables
  
  # PC's for voting - if not correctly specified
  if (is.null(ncomp)) {
    ncomp <- length(s)           # All PC's are included in the voting function
  } else {
    ncomp <- max(ncomp, 2)           # At least 2 PC's must be included in the voting function
  }
  
  T <- U[, seq_len(ncomp), drop = FALSE] * rep(s[seq_len(ncomp)], each = n)  # Non-normalized PCA-scores (including the variance information)
  ssTX <- sum(s^2)           # The total X sum-of-squares
  ssTy <- sum(Y^2)           # The total Y sum-of-squares
  
  a <- 1
  idn <- seq_len(p)          # Indices of candidate variables available for voting/selection
  
  while (a <= nvar) {
    Xidn <- X[, idn, drop = FALSE]
    sX2 <- colSums(Xidn^2)   # The remaining total X-sums-of-squares
    t <- colSums((crossprod(T, Xidn))^2) / sX2  # p-vector with the overall X-variances explained by each available X-variable
    
    # Update the vAll-votes for the X-variables subject to selection
    vAll[idn] <- (colSums((crossprod(Y, Xidn))^2) / sX2) * t
    
    si <- which.max(vAll)    # Identify the variable accounting for the maximal (modified) covariance
    ni[a] <- sqrt(sum(X[, si]^2))  # The norm of selected candidate variable
    
    if (ni[a] > eps) {       # Ignore candidate variables with too small norms
      ids[a] <- si           # Index of the a-th selected variable
      qa <- X[, si] / ni[a]  # Unit vector in the direction of the selected variable
      ssEX[a] <- 100 * sum((crossprod(qa, T))^2) / ssTX  # X-variance accounted for by the selected variable
      ssEY[a] <- 100 * sum((crossprod(qa, Y))^2) / ssTy  # Y-variance accounted for by the selected variable
      R[a, idn] <- crossprod(qa, Xidn)  # Deflation coefficients in the chosen direction
      X[, idn] <- Xidn - outer(qa, R[a, idn])  # Deflate X wrt the chosen variable (modified Gram-Schmidt step)
      Y <- Y - qa*c(crossprod(qa, Y))     # Deflate Y wrt the chosen variable
      Q[, a] <- qa           # Store qa into Q
      a <- a + 1             # Update a to prepare for the subsequent selection
    }
    
    idn <- setdiff(idn, si)  # Eliminate selected/vanishing variable from future voting assignments
    vAll[si] <- 0            # Set the future votes of the selected/vanishing variable to 0
  }
  
  vperm <- c(ids, setdiff(seq_len(p), ids))  # Indices of the selected and unselected variables
  
  # Calculate regression coefficients for using 1, 2, ..., nvar variables
  # Compute inv(R) and project Y onto Q scores
  Rinv <- solve(R[, ids, drop = FALSE])
  YtQ <- crossprod(Y0, Q)
  
  # Multiply each column of Rinv by corresponding element of YtQ, then cumsum across columns for each row
  betas <- t(apply(sweep(Rinv, 2, c(YtQ), "*"), 1, cumsum))
  
  # Add intercept term
  intercepts <- c(mY) - colSums(mX[ids] * betas)
  betas <- rbind(intercepts, betas)
  
  return(list(ids = ids, betas = betas, Q = Q, R = R, vperm = vperm, 
              U = U, s = s, ssEX = ssEX, ssEY = ssEY, ni = ni))
}

# MATLAB original code (commented out):
# function [ids, betas, Q, R, vperm, U, s, ssEX, ssEY, ni] = PVR(X, Y, A, B)
# % Principal Variable Regression
# % Last update: 30.11.2021 (Ulf Indahl)
# % Greedy algorithm for extracting the most dominant variables/columns
# % w.r.t. simultaneous explained X-variance and squared correlation with Y.
# 
# % Inputs:
# % X - the datamatrix
# % A - the required number of selected variables (indicated by "ids")
# % Y - the responses (single or multiple)
# % B - the number of PC's to include in the voting process
# 
# % Outputs:
# % ids  - the indices of the selected variables
# % betas- the regression coefficients (including the constant term) for prediction of y's from the selected variables.
# % Q    - orthonormal scores (associated with the selected variables)
# % R    - corresponding loadings. NOTE: R(:,vperm) is upper triangular.
# % vperm - indices arranged in the order of the A selected- and all non-selected variables. NOTE: R(:,vperm) is upper triangular.
# % U    - the normalized PCA-scores
# % s    - singular values of the mean centered X
# % ssEX - the X-variances explained by the selected variables
# % ssEY - the y-variances explained by the selected variables
# % ni   - the norms of the (residual) selected variables before the score-normalization (Q)
# 
# % Initializations:
# mX = mean(X); X = bsxfun(@minus,X,mX);   % Centering of the data matrix.
# mY = mean(Y); Y = bsxfun(@minus,Y,mY);   % Centering of the responses.
# [U, S, ~] = svd(X, 'econ'); s = diag(S); % PCA of X.
# Y0 = Y;
# if nargin == 2, A = 2; else A = min(A, length(s)); end
# eps  = 1e-10;                            % Tolerance.
# [n, p] = size(X);                        % Data dimensions
# ids  = NaN(A,1);                         % Indices of the selected variables & the explained (X) sum-of-squares accounted for by the selected variables.
# Q    = zeros(n,A); R = zeros(A,p);       % See Output-description above.
# ssEX = zeros(A,1); ssEY = zeros(A,1);    % See Output-description above.
# ni   = zeros(A,1);                       % See Output-description above.
# vAll = zeros(1,p);                       % Vector for holding the votes for all variables.
# if nargin == 4% PC's for voting - if not correctly specified:
#     B = max(B,2);  % At least 2 PC's must be included in the voting function.
# elseif nargin < 4
#     B = length(s); % All PC's are included in the voting function.
# end
# 
# T = U(:,1:B).*s(1:B)';                   % Non-normalized PCA-scores (including the variance information of the Principal components).
# ssTX = sum(s.^2);                        % The total X sum-of-squares.
# ssTy = sum(Y(:).^2);                     % The total Y sum-of-squares.
# 
# a = 1;  idn = 1:p;                       % Book-keeping: a-th variable to be selected & idn - indices of candidate variables available for voting/selection.
# while a <= A
#     Xidn = X(:,idn);
#     sX2  = sum(Xidn.^2);                 % The remaining total X-sums-of-squares.
#     t    = sum((T'*Xidn).^2)./sX2;       % p-vector with the overall X-variances explained by each of the available X-variables.
#     %vAll(idn) = prod(((Y'*Xidn).^2)./sX2,1).*t; % Update the vAll-votes for the X-variables subject to selection.
#     vAll(idn) = (((Y'*Xidn).^2)./sX2).*t; % Update the vAll-votes for the X-variables subject to selection.
# 
#     [~, si]   = max(vAll);               % Identify & store the variable accounting for the maximal (modified) covariance.
#     ni(a)     = norm(X(:,si));           % the norm of selected candiate variable.
#     if ni(a) > eps                       % Ignore for candidate variables with too small norms.
#         ids(a)  = si;                    % Index of the a-th selected variable.
#         qa      = X(:,si)./ni(a);        % Unit vector in the direction of the selected variable.
#         ssEX(a) = 100*sum((qa'*T).^2)/ssTX; % X-variance accouted for by the selected variable.
#         ssEY(a) = 100*sum((qa'*Y).^2)/ssTy; % y-variance accouted for by the selected variable.
#         R(a,idn)= (qa'*Xidn);            % Deflation coeficients in the chosen (v) direction.
#         X(:,idn)= Xidn - qa.*R(a,idn);   % Deflate X wrt the chosen variable (modified Gram-Schmidt step).
#         Y       = Y - qa.*(qa'*Y);       % Deflate Y wrt the chosen variable.
#         Q(:,a)  = qa;                    % Store v into Q.
#         a       = a+1;                   % Update k to prepare for the subsequent selection.
#     end
#     idn = setdiff(idn,si);               % Eliminate selected/vanishing variable from future voting assignments.
#     vAll(si)=0;                          % Set the future votes of the selected/vanishing variable to 0.
# end
# 
# vperm = [ids; setdiff((1:p)', ids)];     % Indices of the selected- and unselected variables.
# % betas = cumsum(bsxfun(@times,eye(A)/R(:,ids), Y0'*Q),2); % betas = cumsum(bsxfun(@times,inv(R(:,ids)), Y0'*Q),2); %
# betas = cumsum((eye(A)/R(:,ids)).*(Y0'*Q),2);
# betas = [mY - mX(ids)*betas; betas];     % The regression coeffs associated with the selected variables.
