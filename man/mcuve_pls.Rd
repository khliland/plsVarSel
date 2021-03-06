% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UVE.R
\name{mcuve_pls}
\alias{mcuve_pls}
\title{Uninformative variable elimination in PLS (UVE-PLS)}
\usage{
mcuve_pls(y, X, ncomp = 10, N = 3, ratio = 0.75, MCUVE.threshold = NA)
}
\arguments{
\item{y}{vector of response values (\code{numeric} or \code{factor}).}

\item{X}{numeric predictor \code{matrix}.}

\item{ncomp}{integer number of components (default = 10).}

\item{N}{number of samples Mone Carlo simulations (default = 3).}

\item{ratio}{the proportion of the samples to use for calibration (default = 0.75).}

\item{MCUVE.threshold}{thresholding separate signal from noise (default = NA creates 
automatic threshold from data).}
}
\value{
Returns a vector of variable numbers corresponding to the model 
having lowest prediction error.
}
\description{
Artificial noise variables are added to the predictor set before the PLSR 
model is fitted. All the original variables having lower "importance" than the artificial 
noise variables are eliminated before the procedure is repeated until a stop criterion is 
reached.
}
\examples{
data(gasoline, package = "pls")
with( gasoline, mcuve_pls(octane, NIR) )

}
\references{
V. Centner, D. Massart, O. de Noord, S. de Jong, B. Vandeginste, C. Sterna, 
Elimination of uninformative variables for multivariate calibration, Analytical Chemistry 
68 (1996) 3851-3858.
}
\seealso{
\code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{shaving}}, 
\code{\link{stpls}}, \code{\link{truncation}},
\code{\link{bve_pls}}, \code{\link{ga_pls}}, \code{\link{ipw_pls}}, \code{\link{mcuve_pls}},
\code{\link{rep_pls}}, \code{\link{spa_pls}},
\code{\link{lda_from_pls}}, \code{\link{lda_from_pls_cv}}, \code{\link{setDA}}.
}
\author{
Tahir Mehmood, Kristian Hovde Liland, Solve Sæbø.
}
