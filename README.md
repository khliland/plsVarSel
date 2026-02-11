# Variable selection methods for Partial Least Squares - plsVarSel

## Installation

``` r
# Install release version from CRAN  
install.packages("plsVarSel")  
# Install development version from GitHub  
devtools::install_github("khliland/plsVarSel")
```

## Contents

- Filter methods
    - VIP - Variable Importance in Projections
    - SR - Selectivity Ratio
    - sMC - Significance Multivariate Correlation
    - LW - Loading Weights
    - RC - Regression Coefficients
    - URC - RC scaled as abs(RC)/max(abs(RC))
    - FRC - URC further scaled as URC/PRESS
    - mRMR - Minimum Redundancy Maximal Relevancy
- Wrapper methods
    - BVE-PLS -	Backward variable elimination PLS
    - GA-PLS - Genetic algorithm combined with PLS regression
    - IPW-PLS - Iterative predictor weighting PLS
    - MCUVE-PLS - Uninformative variable elimination in PLS
    - REP-PLS - Regularized elimination procedure in PLS
    - SPA-PLS - Sub-window permutation analysis coupled with PLS
    - T2-PLS - Hotelling's T^2 based variable selection in PLS
    - WVC-PLS - Weighted Variable Contribution in PLS
- Embedded methods
    - Trunction PLS
    - ST-PLS - Soft-Threshold PLS
    - CovSel - Covariance Selection
    - PVS/PVR - Principal Variable Selection and Regression
- LDA wrappers for PLS classficiations and cross-validation
- Shaving - Repeated shaving of variables using filters (experimental)
- Simulation tools

## Main references (more in package)
- T. Mehmood, K.H. Liland, L. Snipen, S. Sæbø, A review of variable selection methods in Partial Least Squares Regression, Chemometrics and Intelligent Laboratory Systems 118 (2012) 62-69. 
- T. Mehmood, S. Sæbø, K.H. Liland, Comparison of variable selection methods in partial least squares regression, Journal of Chemometrics 34 (2020) e3226.
