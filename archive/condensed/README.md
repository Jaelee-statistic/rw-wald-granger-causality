# Random-weight Wald test for Granger causality

This repository reproduces the supplied empirical R script for a bivariate VAR(1) in first differences. The question is whether lagged changes in the balance on services (BOPBSV) help predict changes in the merchandise balance (BOPBM), conditional on its own lag.

## Run from the repository root

Use base R (tested on Windows with R 4.6.0). No additional packages or downloads are required.

```sh
Rscript --vanilla run.R
Rscript --vanilla tests/verify.R
```

The first command prints the four analysis rows and writes results/summary.csv, results/coefficients.csv (response rows, predictor columns), and results/bootstrap.csv (one row per draw). The second checks the original outputs and indicators, dates, lag alignment, coefficient orientation, covariance and bootstrap tail calculation. The included result files were generated with these settings.

## Data and timing

The two tab-delimited files in data/ are the supplied raw files with shorter names. Each has 176 quarterly levels from 1970 Q1 to 2013 Q4. Input dates are already ISO dates. The loader sorts and validates the quarterly grids and requires the dates to match before differencing. The original constructs a merged table but forms its model matrix by position. Columns are ordered BOPBM, BOPBSV. The 175 responses start in 1970 Q2. Following the source, the first lag row is zero and subsequent rows contain the preceding difference. The first response therefore contributes to residual covariance but carries zero weight in the coefficient normal equations. There is no intercept or trend.

## Statistical audit

* The equation is y_t = A z_t + u_t with column vectors. The data matrices Y and Z are 175 by 2. The original calculates A in response-by-predictor orientation but predicts with Z A. The legacy row reproduces this operation. The transpose_corrected row predicts with Z A' so the fitted residuals satisfy the OLS normal equations. Both rows use the same A and bootstrap estimates.
* Column-major vec(A) is (A11,A21,A12,A22). The one-degree-of-freedom null is A12 = 0, services to merchandise. The observed test uses unrestricted OLS. The bootstrap approximates the null fluctuation by centering A*12 at the estimated A12. It holds the lagged design and responses fixed and reweights the estimating equations. It does not simulate a restricted VAR.
* With G = Z'Z/n and Omega = U'U/n, the source uses J = G tensor inverse(Omega). The contrast variance factor d is inverse(G)[2,2] times Omega[1,1]. Both observed and bootstrap statistics use the same observed d. The source computes a bootstrap-specific J but never uses it. Half-sample covariances likewise do not enter its reported tests. The public implementation omits those unused calculations. The chi-square column is the source's conventional reference, not a separately validated heteroskedasticity-robust calibration.
* Each of 199 draws uses independent unit-mean exponential weights (block length one), seed 165. The unused two-point draws are retained in their original sequence. The source uses the strict upper tail, count(Q* > Q)/199, with no plus-one adjustment. Because the common positive scale cancels in each comparison, correcting only the residual transpose changes both Wald scales but preserves every bootstrap exceedance mathematically.
* Rows with fixed=TRUE reproduce a separate source diagnostic that substitutes (0.43,0,-1.02,-0.01) for vec(A) in the observed numerator and bootstrap center. Its origin is unspecified. Use fixed=FALSE for the estimated-coefficient test. The last column named last in the summary is the final bootstrap statistic, a diagnostic only.

## Verification scope

The reference CSVs under tests/ were generated from the supplied script's original loops. The reference run redirected the two input paths, replaced the sorting/differencing pipeline with equivalent base R operations for these complete data, removed package imports and skipped the unused external wald.test diagnostic. The reference preserved the empirical statistic formulas and random draws. Both branches match the original Wald values, chi-square tails, final draw statistics, exceedance counts and empirical p-values within numeric tolerance, with exact matching indicators for all draws in legacy mode. The corrected observed Wald is separately checked against equation-wise least squares and the full Kronecker covariance. This comparison establishes reproduction of the supplied script. Agreement with a published table and data redistribution terms require author confirmation.
