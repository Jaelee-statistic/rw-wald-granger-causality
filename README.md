# RW-Wald Granger-causality analysis: original research script

`original_analysis.R` is the primary script. It preserves the author's calculation order, functions, preprocessing, intermediate calculations, residual conventions, fixed-coefficient analysis, and bootstrap weight distributions and draw order. Repository changes cover relative file paths, English comments, and removal of the fixed RNG seed. Inactive plotting examples remain available.

## Related paper

This repository provides the empirical analysis code used in the author's paper:

**[Random-weighting bootstrap correction for Granger causality tests in VAR models with time-varying variance of unknown form (2026)](https://www.sciencedirect.com/science/article/abs/pii/S016517652600234X?via%3Dihub)**

The primary script, `original_analysis.R`, contains the empirical analysis code used in this study.

## Run

Open a fresh R session with this repository as the working directory. The original program clears the workspace. Install its original dependencies once:

```sh
Rscript --vanilla install_packages.R
Rscript --vanilla run.R
```

Alternatively, in R from the repository root:

```r
source("original_analysis.R", echo = TRUE, encoding = "UTF-8")
```

The wrapper echoes the original expressions so the final results are visible. It does not replace the analysis with functions or change its calculations. Package installation is separate from the analysis. All 19 original imports remain, including packages not directly used by the final statistics. Internet access and, depending on platform/package availability, build tools may be needed for installation.

## Preserved behavior

The analysis uses the current R random-number state without setting a seed. Each run draws new bootstrap weights, so bootstrap statistics and p-values can vary. To repeat the same draws, call `set.seed()` with your chosen seed before sourcing the script in R.

- `199` replications, sample size `175`, VAR order `1`, and original initial zero lag row.
- Original date conversion, sorting, numeric-column differencing, NA removal, merge, and positional construction of the model matrix.
- Original untransposed fitted-value multiplication in both OLS and bootstrap calculations. No transpose correction is applied.
- First-half and second-half residual means and covariance matrices, with their original loops and divisors.
- Full fixed coefficient vector `(0.43, 0, -1.02, -0.01)` and its separate statistics.
- Two-point random draw followed by exponential draw in every weight block, including the unused two-point weights.
- Bootstrap residual covariance and information-matrix calculations, even where not used by subsequent test statistics.
- Original matrix inversion and multiplication order, strict `>` exceedance comparisons, and division by `199` without a plus-one adjustment.

The final eight expressions display `Q_s_real`, `Q_s_real_pv`, `Q_WB_s_real_199`, `Q_WB_s_pvalue_real_199`, and their four fixed-coefficient counterparts. `Q_WB_s_real_199` is the last bootstrap statistic, not the entire distribution. Half-sample diagnostics remain available as `Omega_front` and `Omega_back` in an interactive session.

## Main results and supplementary reference results

The main empirical results are `Q_s_real`, `Q_s_real_pv`, and `Q_WB_s_pvalue_real_199`, based on the OLS estimator from the actual sample. These are the primary RW-Wald inference results.

Results with the `_fix` suffix are supplementary reference calculations using the prespecified coefficient vector `(0.43, 0, -1.02, -0.01)`. They are not the main RW-Wald inference results and do not represent a separately estimated model.

This clarification does not change the calculation methods, data, or original R code.

## Data

`data/bopbm_rawdata_2.txt` and `data/bopbsv_rawdata_2.txt` are byte-for-byte copies of the supplied attachments; only attachment filename suffixes were removed. No dates, values, delimiters, or line endings were rewritten. SHA-256 hashes are recorded in `data/SHA256SUMS`. Use these supplied snapshots for reproduction rather than downloading updated series.

The original script's private absolute paths have been replaced with `data/...`. Inactive plot examples use local `plot_1.png` through `plot_4.png` paths. Plotting remains commented out as in the original.

## Earlier condensed implementation

The previous repository contents are retained under `archive/condensed/` for history. They are **not the primary reproduction script**. That implementation refactors the original and includes a `transpose_corrected` branch; it omits original intermediate calculations. Its README, tests, and saved results refer only to that archived implementation. To inspect its old commands, first change the working directory to `archive/condensed/`.

## Verification

See `VERIFICATION.md` for the checks and their scope. To compare bootstrap outputs with the original script, use a common initial RNG state and the same R and package environment. Matching a published table has not been independently established. Bibliographic details and a software/data license are not invented here.
