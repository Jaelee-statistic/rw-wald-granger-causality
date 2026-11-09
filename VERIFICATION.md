# Preservation verification — 2026-11-09

## Source and changes

The baseline is the attached original R research script supplied by the author, not the earlier condensed implementation. A local verification copy redirects only its private paths. This private baseline is not published because the supplied original contains personal directory names.

The primary public script translates the original comments into English, adds English explanatory comments, and redirects the two data paths. Private paths in inactive plotting comments are also sanitized. No executable statements are removed, reordered, refactored, or corrected. The original package imports and inactive code remain. All original Korean comments have been replaced with English explanations. Comment-only corrections clarify array dimensions and the stored exceedance indicators without changing executable statements.

## Completed checks

Using R 4.6.0 on Windows, `parse(..., keep.source = FALSE)` for the entire path-normalized original and primary script returned exactly identical expression trees. This covers all executable formulas, control flow, constants, imports, preprocessing, covariance calculations, diagnostics, and output expressions, rather than a selection of final formulas.

The original and revised weight-generation loops were separately executed for 199 iterations with 175 blocks and seed 165. All 34,825 two-point weights, all 34,825 exponential weights, and all 199 iteration-end `.Random.seed` states matched exactly using `identical()`. RNG kinds were `Mersenne-Twister`, `Inversion`, and `Rejection`. This is a scoped RNG check, not a claim that the complete analysis was run.

Both published data files were copied as bytes from the supplied attachments and checked against them. Their SHA-256 checksums are in `data/SHA256SUMS`. A scan of the proposed public files found no remaining personal absolute filesystem paths.

To repeat the structural and scoped RNG checks, prepare a local copy of the author's original with only its file paths redirected, then run from the repository root:

```sh
Rscript --vanilla tests/verify_preservation.R /path/to/path_only_original.R
```

## Not yet verified

The installed R environment lacks 18 of the 19 imported packages (only MASS is present). Consequently, neither full original nor full revised analysis was executed in this verification. Numerical equality of final Wald statistics, p-values, intermediate covariance matrices, and the complete-run RNG state remains unverified. Dependencies were not removed or replaced to manufacture a successful run. The startup also reported Windows locale-setting warnings.

For full verification, install the dependencies and run each complete script in a separate fresh R process with identical R/package versions and RNG settings. Compare all original result objects, both covariance diagnostics, both 199-element exceedance vectors, and the final `.Random.seed` exactly. Record `sessionInfo()` for both runs. The old archived CSVs are not new validation evidence for this primary script.

No comparison against a published paper table has been performed. GitHub publication is a separate action requiring the author's approval.

## English-only documentation update

After translating the comments, the complete parsed-program equality and scoped RNG checks were rerun and passed. All tracked public text files were scanned for Hangul; no Korean text remains in the current files. The full-analysis verification limitations above remain unchanged.
