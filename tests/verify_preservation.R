# Compare against the supplied original after ONLY its private paths are redirected.
# Usage from repository root: Rscript --vanilla tests/verify_preservation.R PATH_ONLY_ORIGINAL
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
original <- parse(args[1], keep.source = FALSE, encoding = "UTF-8")
revised <- parse("original_analysis.R", keep.source = FALSE, encoding = "UTF-8")
stopifnot(identical(original, revised))
cat("PASS: complete parsed programs are identical after path normalization.\n")

# Execute only the unchanged weight-generation loop, without substituting
# dependencies or preprocessing in the full analysis. This is a scoped RNG check.
find_weight_loop <- function(expr) {
  if (!is.call(expr)) return(NULL)
  if (length(expr) < 2L) return(NULL)
  if (identical(expr[[1]], as.name("for")) &&
      identical(expr[[2]], as.name("i"))) return(expr)
  for (j in seq.int(2L, length(expr))) {
    if (identical(expr[[j]], quote(expr = ))) next
    found <- find_weight_loop(expr[[j]])
    if (!is.null(found)) return(found)
  }
  NULL
}
trace_weights <- function(program) {
  loop <- NULL
  for (expr in program) {
    loop <- find_weight_loop(expr)
    if (!is.null(loop)) break
  }
  stopifnot(!is.null(loop))
  set.seed(165)
  e <- new.env()
  e$block_size <- 1
  e$block_number <- 175
  states <- vector("list", 199)
  for (k in 1:199) {
    e$block_w_be <- matrix(0, 175, 1)
    e$block_w_ex <- matrix(0, 175, 1)
    eval(loop, envir = e)
    states[[k]] <- list(be = e$block_w_be, ex = e$block_w_ex,
                       seed = .Random.seed)
  }
  states
}
stopifnot(identical(trace_weights(original), trace_weights(revised)))
cat("PASS: 34,825 two-point weights, 34,825 exponential weights, and all 199 iteration-end RNG states match exactly.\n")
cat("This scoped check does not execute the full package-dependent analysis.\n")
print(R.version.string)
print(RNGkind())
