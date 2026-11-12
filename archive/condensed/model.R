read_series <- function(name) {
  x <- read.table(file.path("data", paste0(tolower(name), ".txt")), header = TRUE)
  stopifnot(identical(names(x), c("date", name)), all(is.finite(x[[2]])))
  x$date <- as.Date(x$date, "%Y-%m-%d")
  x <- x[order(x$date), ]
  stopifnot(!anyNA(x$date), !anyDuplicated(x$date),
            identical(x$date, seq(min(x$date), max(x$date), by = "quarter")))
  x
}

make_design <- function() {
  m <- read_series("BOPBM")
  s <- read_series("BOPBSV")
  stopifnot(identical(m$date, s$date))
  y <- cbind(BOPBM = diff(m$BOPBM), BOPBSV = diff(s$BOPBSV))
  list(y = y, z = rbind(c(0, 0), head(y, -1)), date = m$date[-1])
}

fit_var <- function(y, z, weights) {
  # Rows of A index responses. Columns index lagged predictors.
  t(solve(crossprod(z, weights * z), crossprod(z, weights * y)))
}

run_analysis <- function(design, B = 199L) {
  y <- design$y
  z <- design$z
  n <- nrow(y)
  A <- fit_var(y, z, rep(1, n))
  # For vec(A) stacked by column, R = (0,0,1,0) selects A[1,2].
  # R J^-1 R' reduces to inv(Z'Z/n)[2,2] * Omega[1,1].
  lag_factor <- solve(crossprod(z) / n)[2, 2]
  modes <- c("legacy", "transpose_corrected")
  denom <- sapply(modes, function(mode) {
    fitted <- if (mode == "legacy") z %*% A else z %*% t(A)
    lag_factor * mean((y - fitted)[, 1]^2)
  })
  cases <- expand.grid(fixed = c(FALSE, TRUE), mode = modes, stringsAsFactors = FALSE)
  # The fixed vector is a source-script diagnostic, not a second fitted model.
  center <- ifelse(cases$fixed, -1.02, A[1, 2])
  scale <- n / denom[cases$mode]
  observed <- scale * center^2
  draws <- matrix(NA_real_, B, nrow(cases))
  for (k in seq_len(B)) {
    w <- numeric(n)
    for (i in seq_len(n)) {
      # Keep the unused two-point draw to reproduce the original RNG stream.
      sample(c((3-sqrt(5))/2, (3+sqrt(5))/2), 1,
             prob = c(1+sqrt(5), sqrt(5)-1)/(2*sqrt(5)))
      w[i] <- rexp(1)
    }
    star <- fit_var(y, z, w)
    draws[k, ] <- scale * (star[1, 2] - center)^2
  }
  exceed <- sweep(draws, 2, observed, ">")
  summary <- data.frame(cases, wald = observed, chi_p = pchisq(observed, 1, lower.tail = FALSE),
                        last = draws[B, ], exceed = colSums(exceed), rw_p = colMeans(exceed))
  colnames(draws) <- paste(cases$mode, ifelse(cases$fixed, "fixed", "estimated"), sep = "_")
  list(summary = summary, draws = draws, exceed = exceed, A = A)
}
