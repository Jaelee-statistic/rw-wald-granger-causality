source("model.R")
d <- make_design()
stopifnot(nrow(d$y) == 175, all(d$z[1, ] == 0),
          all(d$z[-1, ] == head(d$y, -1)),
          identical(as.character(range(d$date)), c("1970-04-01", "2013-10-01")))
r <- run_analysis(d)
ref <- read.csv("tests/original-reference.csv")
stopifnot(isTRUE(all.equal(r$summary[1:2, names(ref)], ref, check.attributes = FALSE)),
          all(r$exceed[, 1:2] == as.matrix(read.csv("tests/original-indicators.csv"))))
# Independently estimate both equations and test the services coefficient.
ols <- lm.fit(d$z, d$y)
stopifnot(max(abs(t(ols$coefficients) - r$A)) < 1e-10,
          max(abs(crossprod(d$z, d$y - d$z %*% t(r$A)))) < 1e-8)
covariance <- kronecker(solve(crossprod(d$z)), crossprod(ols$residuals) / nrow(d$y))
R <- matrix(c(0, 0, 1, 0), 1)
wald <- drop((R %*% as.vector(r$A))^2 / (R %*% covariance %*% t(R)))
stopifnot(abs(wald - r$summary$wald[3]) < 1e-10,
          all(r$summary$rw_p == colMeans(sweep(r$draws, 2, r$summary$wald, ">"))))
cat("PASS: time alignment, both original branches, all 398 indicators, OLS and Wald direction, bootstrap tails.\n")
