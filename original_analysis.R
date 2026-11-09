# PRIMARY REPRODUCTION SCRIPT
# Only file paths and explanatory comments differ from the supplied original.
# Run from the repository root; see README.md for packages and execution.
# Original comments are translated into English; inactive plotting code is retained.

#=================================================================================================|
#
# Wild bootstrap based Wald test for VAR model with time-varying variance of unknown form
#
#=================================================================================================|

# Start in a fresh R session. This original statement clears the current workspace.
rm(list=ls()) # reset

# [library] ======================================================================================|
# Original package imports are retained in their original order, including unused imports.
library(vars)
library(tsDyn)
library(MTS)
library(TSA)
library(mvtnorm)
library(matrixcalc)
library(fastmatrix)
library(MASS)
library(dplyr)
library(portes)
library(urca)
library(rugarch)
library(aod)
library(readxl)
library(ggplot2)
library(fBasics)
library(fUnitRoots)
library(timeSeries)
library(svars)
#=================================================================================================|



#------------------------------------------------------------#
# 1. Load the empirical data
#------------------------------------------------------------#
# Original RNG seed. Do not move this statement or insert random draws below it.
set.seed(165)

# Read the raw data from a text file.
# Read the supplied merchandise-balance file without changing its contents.
bopbm_rawdata <- read.table("data/bopbm_rawdata_2.txt",header=T)
# "data/
# Preserve the original date conversion, including the appended day string.
bopbm_rawdata$date <- as.Date(paste(bopbm_rawdata$date, "-01", sep = ""), format = "%Y-%m-%d")

# # Create a time-series plot.
# png("plot_1.png", width=6000, height=4000, res=500, type="cairo")
# ggplot(bopbm_rawdata, aes(x = date, y = BOPBM)) +
#   geom_line() + # Draw the time-series line.
#   labs(x = "Date", y = "BOPBM") + # Set the x-axis and y-axis labels.
#   scale_x_date(date_breaks = "5 year", date_labels = "%Y-%m") + # Set the date format on the x-axis.
#   theme_test() # Apply the plotting theme.

# Calculate first differences.
# Original preprocessing: sort by date, difference numeric columns, then omit NA rows.
bopbm_diff <- bopbm_rawdata %>%
  arrange(date) %>%
  mutate(across(where(is.numeric), ~ c(NA, diff(.)))) %>%
  na.omit()  # Remove NA rows introduced by differencing.


# # Create a time-series plot.
# png("plot_2.png", width=6000, height=4000, res=500, type="cairo")
# ggplot(bopbm_diff, aes(x = date, y = BOPBM)) +
#   geom_line() + # Draw the time-series line.
#   labs(x = "Date", y = "BOPBM") + # Set the x-axis and y-axis labels.
#   scale_x_date(date_breaks =  "5 year", date_labels = "%Y-%m") + # Set the date format on the x-axis.
#   theme_test() # Apply the plotting theme.


# Read the raw data from a text file.
# Read the supplied services-balance file without changing its contents.
bopbsv_rawdata <- read.table("data/bopbsv_rawdata_2.txt",header=T)
# Preserve the original date conversion exactly as supplied.
bopbsv_rawdata$date <- as.Date(paste(bopbsv_rawdata$date, "-01", sep = ""), format = "%Y-%m-%d")

# # Create a time-series plot.
# png("plot_3.png", width=6000, height=4000, res=500, type="cairo")
# ggplot(bopbsv_rawdata, aes(x = date, y = BOPBSV)) +
#   geom_line() + # Draw the time-series line.
#   labs(x = "Date", y = "BOPBSV") + # Set the x-axis and y-axis labels.
#   scale_x_date(date_breaks ="5 year", date_labels = "%Y-%m") + # Set the date format on the x-axis.
#   theme_test() # Apply the plotting theme.


# Calculate first differences.
# Apply the original services-series preprocessing without rewriting the pipeline.
bopbsv_diff <- bopbsv_rawdata %>%
  arrange(date) %>%
  mutate(across(where(is.numeric), ~ c(NA, diff(.)))) %>%
  na.omit()  # Remove NA rows introduced by differencing.

# # Create a time-series plot.
# png("plot_4.png", width=6000, height=4000, res=500, type="cairo")
# ggplot(bopbsv_diff, aes(x = date, y = BOPBSV)) +
#   geom_line() + # Draw the time-series line.
#   labs(x = "Date", y = "BOPBSV") + # Set the x-axis and y-axis labels.
#   scale_x_date(date_breaks ="5 year", date_labels = "%Y-%m") + # Set the date format on the x-axis.
#   theme_test() # Apply the plotting theme.


# bopbm_diff: first differences of the merchandise balance.
# bopbsv_diff: first differences of the services balance.

# Construct the empirical dataset.
# real_data <- merge(bopbsv_diff, bopbm_diff, by = "date")
# a <- data.frame(bopbsv_diff$BOPBSV,bopbm_diff$BOPBM)
# Retain both the merged table and the positional construction of the model matrix below.
real_data <- merge(bopbm_diff,bopbsv_diff, by = "date")
a <- data.frame(bopbm_diff$BOPBM,bopbsv_diff$BOPBSV)
var_real <- as.matrix(a)



#----------------------------------------------------------------------------------------------#
#  1. Parameter definition
#----------------------------------------------------------------------------------------------#


mean_resid    <- c(0, 0)
# Original fixed sample size: 175 differenced observations from 176 input rows.
number        <- 175                                       # 175 differenced observations from 176 raw observations.
block_size    <- floor((number)^(0))                       # block size
p             <- 1                                         # Model order: VAR(p).
df            <- 1                                         # Degrees of freedom: number of restrictions (1 or 2).
cv            <- qchisq(0.95, df)                          # critical value
rep_WB_199    <- 199                                       # B=199
Q_wald_WB_real_199      <- matrix(0, nrow = rep_WB_199, ncol = 1)      # Store the exceedance indicator for each of the B bootstrap replications.
Q_wald_WB_real_199_fix  <- matrix(0, nrow = rep_WB_199, ncol = 1)      # Store the fixed-coefficient exceedance indicator for each replication.
# Column-major coefficient order is (a11, a21, a12, a22); the restriction selects a12.
R             <- matrix(c(0,0,1,0), nrow=1)                           # Restriction a12 = 0; column-major order: a11, a21, a12, a22.
half_number   <- floor(number/2)                    # Index separating the first and second halves.


#------------------------------------------------------------#
# 2. Calculate OLS estimates, residuals, and the residual covariance matrix
#------------------------------------------------------------#

# OLS estimate
# Construct the original response and lag arrays, retaining the initial zero lag row.
X_real <- var_real[1:number, ]                            # Construct the number-by-2 response matrix X_t.
X_tilde_real <- var_real[1:(number - 1), ]
X_tilde_real <- rbind(c(0,0),X_tilde_real)                # Construct X_(t-1): prepend a zero row and omit the last observation.

Sigma_Xt_real        <- matrix(0, nrow = 2, ncol = 2)     # Initialize a 2-by-2 matrix.
Sigma_Xt_tilde_real  <- matrix(0, nrow = 2, ncol = 2)     # Initialize a 2-by-2 matrix.


# Accumulate X_t %*% t(X_tilde) over all observations.
for (t in 1:number) {
  xt_real              <- matrix(X_real[t, ], ncol = 1)                             # 2x1
  xt_tilde_real        <- matrix(X_tilde_real[t, ], ncol = 1)                       # 2x1
  Sigma_Xt_real        <- Sigma_Xt_real + (xt_real %*% t(xt_tilde_real))              # 2x2
  Sigma_Xt_tilde_real  <- Sigma_Xt_tilde_real + (xt_tilde_real %*% t(xt_tilde_real))  # 2x2
}
#real
Sigma_Xt_real       <- Sigma_Xt_real / number         #sigma_xt
Sigma_Xt_tilde_real <- Sigma_Xt_tilde_real / number   #sigma_xt_tilde

# Calculate the OLS coefficient estimate.
# OLS estimate using the original accumulated matrix products and inversion order.
OLS_theta_real <- Sigma_Xt_real %*% solve(Sigma_Xt_tilde_real)  # 2x2 OLS estimate

# Calculate OLS residuals.
# Reproduction convention: retain the untransposed coefficient matrix exactly as supplied.
X_hat_real <- X_tilde_real %*% (OLS_theta_real)   # Calculate fitted values using the original coefficient orientation.

resi_OLS_real <- X_real - X_hat_real              # Construct the number-by-2 OLS residual matrix.

# Calculate the 2-by-2 contemporaneous residual covariance using u_t %*% t(u_t).
Omega_OLS_real <- list()

S_real <- matrix(0, nrow = ncol(resi_OLS_real), ncol = ncol(resi_OLS_real))  # Initialize a 2-by-2 matrix.

for (t in 1:number) {
  ut_real     <- matrix(resi_OLS_real[t, ], ncol = 1)            # Construct the 2-by-1 residual vector u_t.
  S_real      <- S_real + (ut_real %*% t(ut_real))               # Accumulate the 2-by-2 residual outer products.
}
# Residual variance-covariance matrix.
Omega_OLS_real <- S_real / number


#------------------------------------------------------------#
# 5. Split residuals into the first and second halves
#------------------------------------------------------------#

# Retain the separate first-half and second-half residual covariance diagnostics.
resi_OLS_real_front <- resi_OLS_real[1:half_number, ]
resi_OLS_real_back  <- resi_OLS_real[(half_number+1):number, ]

number_front <- nrow(resi_OLS_real_front)
number_back  <- nrow(resi_OLS_real_back)


#------------------------------------------------------------#
# 6. Calculate the mean residual vector for the first half
#------------------------------------------------------------#
# Calculate the first-half residual mean using the original summation loop.
mean_front <- matrix(0, nrow = 1, ncol = 2)

for (t in 1:number_front) {
  mean_front <- mean_front + matrix(resi_OLS_real_front[t, ], nrow = 1)
}
mean_front <- mean_front / number_front


#------------------------------------------------------------#
# 7. Calculate the mean residual vector for the second half
#------------------------------------------------------------#
# Calculate the second-half residual mean using the original summation loop.
mean_back <- matrix(0, nrow = 1, ncol = 2)

for (t in 1:number_back) {
  mean_back <- mean_back + matrix(resi_OLS_real_back[t, ], nrow = 1)
}
mean_back <- mean_back / number_back


#------------------------------------------------------------#
# 8. Calculate the residual covariance matrix for the first half
#------------------------------------------------------------#
# Center first-half residuals by their own mean; retain the original divisor.
S_front <- matrix(0, nrow = 2, ncol = 2)

for (t in 1:number_front) {
  ut_front <- matrix(resi_OLS_real_front[t, ] - mean_front, ncol = 1)
  S_front  <- S_front + (ut_front %*% t(ut_front))
}

Omega_front <- S_front / number_front


#------------------------------------------------------------#
# 9. Calculate the residual covariance matrix for the second half
#------------------------------------------------------------#
# Center second-half residuals by their own mean; retain the original divisor.
S_back <- matrix(0, nrow = 2, ncol = 2)

for (t in 1:number_back) {
  ut_back <- matrix(resi_OLS_real_back[t, ] - mean_back, ncol = 1)
  S_back  <- S_back + (ut_back %*% t(ut_back))
}

Omega_back <- S_back / number_back




#------------------------------------------------------------#
# 2-1. Calculate the Wald statistic for the OLS estimate
#------------------------------------------------------------#

OLS_theta_real_vec <- as.vector((OLS_theta_real))
X_tiled_real_J <- matrix(0, nrow = 2, ncol = 2)
# Retain the complete fixed coefficient vector and its separate analysis.
OLS_theta_real_vec_fix <- as.vector((c(0.43,0,-1.02,-0.01)))

# Calculate the information matrix J and the Wald statistic.

for (t in 1:number) {
  X_til_real_J    <- matrix(X_tilde_real[t, ], ncol = 1)                  # Construct the 2-by-1 lagged observation vector X_(t-1).
  X_tiled_real_J  <- X_tiled_real_J + (X_til_real_J %*% t(X_til_real_J))    # Accumulate the 2-by-2 lagged cross-product matrix.
}

#real
# Original Kronecker information matrix and Wald formulas; no scalar simplification.
J_hat_real   <- (number^(-1)*X_tiled_real_J) %x% solve(Omega_OLS_real)
Q_s_real     <- number * t(OLS_theta_real_vec)%*%t(R) %*% solve(R %*% solve(J_hat_real) %*% t(R)) %*% R %*% OLS_theta_real_vec

# Wald test
Q_s_real_pv  <- 1 - pchisq(Q_s_real, df=nrow(R))

# Calculate the statistics using the fixed coefficient vector.
Q_s_real_fix <- number * t(OLS_theta_real_vec_fix)%*%t(R) %*% solve(R %*% solve(J_hat_real) %*% t(R)) %*% R %*% OLS_theta_real_vec_fix
Q_s_real_pv_fix <- 1 - pchisq(Q_s_real_fix, df=nrow(R))

#fit <- VAR(var_real, p = 1)
# Retain the original aod diagnostic as well as the manually calculated Wald statistics.
wald_test_func <- wald.test(b = OLS_theta_real_vec, Sigma = solve(J_hat_real), L = R)   # The diagnostic statistic requires multiplication by the sample size T.


# Calculate empirical p-values from random-weight Wald statistics (B = 199).
# Run all 199 original bootstrap iterations in the original order.
for (k in 1:rep_WB_199){

#------------------------------------------------------------#
# 3. Generate random weights {w}
#------------------------------------------------------------#

  block_w_be     <- matrix(0, nrow = nrow(resi_OLS_real), ncol = 1)  # Initialize the two-point (Bernoulli-type) weights {w_t}.
  block_w_ex     <- matrix(0, nrow = nrow(resi_OLS_real), ncol = 1)  # Initialize the standard exponential weights {w_t}.
  block_number   <- floor(nrow(resi_OLS_real)/block_size)            # Number of complete blocks; any incomplete trailing block is excluded.

  # Generate weights for each complete block.
  for (i in 1:block_number) {
    start_index  <- (i - 1) * (block_size) + 1
    end_index    <- i * block_size

    # Bernoulli-type distribution random variables
    # Draw the two-point weight first, even though the exponential weight is selected below.
    random_sign_be <- sample(c((3-sqrt(5))/2, (3+sqrt(5))/2), 1, prob = c((1+sqrt(5))/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5))))

    # standard exponential distribution random variables
    # Draw the exponential weight second. This ordering is essential to RNG reproduction.
    random_sign_ex <- rexp(1, rate = 1)

    # Assign each generated weight to the observations in its block.
    block_w_be[start_index:end_index, ] <- random_sign_be
    block_w_ex[start_index:end_index, ] <- random_sign_ex
  }

  # Select the random-weight distribution.
  #block_w <- block_w_be
  # Keep the original exponential-weight selection and the inactive alternative above.
  block_w <- block_w_ex

  #------------------------------------------------------------#
  # 4. Calculate bootstrap estimates, residuals, and the residual covariance matrix
  #------------------------------------------------------------#

  # Calculate the bootstrap coefficient estimate.
  #real
  X_real_199 <- var_real[1:number, ]                                # Construct the number-by-2 bootstrap response matrix X_t.
  X_tilde_real_199 <- var_real[1:(number - 1), ]
  X_tilde_real_199 <- rbind(c(0,0),X_tilde_real_199)                # Construct X_(t-1): prepend a zero row and omit the last observation.

  Sigma_Xt_real_199        <- matrix(0, nrow = 2, ncol = 2)    # Initialize a 2-by-2 matrix.
  Sigma_Xt_tilde_real_199  <- matrix(0, nrow = 2, ncol = 2)    # Initialize a 2-by-2 matrix.

  # Accumulate w_t * X_t %*% t(X_tilde) over all observations.
  for (t in 1:number) {
    xt_real_199              <- matrix(X_real_199[t, ], ncol = 1)                                                  # 2x1
    xt_tilde_real_199        <- matrix(X_tilde_real_199[t, ], ncol = 1)                                            # 2x1
    Sigma_Xt_real_199        <- Sigma_Xt_real_199 + block_w[t] * (xt_real_199 %*% t(xt_tilde_real_199))              # 2x2
    Sigma_Xt_tilde_real_199  <- Sigma_Xt_tilde_real_199 + block_w[t] * (xt_tilde_real_199 %*% t(xt_tilde_real_199))  # 2x2
  }
  #real
  Sigma_Xt_real_199        <- Sigma_Xt_real_199 / number
  Sigma_Xt_tilde_real_199  <- Sigma_Xt_tilde_real_199 / number

  # Calculate the random-weight bootstrap coefficient estimate.
  WB_theta_real_199 <- Sigma_Xt_real_199 %*% solve(Sigma_Xt_tilde_real_199)  # 2x2 WB estimate

  # Original author note: check whether the bootstrap coefficient matrix needs transposing. The supplied calculation is preserved.
  # Retain the original untransposed bootstrap fitted values without correction.
  X_hat_real_199    <- X_tilde_real_199 %*% (WB_theta_real_199)   # Calculate bootstrap fitted values using the original coefficient orientation.

  resi_WB_real_199  <- X_real_199 - X_hat_real_199         # Construct the number-by-2 bootstrap residual matrix.

  # Calculate the 2-by-2 contemporaneous bootstrap residual covariance using u_t %*% t(u_t).
  Omega_WB_real_199 <- list()

  S_real_199 <- matrix(0, nrow = ncol(resi_WB_real_199), ncol = ncol(resi_WB_real_199))  # Initialize a 2-by-2 matrix.

  for (t in 1:number) {

    ut_real_199     <- matrix(resi_WB_real_199[t, ], ncol = 1)          # Construct the 2-by-1 bootstrap residual vector u_t.
    S_real_199      <- S_real_199 + (ut_real_199 %*% t(ut_real_199))       # Accumulate the 2-by-2 bootstrap residual outer products.
  }
  # Bootstrap residual variance-covariance matrix.
  # Keep the bootstrap residual covariance calculation even where not used by the test.
  Omega_WB_real_199 <- S_real_199 / number

  #------------------------------------------------------------#
  # 5. Calculate the bootstrap Wald statistics
  #------------------------------------------------------------#

  # Initialize objects for the bootstrap Wald calculation.
  WB_theta_real_vec_199  <- as.vector((WB_theta_real_199))
  X_tiled_real_J_199     <- matrix(0, nrow = 2, ncol = 2)

  # Calculate the information matrix J and the Wald statistic.

  for (t in 1:number) {

    X_til_real_J_199    <- matrix(X_tilde_real_199[t, ], ncol = 1)                                       # Construct the 2-by-1 lagged observation vector X_(t-1).
    X_tiled_real_J_199  <- X_tiled_real_J_199 + (X_til_real_J_199 %*% t(X_til_real_J_199))    # Accumulate the 2-by-2 lagged cross-product matrix.
  }

  #real
  # Preserve this intermediate matrix; the original statistics below use J_hat_real.
  J_hat_real_199      <- (number^(-1)*X_tiled_real_J) %x% solve(Omega_WB_real_199)
  Q_WB_s_real_199     <- number * t((WB_theta_real_vec_199-OLS_theta_real_vec))%*%t(R) %*% solve(R %*% solve(J_hat_real) %*% t(R)) %*% R %*% (WB_theta_real_vec_199-OLS_theta_real_vec)
  Q_WB_s_real_199_fix <- number * t((WB_theta_real_vec_199-OLS_theta_real_vec_fix))%*%t(R) %*% solve(R %*% solve(J_hat_real) %*% t(R)) %*% R %*% (WB_theta_real_vec_199-OLS_theta_real_vec_fix)

  # Wald test
  # Retain strict greater-than comparisons; ties do not count as exceedances.
  Q_wald_WB_real_199[k]     <- ifelse(Q_WB_s_real_199 > Q_s_real, 1 , 0)
  Q_wald_WB_real_199_fix[k] <- ifelse(Q_WB_s_real_199_fix > Q_s_real_fix, 1 , 0)
}


# Calculate empirical p-values from B = 199 replications.
# Original empirical p-values: exceedance counts divided by 199, without adjustment.
Q_WB_s_pvalue_real_199     <- sum(Q_wald_WB_real_199)/rep_WB_199
Q_WB_s_pvalue_real_199_fix <- sum(Q_wald_WB_real_199_fix)/rep_WB_199


# Display the original eight result expressions in their original order.
Q_s_real
Q_s_real_pv
Q_WB_s_real_199
Q_WB_s_pvalue_real_199

Q_s_real_fix
Q_s_real_pv_fix
Q_WB_s_real_199_fix
Q_WB_s_pvalue_real_199_fix
