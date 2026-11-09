# PRIMARY REPRODUCTION SCRIPT
# Only file paths and explanatory comments differ from the supplied original.
# Run from the repository root; see README.md for packages and execution.
# Original comments and inactive plotting code are retained for traceability.

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
# 1. real data 불러오기
#------------------------------------------------------------#
# Original RNG seed. Do not move this statement or insert random draws below it.
set.seed(165)

# txt 파일을 불러옵니다. (rawdata)
# Read the supplied merchandise-balance file without changing its contents.
bopbm_rawdata <- read.table("data/bopbm_rawdata_2.txt",header=T)
# "data/
# Preserve the original date conversion, including the appended day string.
bopbm_rawdata$date <- as.Date(paste(bopbm_rawdata$date, "-01", sep = ""), format = "%Y-%m-%d")

# # 시계열 그래프 생성
# png("plot_1.png", width=6000, height=4000, res=500, type="cairo")
# ggplot(bopbm_rawdata, aes(x = date, y = BOPBM)) +
#   geom_line() + # 선 그래프 생성
#   labs(x = "Date", y = "BOPBM") + # x축, y축 라벨 설정
#   scale_x_date(date_breaks = "5 year", date_labels = "%Y-%m") + # x축 날짜 형식 지정
#   theme_test() # 기본 테마 적용

# 1차 차분을 계산합니다.
# Original preprocessing: sort by date, difference numeric columns, then omit NA rows.
bopbm_diff <- bopbm_rawdata %>%
  arrange(date) %>%
  mutate(across(where(is.numeric), ~ c(NA, diff(.)))) %>%
  na.omit()  # 첫 번째 차분에서 발생하는 NA 행을 제거합니다.


# # 시계열 그래프 생성
# png("plot_2.png", width=6000, height=4000, res=500, type="cairo")
# ggplot(bopbm_diff, aes(x = date, y = BOPBM)) +
#   geom_line() + # 선 그래프 생성
#   labs(x = "Date", y = "BOPBM") + # x축, y축 라벨 설정
#   scale_x_date(date_breaks =  "5 year", date_labels = "%Y-%m") + # x축 날짜 형식 지정
#   theme_test() # 기본 테마 적용


# txt 파일을 불러옵니다. (rawdata)
# Read the supplied services-balance file without changing its contents.
bopbsv_rawdata <- read.table("data/bopbsv_rawdata_2.txt",header=T)
# Preserve the original date conversion exactly as supplied.
bopbsv_rawdata$date <- as.Date(paste(bopbsv_rawdata$date, "-01", sep = ""), format = "%Y-%m-%d")

# # 시계열 그래프 생성
# png("plot_3.png", width=6000, height=4000, res=500, type="cairo")
# ggplot(bopbsv_rawdata, aes(x = date, y = BOPBSV)) +
#   geom_line() + # 선 그래프 생성
#   labs(x = "Date", y = "BOPBSV") + # x축, y축 라벨 설정
#   scale_x_date(date_breaks ="5 year", date_labels = "%Y-%m") + # x축 날짜 형식 지정
#   theme_test() # 기본 테마 적용


# 1차 차분을 계산합니다.
# Apply the original services-series preprocessing without rewriting the pipeline.
bopbsv_diff <- bopbsv_rawdata %>%
  arrange(date) %>%
  mutate(across(where(is.numeric), ~ c(NA, diff(.)))) %>%
  na.omit()  # 첫 번째 차분에서 발생하는 NA 행을 제거합니다.

# # 시계열 그래프 생성
# png("plot_4.png", width=6000, height=4000, res=500, type="cairo")
# ggplot(bopbsv_diff, aes(x = date, y = BOPBSV)) +
#   geom_line() + # 선 그래프 생성
#   labs(x = "Date", y = "BOPBSV") + # x축, y축 라벨 설정
#   scale_x_date(date_breaks ="5 year", date_labels = "%Y-%m") + # x축 날짜 형식 지정
#   theme_test() # 기본 테마 적용


# bopbm_diff = 상품수지 (balance on merchandise)
# bopbsv_diff = 서비스스수지 (balance on services)

# Real Dataset 생성
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
number        <- 175                                       # 실제 데이터 개수 176
block_size    <- floor((number)^(0))                       # block size
p             <- 1                                         # 모형적합 차수 : VAR(p)
df            <- 1                                         # 자유도 : 제약식의 수 (1 또는 2)
cv            <- qchisq(0.95, df)                          # critical value
rep_WB_199    <- 199                                       # B=199
Q_wald_WB_real_199      <- matrix(0, nrow = rep_WB_199, ncol = 1)      # B번 계산된 WB-wald통계량 (real)
Q_wald_WB_real_199_fix  <- matrix(0, nrow = rep_WB_199, ncol = 1)      # B번 계산된 WB-wald통계량 (real)
# Column-major coefficient order is (a11, a21, a12, a22); the restriction selects a12.
R             <- matrix(c(0,0,1,0), nrow=1)                           # 제약식 a_12가 0이다. (R은 열 기준, a11,a21,a12,a22)
half_number   <- floor(number/2)                    # 1/2 지점


#------------------------------------------------------------#
# 2. OLS 추정량,잔차,잔차의 자기공분산행렬 구하기
#------------------------------------------------------------#

# OLS estimate
# Construct the original response and lag arrays, retaining the initial zero lag row.
X_real <- var_real[1:number, ]                            # (number)×2, t자료 생성 X_t
X_tilde_real <- var_real[1:(number - 1), ]
X_tilde_real <- rbind(c(0,0),X_tilde_real)                # (number)×2, 맨마지막 자료 제외(처음자료는 0으로 설정정) X_t-1

Sigma_Xt_real        <- matrix(0, nrow = 2, ncol = 2)     # 2*2 matrix 생성
Sigma_Xt_tilde_real  <- matrix(0, nrow = 2, ncol = 2)     # 2*2 matrix 생성


# 각 시점마다 X_t %*% t(Xt_tilde)를 계산해서 누적합 구하기
for (t in 1:number) {
  xt_real              <- matrix(X_real[t, ], ncol = 1)                             # 2x1
  xt_tilde_real        <- matrix(X_tilde_real[t, ], ncol = 1)                       # 2x1
  Sigma_Xt_real        <- Sigma_Xt_real + (xt_real %*% t(xt_tilde_real))              # 2x2
  Sigma_Xt_tilde_real  <- Sigma_Xt_tilde_real + (xt_tilde_real %*% t(xt_tilde_real))  # 2x2
}
#real
Sigma_Xt_real       <- Sigma_Xt_real / number         #sigma_xt
Sigma_Xt_tilde_real <- Sigma_Xt_tilde_real / number   #sigma_xt_tilde

# OLS 추정량 구하기
# OLS estimate using the original accumulated matrix products and inversion order.
OLS_theta_real <- Sigma_Xt_real %*% solve(Sigma_Xt_tilde_real)  # 2x2 OLS estimate

# OLS 잔차 구하기
# Reproduction convention: retain the untransposed coefficient matrix exactly as supplied.
X_hat_real <- X_tilde_real %*% (OLS_theta_real)   # (number-1)×2의 예측값 생성 (마지막자료 제외 * 추정량)

resi_OLS_real <- X_real - X_hat_real              # OLS 잔차 생성 ((number-1)×2)

# OLS 잔차의 분산-공분산행렬 구하기(2*2, 동일 시점 내 두 변수 관계) ut_real %*% t(ut_real)
Omega_OLS_real <- list()

S_real <- matrix(0, nrow = ncol(resi_OLS_real), ncol = ncol(resi_OLS_real))  # 2x2 matrix 생성

for (t in 1:number) {
  ut_real     <- matrix(resi_OLS_real[t, ], ncol = 1)            # 2×1 u_t 구하기
  S_real      <- S_real + (ut_real %*% t(ut_real))               # 2×2 u_t×u_t 분산-공분산행렬 Omega
}
# 잔차의 분산-공분산 행렬
Omega_OLS_real <- S_real / number


#------------------------------------------------------------#
# 5. 잔차를 앞 1/2, 뒤 1/2로 나누기
#------------------------------------------------------------#

# Retain the separate first-half and second-half residual covariance diagnostics.
resi_OLS_real_front <- resi_OLS_real[1:half_number, ]
resi_OLS_real_back  <- resi_OLS_real[(half_number+1):number, ]

number_front <- nrow(resi_OLS_real_front)
number_back  <- nrow(resi_OLS_real_back)


#------------------------------------------------------------#
# 6. 앞 1/2 구간 잔차의 평균 구하기
#------------------------------------------------------------#
# Calculate the first-half residual mean using the original summation loop.
mean_front <- matrix(0, nrow = 1, ncol = 2)

for (t in 1:number_front) {
  mean_front <- mean_front + matrix(resi_OLS_real_front[t, ], nrow = 1)
}
mean_front <- mean_front / number_front


#------------------------------------------------------------#
# 7. 뒤 1/2 구간 잔차의 평균 구하기
#------------------------------------------------------------#
# Calculate the second-half residual mean using the original summation loop.
mean_back <- matrix(0, nrow = 1, ncol = 2)

for (t in 1:number_back) {
  mean_back <- mean_back + matrix(resi_OLS_real_back[t, ], nrow = 1)
}
mean_back <- mean_back / number_back


#------------------------------------------------------------#
# 8. 앞 1/2 구간 분산-공분산행렬 구하기
#------------------------------------------------------------#
# Center first-half residuals by their own mean; retain the original divisor.
S_front <- matrix(0, nrow = 2, ncol = 2)

for (t in 1:number_front) {
  ut_front <- matrix(resi_OLS_real_front[t, ] - mean_front, ncol = 1)
  S_front  <- S_front + (ut_front %*% t(ut_front))
}

Omega_front <- S_front / number_front


#------------------------------------------------------------#
# 9. 뒤 1/2 구간 분산-공분산행렬 구하기
#------------------------------------------------------------#
# Center second-half residuals by their own mean; retain the original divisor.
S_back <- matrix(0, nrow = 2, ncol = 2)

for (t in 1:number_back) {
  ut_back <- matrix(resi_OLS_real_back[t, ] - mean_back, ncol = 1)
  S_back  <- S_back + (ut_back %*% t(ut_back))
}

Omega_back <- S_back / number_back




#------------------------------------------------------------#
# 2-1. OLS 추정량에 대한 Wald 통계량 구하기
#------------------------------------------------------------#

OLS_theta_real_vec <- as.vector((OLS_theta_real))
X_tiled_real_J <- matrix(0, nrow = 2, ncol = 2)
# Retain the complete fixed coefficient vector and its separate analysis.
OLS_theta_real_vec_fix <- as.vector((c(0.43,0,-1.02,-0.01)))

# J행렬과 Wald 통계량 계산하기

for (t in 1:number) {
  X_til_real_J    <- matrix(X_tilde_real[t, ], ncol = 1)                  # 2×1 X_(t-1) 구하기
  X_tiled_real_J  <- X_tiled_real_J + (X_til_real_J %*% t(X_til_real_J))    # 2×2 X_tiled_real_J 행렬
}

#real
# Original Kronecker information matrix and Wald formulas; no scalar simplification.
J_hat_real   <- (number^(-1)*X_tiled_real_J) %x% solve(Omega_OLS_real)
Q_s_real     <- number * t(OLS_theta_real_vec)%*%t(R) %*% solve(R %*% solve(J_hat_real) %*% t(R)) %*% R %*% OLS_theta_real_vec

# Wald test
Q_s_real_pv  <- 1 - pchisq(Q_s_real, df=nrow(R))

# fix 추정량 계산산
Q_s_real_fix <- number * t(OLS_theta_real_vec_fix)%*%t(R) %*% solve(R %*% solve(J_hat_real) %*% t(R)) %*% R %*% OLS_theta_real_vec_fix
Q_s_real_pv_fix <- 1 - pchisq(Q_s_real_fix, df=nrow(R))

#fit <- VAR(var_real, p = 1)
# Retain the original aod diagnostic as well as the manually calculated Wald statistics.
wald_test_func <- wald.test(b = OLS_theta_real_vec, Sigma = solve(J_hat_real), L = R)   # 통계량에서 T를 곱해줘야함


# Random weight-LB 통계량으로 경험적 p-value 구하기기 (B=199)
# Run all 199 original bootstrap iterations in the original order.
for (k in 1:rep_WB_199){

#------------------------------------------------------------#
# 3. Random weight {w} 생성하기
#------------------------------------------------------------#

  block_w_be     <- matrix(0, nrow = nrow(resi_OLS_real), ncol = 1)  # {w_t} 생성 (Bernoulli)
  block_w_ex     <- matrix(0, nrow = nrow(resi_OLS_real), ncol = 1)  # {w_t} 생성 (Std.Exponential)
  block_number   <- floor(nrow(resi_OLS_real)/block_size)            # full block의 개수 (마지막 제외)

  # full block에 대한 {w} 생성  (1st ~ (L_T-1)th) block
  for (i in 1:block_number) {
    start_index  <- (i - 1) * (block_size) + 1
    end_index    <- i * block_size

    # Bernoulli-type distribution random variables
    # Draw the two-point weight first, even though the exponential weight is selected below.
    random_sign_be <- sample(c((3-sqrt(5))/2, (3+sqrt(5))/2), 1, prob = c((1+sqrt(5))/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5))))

    # standard exponential distribution random variables
    # Draw the exponential weight second. This ordering is essential to RNG reproduction.
    random_sign_ex <- rexp(1, rate = 1)

    # 변형된 블록을 block_w에 저장
    block_w_be[start_index:end_index, ] <- random_sign_be
    block_w_ex[start_index:end_index, ] <- random_sign_ex
  }

  #random weight 설정하기
  #block_w <- block_w_be
  # Keep the original exponential-weight selection and the inactive alternative above.
  block_w <- block_w_ex

  #------------------------------------------------------------#
  # 4. Bootstrap 추정량,잔차,잔차의 자기공분산행렬 구하기
  #------------------------------------------------------------#

  # Bootstrap estimate 구하기
  #real
  X_real_199 <- var_real[1:number, ]                                # (number)×2, 맨처음 자료 제외  (t자료 생성) X_t
  X_tilde_real_199 <- var_real[1:(number - 1), ]
  X_tilde_real_199 <- rbind(c(0,0),X_tilde_real_199)                # (number)×2, 맨마지막 자료 제외(처음자료는 0으로 설정) X_t-1

  Sigma_Xt_real_199        <- matrix(0, nrow = 2, ncol = 2)    # 2*2 matrix 생성
  Sigma_Xt_tilde_real_199  <- matrix(0, nrow = 2, ncol = 2)    # 2*2 matrix 생성

  # 각 시점마다 w_t * X_t %*% t(Xt_tilde)를 계산해서 누적합 구하기
  for (t in 1:number) {
    xt_real_199              <- matrix(X_real_199[t, ], ncol = 1)                                                  # 2x1
    xt_tilde_real_199        <- matrix(X_tilde_real_199[t, ], ncol = 1)                                            # 2x1
    Sigma_Xt_real_199        <- Sigma_Xt_real_199 + block_w[t] * (xt_real_199 %*% t(xt_tilde_real_199))              # 2x2
    Sigma_Xt_tilde_real_199  <- Sigma_Xt_tilde_real_199 + block_w[t] * (xt_tilde_real_199 %*% t(xt_tilde_real_199))  # 2x2
  }
  #real
  Sigma_Xt_real_199        <- Sigma_Xt_real_199 / number
  Sigma_Xt_tilde_real_199  <- Sigma_Xt_tilde_real_199 / number

  # WB 추정량 구하기
  WB_theta_real_199 <- Sigma_Xt_real_199 %*% solve(Sigma_Xt_tilde_real_199)  # 2x2 WB estimate

  # WB 잔차 구하기 WB 추정량에 transpose를 해야하는지 확인해볼 것!! (수정필요)
  # Retain the original untransposed bootstrap fitted values without correction.
  X_hat_real_199    <- X_tilde_real_199 %*% (WB_theta_real_199)   # 199×2의 예측값 생성 (마지막자료 제외 * 추정량)

  resi_WB_real_199  <- X_real_199 - X_hat_real_199         # WB 잔차 생성 (199×2)

  # WB 잔차의 분산-공분산행렬 구하기(2*2, 동일 시점 내 두 변수 관계) ut_real %*% t(ut_real)
  Omega_WB_real_199 <- list()

  S_real_199 <- matrix(0, nrow = ncol(resi_WB_real_199), ncol = ncol(resi_WB_real_199))  # 2x2 matrix 생성

  for (t in 1:number) {

    ut_real_199     <- matrix(resi_WB_real_199[t, ], ncol = 1)          # 2×1 u_t 구하기
    S_real_199      <- S_real_199 + (ut_real_199 %*% t(ut_real_199))       # 2×2 u_t×u_t 분산-공분산행렬 Omega
  }
  # WB 잔차의 분산-공분산 행렬
  # Keep the bootstrap residual covariance calculation even where not used by the test.
  Omega_WB_real_199 <- S_real_199 / number

  #------------------------------------------------------------#
  # 5. WB-Wald 통계량 구하기
  #------------------------------------------------------------#

  # WB-Wald 통계량 초기화
  WB_theta_real_vec_199  <- as.vector((WB_theta_real_199))
  X_tiled_real_J_199     <- matrix(0, nrow = 2, ncol = 2)

  # J행렬과 Wald 통계량 계산하기

  for (t in 1:number) {

    X_til_real_J_199    <- matrix(X_tilde_real_199[t, ], ncol = 1)                                       # 2×1 X_(t-1) 구하기
    X_tiled_real_J_199  <- X_tiled_real_J_199 + (X_til_real_J_199 %*% t(X_til_real_J_199))    # 2×2 X_tiled_real_J 행렬
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


# B=199번
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
