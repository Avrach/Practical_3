# Github Repo Link: https://github.com/Avrach/Practical_3

# Name:Yuyang Zhang  ID:S2808557
# Name:Qihong Xie  ID:S2792875
# Name:Yujie Sun  ID:S2805445

# Contribution
# Yuyang Zhang completed % of the work.
# Qihong Xie completed % of the work.
# Yujie Sun completed % of the work.

#
#
#

# Initialization 初始化
#加载必要库和数据集
library(splines)
engcov_data <- read.table("D:/Programming/Practical_3/engcov.txt", header = TRUE)
K_val=80

# Task 1
evaluate_xtilde_x_s <- function(data, K){
  #从感染到死亡的最长时间范围是80天
  max_duration <- 80
  #从首次观测的前30天开始估算感染。如果时间过于靠近可能会出现记录感染却找不到死亡的情况
  lookback_start <- 30
  
  #建立julian一年中第几天和nhs死亡人数的变量
  day_year <- data$julian
  death_num <- data$nhs
  
  #n指有多少条数据
  n <- length(day_year)
  
  # f(t)的时间向量，从julian的最小值-30到julian的最大值
  t_f <- (min(day_year) - lookback_start):max(day_year)
  m <- length(t_f)
  
  #计算π(j) 即pd
  d <- 1:max_duration
  edur <- 3.151
  sdur <- .469
  pd <- dlnorm(d, edur, sdur)
  pd <- pd / sum(pd)
  
  #1.计算矩阵S
  S <- crossprod(diff(diag(K), diff=2))
  
  #2.计算感染模型矩阵X_tilde
  ORD <- 4
  min_t_f <- min(t_f)
  max_t_f <- max(t_f)
  
  #中间的K-2个节点覆盖f(t)的区间
  inner_knots <- seq(min_t_f, max_t_f, length.out=K-2)
  #计算内部节点的间隔
  knot_spacing <- inner_knots[2]-inner_knots[1]
  #添加外部节点
  outer_knots_low <- seq(min_t_f - (ORD-1) * knot_spacing, min_t_f - knot_spacing, by=knot_spacing)
  outer_knots_high <- seq(max_t_f + knot_spacing, max_t_f + (ORD-1) * knot_spacing, by=knot_spacing)
  #全部节点
  all_knots <- c(outer_knots_low, inner_knots, outer_knots_high)
  #计算感染模型矩阵X_tilde
  X_tilde <- splineDesign(knots = all_knots, x = t_f, ord = ORD)
  
  #3.计算死亡模型矩阵X
  #计算观测到死亡的时间与感染时间的差值
  lags <- sapply(t_f, function(inf_day) day_year - inf_day)
  
  #存放概率的矩阵
  P_conv <- matrix(0,n,m)
  limits_matrix <- matrix(29 + 1:n,n,m)
  
  #填充P_conv矩阵
  P_conv[(lags >= 1) & (lags <= max_duration) & (lags <= limits_matrix)] <- pd[lags[(lags >= 1) & (lags <= max_duration) & (lags <= limits_matrix)]]
  
  #计算X
  X <- P_conv %*% X_tilde
  
  return(list(X=X, X_tilde=X_tilde, S=S))
}

# Task 2
#1.计算带惩罚的负对数似然penalized_negative_log_likelihood(pnll) 省略log(yi!)
calculate_pnll <- function(gamma, y, X, S, lambda){
  beta <- exp(gamma)
  mu <- X %*% beta
  
  #计算负对数似然
  nll <- sum(mu - y * log(mu))
  
  #计算惩罚项P
  P <- lambda * crossprod(beta, S %*% beta) / 2
  
  #计算PNLL
  pnll <- nll + P
  return(pnll)
}

#2.计算PNLL的梯度
calculate_pnll_grad <- function(gamma, y, X, S, lambda) {
  beta <- exp(gamma)
  mu <- X %*% beta
  
  #计算NLL的梯度 
  w <- 1 - y / mu
  grad_NLL <- beta * (t(X) %*% w)
  
  #计算惩罚项目P的梯度 
  S_beta <- S %*% beta
  grad_P <- lambda * beta * S_beta
  
  #返回总梯度
  grad <- grad_NLL + grad_P
  return(grad)
}

#3.测试梯度
#设置用于测试的任意起始值
gamma_test <- rep(0, K_val)
lambda_test <- 1.0 #任意选择一个 lambda值进行测试
test <- evaluate_xtilde_x_s(engcov_data, K_val)
y <- engcov_data$nhs
X <- test$X
S <- test$S

#用calculate_pnll_grad计算梯度
grad_analytic <- calculate_pnll_grad(gamma = gamma_test,
                           y = y, 
                           X = X, 
                           S = S, 
                           lambda = lambda_test)

#使用中心差分法计算梯度
h <- 1e-7 #定义一个非常小的步长
grad_numeric_manual <- numeric(K_val) #初始化一个空向量来存储结果

for (i in 1:K_val) {
  # 创建 gamma + h (仅在第i个元素上)
  gamma_plus_h <- gamma_test
  gamma_plus_h[i] <- gamma_plus_h[i] + h
  
  # 创建 gamma - h (仅在第i个元素上)
  gamma_minus_h <- gamma_test
  gamma_minus_h[i] <- gamma_minus_h[i] - h
  
  # 计算 f(gamma + h)
  f_plus_h <- calculate_pnll(gamma = gamma_plus_h, y = y, X = X, S = S, lambda = lambda_test)
  
  # 计算 f(gamma - h)
  f_minus_h <- calculate_pnll(gamma = gamma_minus_h, y = y, X = X, S = S, lambda = lambda_test)
  
  # 应用中心差分公式
  grad_numeric_manual[i] <- (f_plus_h - f_minus_h) / (2 * h)
}

# 3. 比较结果
print("--- 梯度测试结果 ---")
print("使用calculate_pnll_grad计算的梯度(前5个):")
print(head(grad_analytic, 5))

print("使用中心差分法计算的梯度(前5个):")
print(head(grad_numeric_manual, 5))

# 检查差异
diffs <- grad_analytic - grad_numeric_manual
print("差异摘要 (应接近于零):")
print(summary(diffs))

# Task 3

lambda_fixed <- 5e-5         
# Task 1 function to get X, X_tilde, and S
mats <- evaluate_xtilde_x_s(engcov_data, K_val)  
X <- mats$X
X_tilde <- mats$X_tilde
S <- mats$S

# Observed death, julian
y <- engcov_data$nhs
day <- engcov_data$julian
t_f <- (min(day) - 30):max(day)   

# Task 2 function to get obj and grad
obj  <- function(g) calculate_pnll(gamma = g, y = y, X = X, S = S, lambda = lambda_fixed)
grad <- function(g) calculate_pnll_grad(gamma = g, y = y, X = X, S = S, lambda = lambda_fixed)


gamma0 <- rep(0, K_val)
# Use BFGS to optim PNLL(gamma)
fit <- optim(par = gamma0, fn = obj, gr = grad, method = "BFGS", control = list(maxit = 500))

beta_hat <- exp(fit$par)
# Death f(t)
mu_hat <- as.numeric(X %*% beta_hat)   
# Inflection f^(t)
f_hat  <- as.numeric(X_tilde %*% beta_hat)   


par(mfrow = c(2,1), mar = c(4,4,2,1))
# Spot--real death data,line--Estimated death
plot(day, y, pch=16, cex=.6, xlab="Day (2020, julian)", ylab="Deaths (NHS)",
main=paste0("Deaths vs fitted (lambda=", lambda_fixed, ")"))
lines(day, mu_hat, lwd=2)
plot(t_f, f_hat, type="l", lwd=2, xlab="Day (2020, julian)", ylab="Estimated infections f(t)", main="Estimated daily infections")

# Task 4 
                 
# A function: Fit the penalized Poisson model for a given lambda
fit_one_lambda <- function(lambda) {
  # Penalized negative log-likelihood PNLL(gamma)
  obj  <- function(g) calculate_pnll(gamma = g, y = y, X = X, S = S, lambda = lambda)
  # Gradient of PNLL with respect to gamma
  grad <- function(g) calculate_pnll_grad(gamma = g, y = y, X = X, S = S, lambda = lambda)

  gamma0 <- rep(0, K_val)

  # Minimize PNLL using BFGS
  fit <- optim(par = gamma0, fn = obj, gr = grad, method = "BFGS")

  # Convert the estimated gamma back to beta = exp(gamma) 
  beta_hat <- exp(fit$par)

  # Compute the fitted Poisson mean miu = X beta
  mu_hat   <- as.numeric(X %*% beta_hat)

  # Evaluate the log-likelihood
  ll <- sum(y * log(mu_hat) - mu_hat)

  # Return all relevant components for this lambda
  list(lambda = lambda,
       beta = beta_hat,
       mu = mu_hat,
       loglik = ll)
}

# A funxtion: Compute the edf for a given lambda
compute_edf <- function(mu, lambda) {

  w <- y / (mu^2)
  # X^T W X 
  H0   <- t(X) %*% (X * w)

  # Penalized Hessian: H lambda = X^T W X + lambda S
  Hlam <- H0 + lambda * S

  # Effective degrees of freedom: tr( H lambda^{-1} H0 )
  edf  <- sum(diag(solve(Hlam, H0)))
  edf
}

# Define a grid of lambda values on the log scale
lambda_grid <- 10^seq(-13, -7, length.out = 50)

# Storage for BIC values and fitted models
bic_vals  <- numeric(length(lambda_grid))
fits_list <- vector("list", length(lambda_grid))

n <- length(y)   # Sample size 

# Loop over all lambda values in the grid
for (i in seq_along(lambda_grid)) {

  lam <- lambda_grid[i]

  # Fit the model for this lambda
  fit_i <- fit_one_lambda(lam)
  fits_list[[i]] <- fit_i

  # Compute effective degrees of freedom
  edf_i <- compute_edf(fit_i$mu, lam)

  # Compute BIC(lambda) = -2 logLik + log(n) * edf
  bic_vals[i] <- -2 * fit_i$loglik + log(n) * edf_i
}

# Choose the lambda with the smallest BIC
best_i <- which.min(bic_vals)
best_lambda <- lambda_grid[best_i]
best_fit <- fits_list[[best_i]]

cat("Best lambda =", best_lambda, "\n")

# Plot observed vs fitted deaths
plot(day, y)
lines(day, best_fit$mu, col="red")

# Task 5

# Extract matrices and data, and prepare the time axis

# From the Task1 function we obtain three core matrices:
# X: n×K model matrix (death mean μ = X β)
# X_tilde: m×K spline basis matrix (infection curve f(t) = X_tilde β)
# S:       K×K second-difference penalty matrix (smoothness penalty)
mats <- evaluate_xtilde_x_s(engcov_data, K_val)
X <- mats$X
X_tilde <- mats$X_tilde
S <- mats$S

# y: daily NHS deaths (length n)
y <- engcov_data$nhs

# day: corresponding Julian day indices (t_i)
day <- engcov_data$julian

# t_f: time axis on which we estimate f(t)
t_f <- (min(day) - 30):max(day)

# Fix the smoothing parameter at the value selected in Task4
lambda_star <- best_lambda

# Point estimate and warm-start in the γ space
if ("par" %in% names(best_fit)) {
  # Case A: `best_fit` is the raw return from optim()  
  fit0_beta  <- exp(best_fit$par)
  gamma_last <- best_fit$par 
} else if ("beta" %in% names(best_fit)) {
  # Case B: `best_fit` is custom list from Task 4
  fit0_beta  <- best_fit$beta 
  gamma_last <- log(pmax(fit0_beta, 1e-12))
} else {
  stop("best_fit has neither $par nor $beta; check Task 4 output.")
}

#Define the weighted objective and gradient

# Weighted penalized negative log-likelihood:
pnll_wb <- function(gamma, wb) {
  beta <- exp(gamma)
  mu <- as.numeric(X %*% beta)
  
  # numerical guard to avoid log(μ) underflow
  mu <- pmax(mu, 1e-12)
  
  # Weighted negative log-likelihood
  nll <- sum(wb * (mu - y * log(mu)))
  
  # Smoothness penalty
  pen <- as.numeric(lambda_star * crossprod(beta, S %*% beta) / 2)
  
  # Objective value minimized by BFGS
  nll + pen
}

# Gradient of the weighted objective
grad_wb <- function(gamma, wb) {
  beta <- exp(gamma)
  mu <- as.numeric(X %*% beta) 
  
  # numerical guard to avoid log(μ) underflow
  mu <- pmax(mu, 1e-12)
  
  z <- wb * (1 - y / mu)
  g1 <- as.numeric(beta * (t(X) %*% z))
  g2 <- lambda_star * beta * (S %*% beta)
  
  # length-K gradient vector
  as.numeric(g1 + g2) 
}

# Bootstrap main loop: generate B curves f̂(t)

# number of bootstrap replicates
B <- 200

# number of observations (days)
n <- length(y)

# number of time points for f(t)
m <- nrow(X_tilde)

# Preallocate containers:
# f_boot: m×B, each column is one bootstrap f̂ curve
f_boot <- matrix(NA_real_, m, B)

# mu_boot: n×B, each column is the corresponding μ̂ (for plotting death bands)
mu_boot <- matrix(NA_real_, n, B)

# store BFGS convergence codes(0 = success)
conv_flag <- integer(B)

for (b in 1:B) {
  # sample with replacement n indices from {1,…,n}; tabulate counts per index
  wb <- tabulate(sample(n, replace = TRUE), n)
  
  # wrap objective/gradient for this replicate
  obj_b <- function(g) pnll_wb(g, wb)
  grad_b <- function(g) grad_wb(g, wb)
  
  # BFGS with warm start for speed
  fit_b <- try(
    optim(par = gamma_last, fn = obj_b, gr = grad_b,
          method = "BFGS", control = list(maxit = 1000)),
    silent = TRUE
  )
  if (inherits(fit_b, "try-error") || fit_b$convergence != 0) {
    # If warm start fails (rare with extreme wb), retry from zeros and allow more iters
    fit_b <- optim(par = rep(0, K_val), fn = obj_b, gr = grad_b,
                   method = "BFGS", control = list(maxit = 2000))
  }
  # Record convergence; if success, store this replicate and update warm start
  conv_flag[b] <- if (inherits(fit_b, "try-error")) 99L else fit_b$convergence
  if (conv_flag[b] == 0) {
    gamma_last <- fit_b$par
    beta_b <- exp(fit_b$par)
    
    # fitted deaths mean for this replicate
    mu_boot[,   b ] <- drop(X %*% beta_b)
    
    # fitted infection curve for this replicate
    f_boot[,   b ] <- drop(X_tilde %*% beta_b)
  }
}
# Quick diagnostic table: ideally mostly 0 (success)
print(table(convergence = conv_flag))

# Summaries: point estimate and 95% pointwise bands 

# Point estimates (at best λ) for f(t) and μ
mu_point <- drop(X %*% fit0_beta)
f_point <- drop(X_tilde %*% fit0_beta) 

# For each time point, take the 2.5% and 97.5% quantiles across B replicates
# to form a 95% pointwise band
mu_lo <- apply(mu_boot, 1, quantile, probs = 0.025, na.rm = TRUE)
mu_hi <- apply(mu_boot, 1, quantile, probs = 0.975, na.rm = TRUE)
f_lo <- apply(f_boot,  1, quantile, probs = 0.025, na.rm = TRUE)
f_hi <- apply(f_boot,  1, quantile, probs = 0.975, na.rm = TRUE)

# Task 6 (base R)

# Axis ranges: pad by ~5% at the top to avoid clipping the band or lines
xlim_death <- range(day)
ylim_death <- c(0, max(c(y, mu_hi), na.rm = TRUE) * 1.05)

xlim_inf <- range(t_f)
ylim_inf <- c(0, max(f_hi, na.rm = TRUE) * 1.05)

# two panels, readable axes
op <- par(mfrow = c(2,1), mar = c(4.5,4.8,2.2,1.2), las = 1)

# Top panel: observed deaths (points) + fitted mean (line) + 95% shaded band
plot(day, y, pch = 16, cex = 0.6,
     xlab = "Day of 2020 (Julian)", ylab = "Deaths (NHS)",
     xlim = xlim_death, ylim = ylim_death,
     main = sprintf("Daily deaths with model fit (B=%d)", B))

# Draw band first so it stays behind the line and points
xx <- c(day, rev(day))
yy <- c(pmax(mu_lo, 0), rev(pmax(mu_hi, 0)))
polygon(xx, yy, col = rgb(0.2,0.5,0.9,0.20), border = NA)

# Fitted mean (solid line) and observed points on top
lines(day, mu_point, lwd = 2, col = "#1f77b4")
points(day, y, pch = 16, cex = 0.6)

legend("topleft",
       legend = c("Observed deaths", "Fitted mean", "95% band"),
       pch = c(16, NA, 15), pt.cex = c(0.7, NA, 1.2),
       lty = c(NA, 1, NA), lwd = c(NA, 2, NA),
       col = c("black", "#1f77b4", rgb(0.2,0.5,0.9,0.20)),
       bty = "n")

# Bottom panel: f(t) estimate (line) + 95% shaded band
plot(t_f, f_point, type = "n",
     xlab = "Day of 2020 (Julian)", ylab = "Estimated infections f(t)",
     xlim = xlim_inf, ylim = ylim_inf,
     main = "Estimated daily infections with 95% bootstrap band")

xx <- c(t_f, rev(t_f))
yy <- c(pmax(f_lo, 0), rev(pmax(f_hi, 0)))
polygon(xx, yy, col = rgb(0.8,0.2,0.2,0.20), border = NA)

lines(t_f, f_point, lwd = 2, col = "#d62728")

legend("topleft",
       legend = c("f(t) estimate", "95% band"),
       lty = c(1, NA), lwd = c(2, NA),
       pch = c(NA, 15),
       col = c("#d62728", rgb(0.8,0.2,0.2,0.20)),
       bty = "n")

# restore par settings
par(op) 

