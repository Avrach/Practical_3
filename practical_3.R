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