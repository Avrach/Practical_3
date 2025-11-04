# Github Link: https://github.com/Avrach/Practical_3

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
engcov_data <- read.table("engcov.txt", header = TRUE)

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
  
}