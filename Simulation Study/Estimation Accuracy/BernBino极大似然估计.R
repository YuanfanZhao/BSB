# 开始计时
start_time <- Sys.time()

# 导入需要的Package
library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

# ==============================================================================================================================

# 定义Simualtion需要的函数
Simulation <- function(number, size, m, Lambda) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  suppressMessages(suppressWarnings(library(numDeriv)))
  
  # 生成Bern--Bino分布随机数===============================================================
  rBernBino <- function(n, m, lambda) {
    K <- length(lambda) - 1
    B <- matrix(nrow = K + 1, ncol = n)
    for (k in 0:K) {
      B[k + 1,] <- rbeta(n, k + 1, K - k + 1)
    }
    Z <- sample(K + 1, n, replace = TRUE, prob = lambda)
    p <- B[cbind(Z, 1:n)]
    result <- rbinom(n, m, p)
    return(result)
  }
  
  # MM算法估计MLE================================================================
  loglikeli.fun <- function(lambda, x) {
    f <- function(data) {
      k <- 0:K
      c <- beta(k+1+data, K-k+1+m-data) / beta(k+1,K-k+1)
      
      result <- choose(m, data) * sum(lambda * c)
      return(result)
    }
    pdf <- sapply(x, f)
    return(sum(log(pdf)))
  }
  
  MM <- function(m, lambda, x, ep = 1e-6) {
    index <- 0
    number <- 1
    
    K <- length(lambda) - 1
    
    lambda_new <- lambda
    loglikeli_new <- loglikeli.fun(lambda, x)
    while (number <= 1000) {
      lambda <- lambda_new
      loglikeli <- loglikeli_new
      
      C <- sapply(0:K, function(k) beta(k+1+x, K-k+1+m-x) / beta(k+1, K-k+1))
      prop <- apply(C, 1, function(c) c * lambda)
      gamma.t <- apply(prop, 2, function(v) v / sum(v)) %>% rowSums()
      
      lambda_new[-1] <- gamma.t[-1] / sum(gamma.t)
      lambda_new[1] <- 1 - sum(lambda_new[-1])
      cat("Lambda:", lambda_new, "\n")
      
      loglikeli_new <- loglikeli.fun(lambda_new, x)
      cat("Loglikelihood:", loglikeli_new, "\n")
      
      if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
        index <- 1
        break
      }
      cat("Iteration number:", number, "\n")
      number <- number + 1
    }
    
    result <- list(Lambda = lambda_new, 
                   Loglikelihood = loglikeli_new,
                   number = number, 
                   index = index)
    return(result)
  }
  
  score_function <- function(params) {
    grad(function(p) loglikeli.fun(p, x = x), params)
  }
  
  # 计算海森矩阵在MLE处的值
  hessian_matrix <- function(params) {
    hessian(function(p) loglikeli.fun(p, x = x), params)
  }
  
  
  # 开始Simulation
  k <- 1
  
  result <- list(Lambda = matrix(nrow = length(Lambda), ncol = number), 
                 Loglikelihood = numeric(length = number), 
                 Iteration = numeric(length = number), 
                 CR = matrix(nrow = length(Lambda), ncol = number), 
                 Width = matrix(nrow = length(Lambda), ncol = number))
  
  while(k <= number) {
    x <- rBernBino(size, m, Lambda)
    Res <- tryCatch({
      MM(max(x), rep(1 / length(Lambda), length(Lambda)), x)
    }, error = function(e) {
      NA
    }, warning = function(w) {
      NA
    })
    
    
    
    if (!all(is.na(Res))) {
      if (Res$index == 1) {
        result[[1]][, k] <- Res$Lambda
        result[[2]][k] <- Res$Loglikelihood
        result[[3]][k] <- Res$number
        
        s <- score_function(Res$Lambda)
        K.hat <- s %*% t(s)
        H.hat <- hessian_matrix(Res$Lambda)
        Sigma <- solve(H.hat) %*% K.hat %*% solve(H.hat)
        lower.bound <- Res$Lambda - qnorm(0.975) * sqrt(diag(Sigma))
        upper.bound <- Res$Lambda + qnorm(0.975) * sqrt(diag(Sigma))
        result[[4]][, k] <- ifelse(lower.bound <= Lambda & upper.bound >= Lambda, 1, 0)
        result[[5]][, k] <- upper.bound - lower.bound
        
        k <- k + 1
      }
    }
  }
  
  return(result)
}

# ==============================================================================================================================

# 设定参数
num_cores <- detectCores()
all_times <- num_cores * 5
size <- 1000
m <- 8
# Lambda <- c(0.2, 0.5, 0.3)
K <- 5
Lambda <- {x <- rexp(K + 1); x/sum(x)}
# ==============================================================================================================================

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), size, m, Lambda)
stopCluster(cl)


# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)
# ==============================================================================================================================

# 汇总拟合结果并保存
Result_Lambda <- result_part[[1]]$Lambda
Result_Loglikelihood <- result_part[[1]]$Loglikelihood
Result_Iteration <- result_part[[1]]$Iteration

for (i in 2:num_cores) {
  Result_Lambda <- c(Result_Lambda, result_part[[i]]$Lambda)
  Result_Loglikelihood <- c(Result_Loglikelihood, result_part[[i]]$Loglikelihood)
  Result_Iteration <- c(Result_Iteration, result_part[[i]]$Iteration)
}

result <- list(Alpha = Result_Alpha, 
               Beta = Result_Beta, 
               Lambda = Result_Lambda,
               Var.tau = Result_Var.tau, 
               Loglikelihood = Result_Loglikelihood, 
               Iteration = Result_Iteration)

# ==============================================================================================================================

# 计算A-MLE和A-Moment
MLE.Alpha <- mean(result$Alpha)
MLE.Beta <- apply(result$Beta, c(1, 2), mean)
MLE.Lambda <- mean(result$Lambda)
MLE.Var.tau <- mean(result$Var.tau)
MLE.Loglikelihood <- mean(result$Loglikelihood)
MLE.Iteration <- mean(result$Iteration)

# ==============================================================================================================================

# 计算MSE
MLE.Alpha.MSE <- sum((result$Alpha - alpha)^2 / all_times)
MLE.Lambda.MSE <- sum((result$Lambda - lambda)^2 / all_times)
MLE.Var.tau.MSE <- sum((result$Var.tau - Var.tau)^2 / all_times)
b <- array(dim = dim(result$Beta))
for (i in 1:all_times) {
  b[, , i] <- result$Beta[, , i] - Beta
}
MLE.Beta.MSE <- apply(b^2, c(1, 2), mean)
# ==============================================================================================================================

# 总结最终结果并输出
MLE.Alpha %>% round(digits = 4)
MLE.Alpha.MSE %>% round(digits = 4)
MLE.Beta %>% round(digits = 4)
MLE.Beta.MSE %>% round(digits = 4)
MLE.Lambda %>% round(digits = 4)
MLE.Lambda.MSE %>% round(digits = 4)
MLE.Var.tau %>% round(digits = 4)
MLE.Var.tau.MSE %>% round(digits = 4)
# ==============================================================================================================================

