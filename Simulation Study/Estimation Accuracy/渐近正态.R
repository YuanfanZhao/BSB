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
    K <- length(lambda) - 1
    f <- function(data) {
      k <- 0:K
      c <- beta(k+1+data, K-k+1+m-data) / beta(k+1,K-k+1)
      
      result <- choose(m, data) * sum(lambda * c)
      return(result)
    }
    pdf <- sapply(x, f)
    return(sum(log(pdf)))
  }
  
  loglikeli_theta <- function(theta, x) {
    K <- length(theta)
    # 从θ重建完整的λ: λ0 = 1 - sum(θ), λ[1:K] = θ
    lambda <- c(1 - sum(theta), theta)
    return(loglikeli.fun(lambda, x))
  }
  
  MM <- function(m, lambda, x, ep = 1e-8) {
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
  
  # 计算θ的得分函数贡献 (每个观测的梯度)
  score_contributions_theta <- function(theta, x) {
    K <- length(theta)
    n <- length(x)
    
    # 使用数值梯度计算每个观测的得分贡献
    score_mat <- matrix(0, nrow = n, ncol = K)
    
    for(i in 1:n) {
      # 为每个观测计算梯度
      grad_i <- numDeriv::grad(
        func = function(th) {
          # 计算单个观测的对数似然
          lambda <- c(1 - sum(th), th)
          k <- 0:K
          c_val <- beta(k+1+x[i], K-k+1+m-x[i]) / beta(k+1, K-k+1)
          log_val <- log(choose(m, x[i]) * sum(lambda * c_val))
          return(log_val)
        },
        theta,
        method = "Richardson",
        method.args = list(eps = 1e-10, d = 0.0001, r = 6)
      )
      score_mat[i, ] <- grad_i
    }
    
    return(score_mat)
  }
  
  # 计算θ的海森矩阵
  hessian_theta <- function(theta, x) {
    K <- length(theta)
    
    H <- numDeriv::hessian(
      func = function(th) loglikeli_theta(th, x = x),
      theta,
      method = "Richardson",
      method.args = list(eps = 1e-10, d = 0.0001, r = 6)
    )
    
    return(H)
  }
  
  # 使用Delta方法计算λ0的方差
  delta_method_variance <- function(theta, Sigma_theta) {
    # g(θ) = 1 - sum(θ)
    # 梯度: ∇g(θ) = (-1, -1, ..., -1)
    grad_g <- rep(-1, length(theta))
    
    # Delta方法方差: Var(g(θ)) ≈ [∇g(θ)]^T Σ(θ) [∇g(θ)]
    var_lambda0 <- as.numeric(t(grad_g) %*% Sigma_theta %*% grad_g)
    
    return(var_lambda0)
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
        
        theta_hat <- Res$Lambda[-1]
        # 计算三明治协方差矩阵
        # 1. 计算得分贡献矩阵
        S <- score_contributions_theta(theta_hat, x)
        B_matrix <- crossprod(S)
        H <- -hessian_theta(theta_hat, x)
        Sigma_theta <- solve(H) %*% B_matrix %*% solve(H)
        var_lambda0 <- delta_method_variance(theta_hat, Sigma_theta)
        se <- c(var_lambda0, diag(Sigma_theta)) %>% sqrt()
        
        lower.bound <- Res$Lambda - qnorm(0.975) * se
        upper.bound <- Res$Lambda + qnorm(0.975) * se
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
all_times <- num_cores * 105
size <- 200
m <- 15
Lambda <- c(0.15, 0.3, 0.35, 0.2)
# K <- 5
# Lambda <- {x <- rexp(K + 1); x/sum(x)}
# ==============================================================================================================================

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), size, m, Lambda)
stopCluster(cl)

# ==============================================================================================================================

# 汇总拟合结果并保存
Result_Lambda <- result_part[[1]]$Lambda
Result_Loglikelihood <- result_part[[1]]$Loglikelihood
Result_Iteration <- result_part[[1]]$Iteration
Result_CR <- result_part[[1]]$CR
Result_Width <- result_part[[1]]$Width

for (i in 2:num_cores) {
  Result_Lambda <- cbind(Result_Lambda, result_part[[i]]$Lambda)
  Result_Loglikelihood <- c(Result_Loglikelihood, result_part[[i]]$Loglikelihood)
  Result_Iteration <- c(Result_Iteration, result_part[[i]]$Iteration)
  Result_CR <- cbind(Result_CR, result_part[[i]]$CR)
  Result_Width <- cbind(Result_Width, result_part[[i]]$Width)
}

result <- list(Lambda = Result_Lambda,
               Loglikelihood = Result_Loglikelihood, 
               Iteration = Result_Iteration, 
               CR = Result_CR, 
               Width = Result_Width)

# ==============================================================================================================================

# 计算A-MLE
Lambda.MLE <- rowMeans(result$Lambda)
MLE.Loglikelihood <- mean(result$Loglikelihood)
MLE.Iteration <- mean(result$Iteration)

AIC.value <- 2 * (length(Lambda) - MLE.Loglikelihood)

# ==============================================================================================================================

# 计算MSE
Lambda.MSE <- apply(result$Lambda, 2, function(s) (s - Lambda)^2) %>% rowMeans()

# 计算CR与Width
CR <- rowMeans(result$CR)
Width <- rowMeans(result$Width)
# ==============================================================================================================================

# 总结最终结果并输出
Lambda.MLE %>% round(digits = 4)
Lambda.MSE %>% round(digits = 4)
AIC.value %>% round(digits = 4)
CR %>% round(digits = 4)
Width %>% round(digits = 4)

# ==============================================================================================================================

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)