# 开始计时
start_time <- Sys.time()

# 导入需要的Package
library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

# ==============================================================================================================================

# 定义Simualtion需要的函数
Simulation <- function(number, size, m, Lambda, G) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  
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
      # cat("Lambda:", lambda_new, "\n")
      
      loglikeli_new <- loglikeli.fun(lambda_new, x)
      # cat("Loglikelihood:", loglikeli_new, "\n")
      
      if (abs((loglikeli_new - loglikeli) / loglikeli) < ep) {
        index <- 1
        break
      }
      # cat("Iteration number:", number, "\n")
      number <- number + 1
    }
    
    result <- list(Lambda = lambda_new, 
                   Loglikelihood = loglikeli_new,
                   number = number, 
                   index = index)
    return(result)
  }
  
  mean.std.CI <- function(lasample) {
    G <- dim(lasample)[1]
    lasort <- apply(lasample, 2, sort)
    indexx <- floor(c(0.025 * G, 0.975 * G))
    laL <- (lasort[indexx[1], ] + lasort[indexx[1] + 1, ]) / 2
    laU <- (lasort[indexx[2], ] + lasort[indexx[2] + 1, ]) / 2
    results <- cbind(laL, laU)
    return(results)
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
      MM(m, rep(1 / length(Lambda), length(Lambda)), x)
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
        
        g <- 1
        
        boot.estimation <- matrix(nrow = length(Lambda), ncol = G)
        
        while (g <= G) {
          boot.sample <- rBernBino(size, m, Res$Lambda)
          Res.boot <- tryCatch({
            MM(m, rep(1 / length(Lambda), length(Lambda)), boot.sample)
          }, error = function(e) {
            NA
          }, warning = function(w) {
            NA
          })
          
          if (!all(is.na(Res.boot))) {
            boot.estimation[, g] <- Res.boot$Lambda
            # cat("g:", g, "\n")
            g <- g + 1
          }
        }
        BCI <- mean.std.CI(t(boot.estimation))
        
        
        
        result[[4]][, k] <- ifelse(BCI[, 1] <= Lambda & BCI[, 2] >= Lambda, 1, 0)
        result[[5]][, k] <- BCI[, 2] - BCI[, 1]
        
        k <- k + 1
      }
    }
  }
  
  return(result)
}

# ==============================================================================================================================

# 设定参数
num_cores <- detectCores()
all_times <- num_cores * 11
size <- 800
m <- 15
G <- 500
Lambda <- c(0.3, 0.15, 0.1, 0.1, 0.15, 0.2)
# K <- 5
# Lambda <- {x <- rexp(K + 1); x/sum(x)}
# Simulation(1, size, m, Lambda, G)
# ==============================================================================================================================

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Simulation(number = ceiling(all_times / num_cores), size, m, Lambda, G)
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