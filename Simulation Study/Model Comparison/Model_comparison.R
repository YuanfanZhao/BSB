# 开始计时
start_time <- Sys.time()

# 导入需要的Package
library(doParallel)
library(foreach)
library(tidyverse)
library(abind)

Model.comparison <- function(number, size, m, params, Ratio) {
  suppressMessages(suppressWarnings(library(tidyverse)))
  
  # 生成数据
  generate_beta_binomial <- function(n, m, alpha, beta) {
    tau <- rbeta(n, alpha, beta)
    x <- rbinom(n, size = m, prob = tau)
    return(x)
  }
  generate_LN_binomial <- function(n, m, mu, sigma) {
    logit_tau <- rnorm(n, mean = mu, sd = sigma)
    tau <- exp(logit_tau) / (1 + exp(logit_tau))
    x <- rbinom(n, size = m, prob = tau)
    return(x)
  }
  generate_LS_binomial <- function(n, m, mu, sigma) {
    u <- runif(n)
    tau <- 1 / (1 + ((1 - u) / u)^(1 / sigma) * exp(mu / sigma))
    x <- rbinom(n, size = m, prob = tau)
    return(x)
  }
  generate_Kuma_binomial <- function(n, m, a, b) {
    u <- runif(n)
    tau <- (1 - (1 - u)^(1 / b))^(1 / a)
    x <- rbinom(n, size = m, prob = tau)
    return(x)
  }
  
  # 每个分布的密度，似然函数及求解MLE
  # Bernstein prior
  bern.den <- function(lambda, x, m) {
    K <- length(lambda) - 1
    term <- sapply(0:K, function(k) beta(k+1+x, K-k+1+m-x) / beta(k+1, K-k+1))
    choose(m, x) * sum(lambda * term)
  }
  
  bern.ll <- function(lambda, x, m) {
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
  
  bern.mle <- function(m, lambda, x, ep = 1e-8) {
    index <- 0
    number <- 1
    
    K <- length(lambda) - 1
    
    lambda_new <- lambda
    loglikeli_new <- bern.ll(lambda, x, m)
    while (number <= 1000) {
      lambda <- lambda_new
      loglikeli <- loglikeli_new
      
      C <- sapply(0:K, function(k) beta(k+1+x, K-k+1+m-x) / beta(k+1, K-k+1))
      prop <- apply(C, 1, function(c) c * lambda)
      gamma.t <- apply(prop, 2, function(v) v / sum(v)) %>% rowSums()
      
      lambda_new[-1] <- gamma.t[-1] / sum(gamma.t)
      lambda_new[1] <- 1 - sum(lambda_new[-1])
      # cat("Lambda:", lambda_new, "\n")
      
      loglikeli_new <- bern.ll(lambda_new, x, m)
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
  
  # Beta prior
  bb.den <- function(params, x, m) {
    alpha <- params[1]
    beta <- params[2]
    choose(m, x) * beta(alpha + x, m - x + beta) / beta(alpha, beta)
  }
  
  bb.ll <- function(params, x, m) {
    alpha <- params[1]
    beta <- params[2]
    
    if (alpha <= 0 || beta <= 0) return(-Inf)
    
    ll <- 0
    for (xi in x) {
      integrand <- function(tau) {
        dbeta(tau, alpha, beta) * tau^xi * (1 - tau)^(m - xi)
      }
      prob_xi <- integrate(integrand, lower = 0, upper = 1)$value
      ll <- ll + log(choose(m, xi)) + log(prob_xi)
    }
    return(ll)
  }
  
  bb.mle <- function(x, m) {
    neg_loglik <- function(par) -bb.ll(par, x, m)
    
    result <- optim(
      par = c(1, 1),
      fn = neg_loglik,
      method = "L-BFGS-B",
      lower = c(0.01, 0.01),
      hessian = TRUE
    )
    
    est <- result$par
    names(est) <- c("alpha", "beta")
    
    return(list(est = est, 
                loglikelihood = bb.ll(est, x, m)))
  }
  
  # Logit-normal
  f_tau_LN <- function(tau, mu, sigma){
    exp(-(log((tau/(1-tau)))-mu)^2/(2*sigma^2)) / (sqrt(2*pi)*sigma*tau*(1-tau))
  }
  
  ln.den <- function(params, x, m) {
    mu <- params[1]
    sigma <- params[2]
    int <- function(tau) {
      f_tau_LN(tau, mu, sigma) * tau^x * (1 - tau)^(m - x)
    }
    choose(m, x) * integrate(int, lower = 0, upper = 1)$value
  }
  
  ln.ll <- function(params, x, m) {
    mu <- params[1]
    sigma <- params[2]
    
    if (sigma <= 0) return(-Inf)
    
    ll <- 0
    for (xi in x) {
      integrand <- function(tau) {
        f_tau_LN(tau, mu, sigma) * tau^xi * (1 - tau)^(m - xi)
      }
      prob_xi <- integrate(integrand, lower = 0, upper = 1)$value
      ll <- ll + log(choose(m, xi)) + log(prob_xi)
    }
    return(ll)
  }
  
  ln.mle <- function(x, m) {
    neg_loglik <- function(par) -ln.ll(par, x, m)
    
    result <- optim(
      par = c(1, 1),
      fn = neg_loglik
    )
    
    est <- result$par
    names(est) <- c("mu", "sigma")
    
    return(list(est = est, 
                loglikelihood = ln.ll(est, x, m)))
  }
  
  # Logistic-Stdlogistic prior
  f_tau_LS <- function(tau, mu, sigma) {
    tau <- pmin(pmax(tau, 1e-6), 1 - 1e-6)
    num <- sigma * exp(-mu) * (tau * (1 - tau))^(sigma - 1)
    denom <- (tau^sigma + (1 - tau)^sigma * exp(-mu))^2
    num / denom
  }
  
  ls.den <- function(params, x, m) {
    mu <- params[1]
    sigma <- params[2]
    int <- function(tau) {
      f_tau_LS(tau, mu, sigma) * tau^x * (1 - tau)^(m - x)
    }
    choose(m, x) * integrate(int, lower = 0, upper = 1)$value
  }
  
  ls.ll <- function(params, x, m) {
    mu <- params[1]
    sigma <- params[2]
    
    if (sigma <= 0) return(-Inf)
    
    ll <- 0
    for (xi in x) {
      integrand <- function(tau) {
        f_tau_LS(tau, mu, sigma) * tau^xi * (1 - tau)^(m - xi)
      }
      prob_xi <- integrate(integrand, lower = 0, upper = 1)$value
      ll <- ll + log(choose(m, xi)) + log(prob_xi)
    }
    return(ll)
  }
  
  ls.mle <- function(x, m) {
    neg_loglik <- function(par) -ls.ll(par, x, m)
    
    result <- optim(
      par = c(1, 1),
      fn = neg_loglik
    )
    
    est <- result$par
    names(est) <- c("mu", "sigma")
    
    return(list(est = est, 
                loglikelihood = ls.ll(est, x, m)))
  }
  
  # Kumaraswany prior
  f_tau_Kuma <- function(tau, a, b) {
    a * b * tau^(a - 1) * (1 - tau^a)^(b - 1)
  }
  
  km.den <- function(params, x, m) {
    a <- params[1]
    b <- params[2]
    int <- function(tau) {
      f_tau_Kuma(tau, a, b) * tau^x * (1 - tau)^(m - x)
    }
    choose(m, x) * integrate(int, lower = 0, upper = 1)$value
  }
  
  km.ll <- function(params, x, m) {
    a <- params[1]
    b <- params[2]
    
    if (a <= 0 || b <= 0) return(-Inf)
    
    ll <- 0
    for (xi in x) {
      integrand <- function(tau) {
        f_tau_Kuma(tau, a, b) * tau^xi * (1 - tau)^(m - xi)
      }
      prob_xi <- integrate(integrand, lower = 0, upper = 1)$value
      ll <- ll + log(choose(m, xi)) + log(prob_xi)
    }
    return(ll)
  }
  
  km.mle <- function(x, m) {
    neg_loglik <- function(par) -km.ll(par, x, m)
    
    result <- optim(
      par = c(1, 1),
      fn = neg_loglik,
      method = "L-BFGS-B",
      lower = c(0.01, 0.01),
      hessian = TRUE
    )
    
    est <- result$par
    names(est) <- c("a", "b")
    
    return(list(est = est, 
                loglikelihood = km.ll(est, x, m)))
  }
  
  num <- 1
  Result <- array(dim = c(5, 5, number), dimnames = list(c("Beta", "LN", "LS", "KM", "Bern"), c("Loglikelihood", "AIC", "Rank1", "P.value", "Rank2"), 1:number))
  
  a1 <- params$Beta[1]
  b1 <- params$Beta[2]
  a2 <- params$LN[1]
  b2 <- params$LN[2]
  a3 <- params$LS[1]
  b3 <- params$LS[2]
  a4 <- params$KM[1]
  b4 <- params$KM[2]
  
  while (num <= number) {
    x1 <- generate_beta_binomial(size, m, a1, b1)
    x2 <- generate_LN_binomial(size, m, a2, b2)
    x3 <- generate_LS_binomial(size, m, a3, b3)
    x4 <- generate_Kuma_binomial(size, m, a4, b4)
    z <- rmultinom(size, 1, Ratio)
    x <- x1 * z[1, ] + x2 * z[2, ] + x3 * z[3, ] + x4 * z[4, ]
    
    obser.x <- rep(0, m + 1)
    obser.x <- sapply(0:m, function(i) obser.x[i+1] <- sum(x == i))
    
    # 各分布拟合
    bb.fit <- bb.mle(x, m)
    ln.fit <- ln.mle(x, m)
    ls.fit <- ls.mle(x, m)
    km.fit <- km.mle(x, m)
    
    K <- 1:20
    candidate <- matrix(0, nrow = length(K), ncol = 2)
    candidate <- t(
      vapply(K, function(k) {
        bern.fit <- bern.mle(m, rep(1/(k+1), k+1), x)
        bern.density <- sapply(0:m, function(x) bern.den(bern.fit$Lambda, x, m))
        bern.chi <- chisq.test(x = obser.x, p = bern.density, simulate.p.value = TRUE)$p.value
        # cat("K:", k, "\n")
        
        c(K = k, AIC = 2 * (k + 1 - bern.fit$Loglikelihood), p.value = bern.chi)
      }, numeric(3))
    )
    
    # 卡方拟合优度检验
    bb.density <- sapply(0:m, function(x) bb.den(bb.fit$est, x, m))
    bb.chi <- chisq.test(x = obser.x, p = bb.density, simulate.p.value = TRUE)$p.value
    
    ln.density <- sapply(0:m, function(x) ln.den(ln.fit$est, x, m))
    ln.density <- ln.density / sum(ln.density)
    ln.chi <- chisq.test(x = obser.x, p = ln.density, simulate.p.value = TRUE)$p.value
    
    ls.density <- sapply(0:m, function(x) ls.den(ls.fit$est, x, m))
    ls.density <- ls.density / sum(ls.density)
    ls.chi <- chisq.test(x = obser.x, p = ls.density, simulate.p.value = TRUE)$p.value
    
    km.density <- sapply(0:m, function(x) km.den(km.fit$est, x, m))
    km.density <- km.density / sum(km.density)
    km.chi <- chisq.test(x = obser.x, p = km.density, simulate.p.value = TRUE)$p.value
    
    Result[, 1, num] <- c(bb.fit$loglikelihood, ln.fit$loglikelihood, ls.fit$loglikelihood, km.fit$loglikelihood, candidate[which.min(candidate[, 2]), 1] + 1 - candidate[which.min(candidate[, 2]), 2] / 2) %>% round(digits = 2)
    Result[, 2, num] <- c(4 - 2 * bb.fit$loglikelihood, 4 - 2 * ln.fit$loglikelihood, 4 - 2 * ls.fit$loglikelihood, 4 - 2 * km.fit$loglikelihood, candidate[which.min(candidate[, 2]), 2]) %>% round(digits = 2)
    Result[, 3, num] <- rank(-Result[, 2, num])
    Result[, 4, num] <- c(bb.chi, ln.chi, ls.chi, km.chi, candidate[which.min(candidate[, 2]), 3]) %>% round(digits = 4)
    Result[, 5, num] <- rank(Result[, 4, num])
    num <- num + 1
  }
  return(Result)
}

# 设定参数
num_cores <- detectCores()
all_times <- num_cores * 11
m <- 30
size <- 500
# 第一种情况
# params <- list(Beta = c(1, 3),
#                LN = c(-1, 0.5),
#                LS = c(2, 5),
#                KM = c(1, 2))
# Ratio <- c(0.25, 0.25, 0.25, 0.25)
# Ratio <- c(0.4, 0.2, 0.2, 0.2)
# Ratio <- c(0.2, 0.4, 0.2, 0.2)
# Ratio <- c(0.2, 0.2, 0.4, 0.2)
# Ratio <- c(0.2, 0.2, 0.2, 0.4)

# 第二种情况，极端好
params <- list(Beta = c(1, 0.5),
               LN = c(0, 0.25),
               LS = c(2, 5),
               KM = c(1, 5))
Ratio <- c(0.25, 0.25, 0.25, 0.25)

# 进行并行运算
cl <- makeCluster(num_cores)
registerDoParallel(cl)
result_part <- foreach(x = 1:num_cores) %dopar% Model.comparison(number = ceiling(all_times / num_cores), size, m, params, Ratio)
stopCluster(cl)

Result <- result_part[[1]]
for (i in 2:num_cores) {
  Result <- abind(Result, result_part[[i]])
}

Final <- apply(Result, c(1, 2), mean) %>% round(digits = 4)
knitr::kable(Final)

# 结束运算，输出时间
end_time <- Sys.time()
print(end_time - start_time)