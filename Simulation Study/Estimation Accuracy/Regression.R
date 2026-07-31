lambda.fun <- function(z, beta) {
  eta <- z %*% beta
  eta_full <- cbind(0, eta)
  exp_eta <- exp(eta_full)
  lambda <- exp_eta / rowSums(exp_eta)
  
  return(lambda)
}

rBernBino <- function(n, m, z, beta) {
  lambda <- lambda.fun(z, beta)
  
  rbb <- function(n, m, lambda) {
    K <- length(lambda) - 1
    B <- sapply(0:K, function(k) rbeta(n, k + 1, K - k + 1))
    Z <- sample(K + 1, n, replace = TRUE, prob = lambda)
    p <- B[cbind(1:n, Z)]
    result <- rbinom(n, m, p)
    return(result)
  }
  
  x <- apply(lambda, 1, function(la) {rbb(1, m, la)})
  
  return(x)
}

n <- 1000
z <- cbind(1, rexp(n, 1), rnorm(n, 0, 2))
beta <- matrix(c(0.5, 0, -1, 3, -2, 5, -0.6, 0.5, 1), nrow = 3, ncol = 3, byrow = TRUE)
m <- 10

x <- rBernBino(n, m, z, beta)

Regression <- function(x, z, m, beta, ep = 1e-6, maxit = 1000) {
  loglikeli.fun <- function(x, z, beta, m) {
    n <- length(x)
    K <- ncol(beta)
    
    lambda <- lambda.fun(z, beta)
    
    cmat <- matrix(0, n, K + 1)
    for (k in 0:K) {
      cmat[, k + 1] <- beta(x + k + 1, m - x + (K - k) + 1) / beta(k + 1, K - k + 1)
    }
    
    inside <- rowSums(cmat * lambda)
    ll <- sum(log(choose(m, x) * inside))
    return(ll)
  }
  
  gamma.fun <- function(x, z, beta, m) {
    n <- length(x)
    K <- ncol(beta)
    
    lambda <- lambda.fun(z, beta)
    
    cmat <- matrix(0, n, K + 1)
    for (k in 0:K) {
      cmat[, k+1] <- beta(x + k + 1, m - x + (K - k) + 1) / beta(k + 1, K - k + 1)
    }
    
    denom <- rowSums(cmat * lambda)
    gat <- (cmat * lambda) / denom
    
    return(gat)
  }
  
  Q.fun <- function(x, z, beta, gat, m) {
    lambda <- lambda.fun(z, beta)
    Qval <- sum(gat * log(lambda))
    return(Qval)
  }
  
  IRLS <- function(gat, z, beta_ini, ep = 1e-6, maxit = 1000) {
    n <- nrow(z)
    p <- ncol(z)
    K <- ncol(beta_ini)  # 类别数（不含基类0）
    
    beta <- beta_ini
    Q.val <- Q.fun(x, z, beta, gat, m)
    
    for (iter in 1:maxit) {
      # 1. eta & lambda
      eta_mat <- z %*% beta         # n x K
      eta_full <- cbind(0, eta_mat) # n x (K+1)
      row_max <- apply(eta_full, 1, max)
      exp_eta <- exp(eta_full - row_max)
      lambda_full <- exp_eta / rowSums(exp_eta)  # n x (K+1)
      lambda_no0 <- lambda_full[, -1, drop=FALSE]  # n x K
      
      # 2. 构造设计矩阵 Ztilde (nK x pK)
      Ztilde <- matrix(0, n*K, p*K)
      for (k in 1:K) {
        row_start <- (k - 1) * n + 1
        row_end <- k * n
        col_start <- (k - 1) * p + 1
        col_end <- k * p
        Ztilde[row_start:row_end, col_start:col_end] <- z
      }
      
      # 3. 构造 W
      W <- matrix(0, n*K, n*K)
      u <- rep(0, n*K)
      eps_W <- 1e-8
      for (i in 1:n) {
        li <- lambda_no0[i, ]
        Wi <- diag(li) - outer(li, li)
        
        gamma_i <- gat[i, -1]
        eta_i <- eta_mat[i, ]
        
        delta <- tryCatch(solve(Wi, gamma_i - li),
                          error = function(e) solve(Wi + diag(eps_W, K), gamma_i - li))
        ui <- eta_i + delta
        
        row_start <- (i-1)*K + 1
        row_end <- i*K
        W[row_start:row_end, row_start:row_end] <- Wi
        u[row_start:row_end] <- ui
      }
      
      # 4. IRLS 更新
      lhs <- crossprod(Ztilde, W %*% Ztilde)
      rhs <- crossprod(Ztilde, W %*% u)
      
      beta_vec_new <- solve(lhs, rhs)
      beta_new <- matrix(beta_vec_new, nrow=p, ncol=K)
      # cat("Beta:", beta_new, "\n")
      Q.val.new <- Q.fun(x, z, beta_new, gat, m)
      # cat("Q:", Q.val.new, "\n")
      
      # 5. 收敛判定
      diff <- abs((Q.val.new - Q.val) / Q.val)
      if (diff < ep) {
        message("IRLS converged in ", iter, " iterations")
        return(beta_new)
      }
      beta <- beta_new
      Q.val <- Q.val.new
    }
    
    warning("IRLS did not converge in maxit iterations")
    return(beta)
  }
  
  n <- length(x)
  p <- ncol(z)
  K <- ncol(beta)
  beta_curr <- beta
  
  # 初始对数似然
  ll_curr <- loglikeli.fun(x, z, beta_curr, m)
  loglik_trace <- numeric(0)
  loglik_trace <- c(loglik_trace, ll_curr)
  
  for (iter in 1:maxit) {
    # 1) 计算当前 gamma
    gat <- gamma.fun(x, z, beta_curr, m)   # n x (K+1)
    
    # 2) 内层 IRLS（最大化 Q）
    beta_new <- IRLS(gat = gat, z = z, beta_ini = beta_curr)
    
    # 3) 计算新的对数似然
    ll_new <- loglikeli.fun(x, z, beta_new, m)
    loglik_trace <- c(loglik_trace, ll_new)
    
    # 4) 收敛判定（对数似然变化率，稳健形式）
    diff_ll <- abs((ll_new - ll_curr) / ll_curr)
    if (diff_ll < ep) {
      return(list(beta = beta_new, iterations = iter, converged = TRUE,
                  loglik = ll_new, loglik_trace = loglik_trace))
    }
    
    # 5) 更新
    beta_curr <- beta_new
    ll_curr <- ll_new
  }
  
  warning("Regression: 外层迭代达到最大次数，但未满足收敛条件")
  return(list(beta = beta_curr, iterations = maxit, converged = FALSE,
              loglik = ll_curr, loglik_trace = loglik_trace))
}

Regression(x, z, m, beta)
