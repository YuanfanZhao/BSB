# ============================================================================
# 00_regression_lib.R —— 二项响应回归模型函数库（第 2.5 节模型 + 对比模型）
# ----------------------------------------------------------------------------
# 数据形态：对每个个体 i，观测 (m_i, x_i, z_i)：
#   - x_i = 成功次数（0..m_i），m_i = 试验次数（伯努利重复试验）；
#   - z_i = 协变量向量（第 1 个元素为 1，表示截距）。
# 模型：
#   1) 逻辑回归（Logistic）：X_i ~ Binomial(m_i, p_i)，logit(p_i) = z_i^T alpha
#   2) Beta-Binomial 回归：X_i ~ BetaBin(m_i, a_i, b_i)，均值 logit 链接 + 离散参数 rho
#   3) Logit-Normal-Binomial 回归：X_i | p_i ~ Binomial(m_i, p_i)，
#      logit(p_i) ~ N(z_i^T alpha, sigma^2)
#   4) Bern-Bino 回归（论文第 2.5 节）：p_i ~ Bernstein(K, lambda_i)，
#      lambda_ik = exp(eta_ik)/sum_l exp(eta_il)，eta_ik = phi_k + gamma_k * z_i^T beta，
#      约束 phi_0 = 0，gamma_0 = 1（保证可识别）。
# 本文件只定义函数，不执行分析；由 02 脚本调用。
# ============================================================================

# ----------------------------------------------------------------------------
# 0. 工具函数
# ----------------------------------------------------------------------------
lchoose_mx <- function(m, x) lgamma(m + 1) - lgamma(x + 1) - lgamma(m - x + 1)

# Logit-Normal 积分的核心：∫ plogis(u)^(x+ppow) (1-plogis(u))^(m-x) * dnorm(u, mu, sigma) du
int_ln <- function(x, m, mu, sigma, ppow = 0) {
  f <- function(u) plogis(u)^(x + ppow) * (1 - plogis(u))^(m - x) * dnorm(u, mu, sigma)
  v <- tryCatch(integrate(f, -20, 20, rel.tol = 1e-8, subdivisions = 400,
                          stop.on.error = FALSE)$value, error = function(e) NA)
  if (!is.finite(v) || v <= 0) 1e-300 else v
}

# Bern-Bino 的系数 a_k(x) = B(k+1+x, K-k+1+m-x) / B(k+1, K-k+1)
bern_coef <- function(m, x, K) {
  sapply(0:K, function(k)
    exp(lbeta(k + 1 + x, K - k + 1 + m - x) - lbeta(k + 1, K - k + 1)))
}

# ----------------------------------------------------------------------------
# 1. 逻辑回归（Logistic）
#    用 glm 拟合，返回可计算测试集 logLik 的预测概率函数。
# ----------------------------------------------------------------------------
fit_logistic <- function(X_tr, m_tr, x_tr) {
  fit <- glm(cbind(x_tr, m_tr - x_tr) ~ X_tr[, -1, drop = FALSE], family = binomial)
  # 预测概率（对任意 X 矩阵）
  pred <- function(X) {
    eta <- cbind(1, X[, -1, drop = FALSE]) %*% coef(fit)
    plogis(eta)
  }
  list(npar = length(coef(fit)), pred = pred)
}

# ----------------------------------------------------------------------------
# 2. Beta-Binomial 回归
#    参数：alpha（logit 均值的系数）+ rho（组内相关/离散参数，logit 链接保证在 (0,1)）
#    a_i = mu_i (1-rho)/rho，b_i = (1-mu_i)(1-rho)/rho
# ----------------------------------------------------------------------------
bbreg_ll <- function(theta, X, m, x) {
  p <- ncol(X)
  alpha <- theta[1:p]
  rho <- plogis(theta[p + 1])
  mu <- plogis(X %*% alpha)
  a <- mu * (1 - rho) / rho
  b <- (1 - mu) * (1 - rho) / rho
  sum(lchoose_mx(m, x) + lbeta(a + x, b + m - x) - lbeta(a, b))
}
fit_bbreg <- function(X_tr, m_tr, x_tr, start = NULL) {
  p <- ncol(X_tr)
  # 初始值：用逻辑回归系数作 alpha，rho 初始 0.2
  glm0 <- glm(cbind(x_tr, m_tr - x_tr) ~ X_tr[, -1, drop = FALSE], family = binomial)
  if (is.null(start)) start <- c(coef(glm0), qlogis(0.2))
  opt <- optim(start, function(th) -bbreg_ll(th, X_tr, m_tr, x_tr),
               method = "BFGS", control = list(maxit = 2000))
  theta <- opt$par
  list(npar = length(theta), theta = theta, ll = -opt$value,
       pred = function(X) {
         p <- ncol(X)
         mu <- plogis(X %*% theta[1:p])
         rho <- plogis(theta[p + 1])
         a <- mu * (1 - rho) / rho
         b <- (1 - mu) * (1 - rho) / rho
         # 返回每个观测的 pmf 向量函数，便于计算测试 logLik
         list(mu = mu, a = a, b = b)
       })
}
bbreg_pmf <- function(a, b, m, x) exp(lchoose_mx(m, x) + lbeta(a + x, b + m - x) - lbeta(a, b))

# ----------------------------------------------------------------------------
# 3. Logit-Normal-Binomial 回归
#    参数：alpha（logit 均值的系数）+ log(sigma)
# ----------------------------------------------------------------------------
lnreg_ll <- function(theta, X, m, x) {
  p <- ncol(X)
  alpha <- theta[1:p]
  sigma <- exp(theta[p + 1])
  mu <- plogis(X %*% alpha)
  sum(log(sapply(seq_along(x), function(i)
    exp(lchoose_mx(m[i], x[i])) * int_ln(x[i], m[i], logit(mu[i]), sigma))))
}
fit_lnreg <- function(X_tr, m_tr, x_tr, start = NULL) {
  p <- ncol(X_tr)
  glm0 <- glm(cbind(x_tr, m_tr - x_tr) ~ X_tr[, -1, drop = FALSE], family = binomial)
  if (is.null(start)) start <- c(coef(glm0), log(0.5))
  opt <- optim(start, function(th) -lnreg_ll(th, X_tr, m_tr, x_tr),
               method = "BFGS", control = list(maxit = 2000))
  theta <- opt$par
  list(npar = length(theta), theta = theta, ll = -opt$value,
       pred = function(X) {
         p <- ncol(X)
         list(mu = plogis(X %*% theta[1:p]), sigma = exp(theta[p + 1]))
       })
}
lnreg_pmf <- function(mu, sigma, m, x) exp(lchoose_mx(m, x)) * int_ln(x, m, logit(mu), sigma)

# ----------------------------------------------------------------------------
# 4. Bern-Bino 回归（论文第 2.5 节）
#    参数 theta = (phi_1..phi_K, gamma_1..gamma_K, beta_1..beta_p)
#    约束 phi_0 = 0，gamma_0 = 1。
# ----------------------------------------------------------------------------
bernreg_ll <- function(theta, X, m, x, K) {
  p <- ncol(X)
  phi <- c(0, theta[1:K])
  gam <- c(1, theta[(K + 1):(2 * K)])
  beta <- theta[(2 * K + 1):(2 * K + p)]
  Zb <- as.numeric(X %*% beta)               # n 维
  n <- length(m)
  # eta_ik = phi_k + gamma_k * Zb_i
  eta <- outer(Zb, gam) + matrix(rep(phi, each = n), n, K + 1)
  maxe <- apply(eta, 1, max)
  lam <- exp(eta - maxe)
  lam <- lam / rowSums(lam)
  # a_ik
  A <- t(sapply(seq_len(n), function(i) bern_coef(m[i], x[i], K)))
  val <- rowSums(A * lam)
  sum(log(val)) + sum(lchoose_mx(m, x))
}

# Bern-Bino 回归的对数似然及其解析梯度（供 BFGS 使用，收敛更快更稳）
# 记 w_ik = a_ik * lambda_ik / den_i（观测 i 属于成分 k 的后验概率），
# 则 d ll / d eta_ik = w_ik - lambda_ik，由此链式法则得到各参数梯度。
bernreg_ll_grad <- function(theta, X, m, x, K) {
  p <- ncol(X)
  phi <- c(0, theta[1:K]); gam <- c(1, theta[(K + 1):(2 * K)])
  beta <- theta[(2 * K + 1):(2 * K + p)]
  Zb <- as.numeric(X %*% beta)
  n <- length(m)
  eta <- outer(Zb, gam) + matrix(rep(phi, each = n), n, K + 1)
  maxe <- apply(eta, 1, max)
  lam <- exp(eta - maxe); lam <- lam / rowSums(lam)
  A <- t(sapply(seq_len(n), function(i) bern_coef(m[i], x[i], K)))
  den <- rowSums(A * lam)                    # 每行 = sum_k a_ik lambda_ik（n 维）
  w <- A * lam / den                          # 后验成分权重 w_ik
  diff <- w - lam                             # d ll / d eta_ik
  ll <- sum(log(den)) + sum(lchoose_mx(m, x))
  # 梯度
  g_phi <- colSums(diff[, 2:(K + 1), drop = FALSE])           # j=1..K
  g_gam <- colSums(diff[, 2:(K + 1), drop = FALSE] * Zb)      # j=1..K
  g_beta <- colSums(as.numeric(diff %*% gam) * X)              # p 维
  list(ll = ll, grad = c(g_phi, g_gam, g_beta))
}
fit_bernreg <- function(X_tr, m_tr, x_tr, K, start = NULL, maxit = 800) {
  p <- ncol(X_tr)
  fn <- function(th) -bernreg_ll_grad(th, X_tr, m_tr, x_tr, K)$ll
  gr <- function(th) -bernreg_ll_grad(th, X_tr, m_tr, x_tr, K)$grad
  if (is.null(start)) {
    glm0 <- glm(cbind(x_tr, m_tr - x_tr) ~ X_tr[, -1, drop = FALSE], family = binomial)
    start <- c(rep(0, K), rep(1, K), coef(glm0))
  }
  opt <- optim(start, fn, gr = gr, method = "BFGS",
               control = list(maxit = maxit, reltol = 1e-10))
  if (opt$convergence != 0) {               # 未收敛则扰动重启一次，取更优
    opt2 <- optim(start + c(rnorm(2 * K, 0, 0.5), rnorm(p, 0, 0.1)), fn, gr = gr,
                  method = "BFGS", control = list(maxit = maxit, reltol = 1e-10))
    if (-opt2$value > -opt$value) opt <- opt2
  }
  theta <- opt$par
  list(npar = 2 * K + p, K = K, theta = theta, ll = -opt$value, conv = opt$convergence,
       pred = function(X) {
         p <- ncol(X)
         phi <- c(0, theta[1:K]); gam <- c(1, theta[(K + 1):(2 * K)])
         beta <- theta[(2 * K + 1):(2 * K + p)]
         Zb <- as.numeric(X %*% beta)
         n <- nrow(X)
         eta <- outer(Zb, gam) + matrix(rep(phi, each = n), n, K + 1)
         maxe <- apply(eta, 1, max)
         lam <- exp(eta - maxe); lam <- lam / rowSums(lam)
         list(lambda = lam, phi = phi, gam = gam, beta = beta)
       })
}

bernreg_pmf <- function(lambda, m, x, K) {
  A <- bern_coef(m, x, K)
  exp(lchoose_mx(m, x)) * sum(lambda * A)
}

# ----------------------------------------------------------------------------
# 5. K 选择（在训练集上，用 BIC）：K = 1..Kmax
# ----------------------------------------------------------------------------
bernreg_fit_bic <- function(X_tr, m_tr, x_tr, Kmax = 8) {
  n1 <- length(m_tr)
  best <- NULL
  start <- NULL
  for (K in 1:Kmax) {
    fit <- fit_bernreg(X_tr, m_tr, x_tr, K, start = start)   # warm start
    bic <- fit$npar * log(n1) - 2 * fit$ll
    if (is.null(best) || bic < best$bic) best <- list(K = K, fit = fit, bic = bic)
    # 下一个 K 的 warm start：新类别参数补 0/1，其余沿用
    start <- c(fit$theta[1:(K)], 0, fit$theta[(K + 1):(2 * K)], 1, fit$theta[(2 * K + 1):length(fit$theta)])
  }
  best
}

# 训练集内部 5 折 x 3 次交叉验证选 K（验证准则为对数似然）
bernreg_fit_cv <- function(X_tr, m_tr, x_tr, Kmax = 8, folds = 5, reps = 3, seed = 6789) {
  n1 <- length(m_tr)
  set.seed(seed)
  mean_val <- numeric(Kmax)
  for (K in 1:Kmax) {
    v <- c()
    start <- NULL
    for (r in 1:reps) {
      foldid <- sample(rep(1:folds, length.out = n1))
      for (f in 1:folds) {
        va <- which(foldid == f); tr <- which(foldid != f)
        fit <- fit_bernreg(X_tr[tr, , drop = FALSE], m_tr[tr], x_tr[tr], K, start = start)
        # 验证折对数似然
        v <- c(v, test_loglik('bernreg', fit, X_tr[va, , drop = FALSE], m_tr[va], x_tr[va]))
        start <- fit$theta     # 同一 K 的下一折：直接用上折的解作初值（长度一致）
      }
    }
    mean_val[K] <- mean(v)
  }
  K <- which.max(mean_val)
  list(K = K, fit = fit_bernreg(X_tr, m_tr, x_tr, K), mean_val = mean_val)
}

# ----------------------------------------------------------------------------
# 6. 测试集评价
#    test logLik = sum_i log P(x_i | 训练拟合模型, z_i)
#    test MAE    = mean | x_i/m_i - 模型预测均值比例 |
# ----------------------------------------------------------------------------
test_loglik <- function(model, fit, X_te, m_te, x_te) {
  switch(model,
    logistic = {
      p <- fit$pred(X_te)
      sum(dbinom(x_te, m_te, p, log = TRUE))
    },
    bbreg = {
      pr <- fit$pred(X_te)
      sum(log(bbreg_pmf(pr$a, pr$b, m_te, x_te)))
    },
    lnreg = {
      pr <- fit$pred(X_te)
      sum(log(sapply(seq_along(x_te), function(i)
        lnreg_pmf(pr$mu[i], pr$sigma, m_te[i], x_te[i]))))
    },
    bernreg = {
      K <- fit$K
      pr <- fit$pred(X_te)
      sum(log(sapply(seq_along(x_te), function(i)
        bernreg_pmf(pr$lambda[i, ], m_te[i], x_te[i], K))))
    }
  )
}

test_mae <- function(model, fit, X_te, m_te, x_te) {
  switch(model,
    logistic = mean(abs(x_te / m_te - fit$pred(X_te))),
    bbreg = { pr <- fit$pred(X_te); mean(abs(x_te / m_te - pr$mu)) },
    lnreg = { pr <- fit$pred(X_te); mean(abs(x_te / m_te - pr$mu)) },
    bernreg = {
      K <- fit$K; pr <- fit$pred(X_te)
      # 模型预测均值 E(X|z) = m_i * sum_k lambda_ik (k+1)/(K+2)
      mu1 <- as.numeric(pr$lambda %*% ((0:K + 1) / (K + 2)))
      mean(abs(x_te / m_te - mu1))
    }
  )
}
