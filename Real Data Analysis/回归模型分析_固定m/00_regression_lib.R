# ============================================================================
# 00_regression_lib.R --- 二项响应回归模型函数库（固定 m 版本，第 2.5 节模型 + 对比模型）
# ----------------------------------------------------------------------------
# 数据形态：对每个个体 i，观测 (x_i, m, z_i)，其中 m 对所有个体相同（固定），
#   x_i = 成功次数（取值 0..m），z_i = 协变量向量（第 1 个元素为 1，表示截距）。
#
# 模型（用于 70/30 训练/测试比较）：
#   1) 逻辑回归（Logistic）：X_i ~ Binomial(m, p_i)，logit(p_i) = z_i^T alpha
#   2) Beta-Binomial 回归：  X_i ~ BetaBin(m, a_i, b_i)，均值 logit 链接 + 离散参数 rho
#   3) Bern-Bino 回归（论文第 2.5 节）：X_i | p_i ~ Binomial(m, p_i)，
#      p_i ~ Bernstein(K, lambda_i)，lambda_ik = softmax(eta_ik)，
#      eta_ik = phi_k + gamma_k * z_i^T beta，约束 phi_0 = 0、gamma_0 = 1 保证可识别。
#
# 固定 m 的优化：a_k(x) = B(k+1+x, K-k+1+m-x) / B(k+1, K-k+1) 只依赖 x，
#   因此预先计算 (K+1) x (m+1) 的表 A，之后所有个体直接查表，避免逐观测重复计算。
# 本文件只定义函数，不执行分析；由 02 脚本调用。
# ============================================================================

# ----------------------------------------------------------------------------
# 0. 工具函数
# ----------------------------------------------------------------------------
lchoose_mx <- function(m, x) lgamma(m + 1) - lgamma(x + 1) - lgamma(m - x + 1)

# 预先计算固定 m 下 Bern-Bino 的系数表 A[k, x]（k = 0..K, x = 0..m）
#   A[k, x] = B(k+1+x, K-k+1+m-x) / B(k+1, K-k+1)
bern_coef_table <- function(m, K) {
  ks <- 0:K
  sapply(0:m, function(x)
    exp(lbeta(ks + 1 + x, K - ks + 1 + m - x) - lbeta(ks + 1, K - ks + 1)))  # (K+1) x (m+1)
}

# ----------------------------------------------------------------------------
# 1. 逻辑回归（Logistic）
#    用 glm 拟合，返回可计算测试集 logLik 的预测概率函数。
# ----------------------------------------------------------------------------
fit_logistic <- function(X_tr, m, x_tr) {
  fit <- glm(cbind(x_tr, m - x_tr) ~ X_tr[, -1, drop = FALSE], family = binomial)
  pred <- function(X) {
    eta <- cbind(1, X[, -1, drop = FALSE]) %*% coef(fit)
    plogis(eta)
  }
  list(npar = length(coef(fit)), pred = pred)
}

# ----------------------------------------------------------------------------
# 2. Beta-Binomial 回归
#    参数：alpha（logit 均值的系数）+ rho（离散参数，logit 链接保证在 (0,1)）
#    a_i = mu_i (1-rho)/rho，b_i = (1-mu_i)(1-rho)/rho
# ----------------------------------------------------------------------------
# 数值稳定的 Beta-Binomial 对数概率（上升阶乘形式，避免 lbeta 大参数相消误差）：
#   log[ B(a+x, b+m-x)/B(a,b) ] = sum_{j=0}^{x-1} log(a+j) + sum_{j=0}^{m-x-1} log(b+j)
#                                 - sum_{j=0}^{m-1} log(a+b+j)
# 当 rho -> 0（a,b 巨大）时自动退化为二项概率，数值稳健。
bbreg_ll <- function(theta, X, m, x) {
  p <- ncol(X)
  alpha <- theta[1:p]
  rho <- plogis(theta[p + 1])
  mu <- plogis(X %*% alpha)
  a <- mu * (1 - rho) / rho
  b <- (1 - mu) * (1 - rho) / rho
  # 逐观测稳定计算 log pmf
  n <- length(x)
  ll <- numeric(n)
  for (i in seq_len(n)) {
    xi <- x[i]
    s <- lchoose_mx(m, xi)
    if (xi > 0)          s <- s + sum(log(a[i] + 0:(xi - 1)))
    if (m - xi > 0)      s <- s + sum(log(b[i] + 0:(m - xi - 1)))
    if (m > 0)           s <- s - sum(log(a[i] + b[i] + 0:(m - 1)))
    ll[i] <- s
  }
  sum(ll)
}
fit_bbreg <- function(X_tr, m, x_tr, start = NULL) {
  p <- ncol(X_tr)
  glm0 <- glm(cbind(x_tr, m - x_tr) ~ X_tr[, -1, drop = FALSE], family = binomial)
  if (is.null(start)) start <- c(coef(glm0), qlogis(0.2))
  opt <- optim(start, function(th) -bbreg_ll(th, X_tr, m, x_tr),
               method = "BFGS", control = list(maxit = 2000))
  theta <- opt$par
  list(npar = length(theta), theta = theta, ll = -opt$value,
       pred = function(X) {
         p <- ncol(X)
         mu <- plogis(X %*% theta[1:p])
         rho <- plogis(theta[p + 1])
         a <- mu * (1 - rho) / rho
         b <- (1 - mu) * (1 - rho) / rho
         list(mu = mu, a = a, b = b)
       })
}
bbreg_pmf <- function(a, b, m, x) {
  n <- length(x)
  out <- numeric(n)
  for (i in seq_len(n)) {
    xi <- x[i]
    s <- lchoose_mx(m, xi)
    if (xi > 0)          s <- s + sum(log(a[i] + 0:(xi - 1)))
    if (m - xi > 0)      s <- s + sum(log(b[i] + 0:(m - xi - 1)))
    if (m > 0)           s <- s - sum(log(a[i] + b[i] + 0:(m - 1)))
    out[i] <- exp(s)
  }
  out
}

# ----------------------------------------------------------------------------
# 3. Bern-Bino 回归（论文第 2.5 节，固定 m 版本）
#    对数似然与解析梯度（已与数值梯度核对，误差 ~1e-8，BFGS 稳定收敛）
# ----------------------------------------------------------------------------
bernreg_ll_grad <- function(theta, X, m, x, A_tab) {
  K <- nrow(A_tab) - 1
  n <- length(x); p <- ncol(X)
  phi <- c(0, theta[1:K]); gam <- c(1, theta[(K + 1):(2 * K)])
  beta <- theta[(2 * K + 1):(2 * K + p)]
  Zb <- as.numeric(X %*% beta)
  eta <- outer(Zb, gam) + matrix(rep(phi, each = n), n, K + 1)
  maxe <- apply(eta, 1, max)
  lam <- exp(eta - maxe); lam <- lam / rowSums(lam)         # softmax，n x (K+1)
  A <- A_tab[, x + 1, drop = FALSE]                          # 查表：(K+1) x n，A[k,i] = a_k(x_i)
  P <- t(lam) * A                                            # (K+1) x n，P[k,i] = lam_ik * a_k(x_i)
  dens <- colSums(P)                                         # n：sum_k lam_ik a_k(x_i)
  diff <- t(P) / dens - lam                                  # n x (K+1)：d logL_i / d eta_ik
  ll <- sum(lchoose_mx(m, x) + log(dens))
  g_phi  <- colSums(diff[, 2:(K + 1), drop = FALSE])         # j = 1..K
  g_gam  <- colSums(diff[, 2:(K + 1), drop = FALSE] * Zb)    # j = 1..K
  g_beta <- colSums(as.numeric(diff %*% gam) * X)            # p 维
  list(ll = ll, grad = c(g_phi, g_gam, g_beta))
}

fit_bernreg <- function(X_tr, m, x_tr, K, start = NULL, maxit = 800, nstart = 8,
                         seed = 123) {
  p <- ncol(X_tr)
  A_tab <- bern_coef_table(m, K)
  fn <- function(th) -bernreg_ll_grad(th, X_tr, m, x_tr, A_tab)$ll
  gr <- function(th) -bernreg_ll_grad(th, X_tr, m, x_tr, A_tab)$grad
  if (is.null(start)) {
    glm0 <- glm(cbind(x_tr, m - x_tr) ~ X_tr[, -1, drop = FALSE], family = binomial)
    start <- c(rep(0, K), rep(1, K), coef(glm0))
  }
  # 多起点（multi-start）：基础起点（逻辑回归系数构造）恒在首位，再补若干随机扰动
  # 起点（尺度 0.6 / 1.5 / 2.5 交替），取最优者，避免局部最优；固定种子保证可复现。
  # 注意：不使用"上一 K 的解"作唯一锚点，因为劣质局部解会沿 K 链传播。
  starts <- list(start)
  if (nstart > 1) {
    set.seed(2024)
    scales <- c(0.6, 1.5, 2.5)
    for (s in 2:nstart) {
      sc <- scales[(s - 2) %% 3 + 1]
      starts[[s]] <- start + c(rnorm(2 * K, 0, sc), rnorm(p, 0, 0.3 * sc))
    }
  }
  best <- NULL
  for (st in starts) {
    opt <- optim(st, fn, gr = gr, method = "BFGS",
                 control = list(maxit = maxit, reltol = 1e-10))
    if (opt$convergence != 0) {            # 未收敛则扰动重启一次，取更优
      opt2 <- optim(st + c(rnorm(2 * K, 0, 0.5), rnorm(p, 0, 0.1)), fn, gr = gr,
                    method = "BFGS", control = list(maxit = maxit, reltol = 1e-10))
      if (-opt2$value > -opt$value) opt <- opt2
    }
    if (is.null(best) || -opt$value > best$ll) {
      best <- list(theta = opt$par, ll = -opt$value, conv = opt$convergence)
    }
  }
  theta <- best$theta
  list(npar = 2 * K + p, K = K, theta = theta, ll = best$ll, conv = best$conv,
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

bernreg_pmf <- function(lambda, m, x, A_tab) {
  exp(lchoose_mx(m, x)) * as.numeric(lambda %*% A_tab[, x + 1])
}

# ----------------------------------------------------------------------------
# 4. K 选择（只在训练集上进行）
#    (a) 训练集 BIC：K = 1..Kmax，取 BIC 最小者
#    (b) 训练集内部 5 折 x 3 次交叉验证：验证准则为对数似然，取均值最大者
# ----------------------------------------------------------------------------
bernreg_fit_bic <- function(X_tr, m, x_tr, Kmax = 10, nstart = 12) {
  n1 <- length(x_tr)
  best <- NULL
  curve <- data.frame(K = 1:Kmax, train_ll = NA, npar = NA, BIC = NA)
  for (K in 1:Kmax) {
    # 每个 K 都用逻辑回归系数构造的基础起点（fit_bernreg 内部会多起点），
    # 不使用"上一 K 的解"作锚点，避免劣质局部解沿 K 链传播。
    fit <- fit_bernreg(X_tr, m, x_tr, K, start = NULL, nstart = nstart)
    bic <- fit$npar * log(n1) - 2 * fit$ll
    curve$train_ll[K] <- fit$ll; curve$npar[K] <- fit$npar; curve$BIC[K] <- bic
    if (is.null(best) || bic < best$bic) best <- list(K = K, fit = fit, bic = bic)
  }
  best$curve <- curve
  best
}

bernreg_fit_cv <- function(X_tr, m, x_tr, Kmax = 10, folds = 5, reps = 1,
                         seed = 6789, nstart = 5) {
  n1 <- length(x_tr)
  set.seed(seed)
  mean_val <- numeric(Kmax)
  for (K in 1:Kmax) {
    v <- c(); start <- NULL
    for (r in 1:reps) {
      foldid <- sample(rep(1:folds, length.out = n1))
      for (f in 1:folds) {
        va <- which(foldid == f); tr <- which(foldid != f)
        fit <- fit_bernreg(X_tr[tr, , drop = FALSE], m, x_tr[tr], K,
                           start = start, nstart = nstart)
        v <- c(v, test_loglik("bernreg", fit, X_tr[va, , drop = FALSE], m, x_tr[va], K))
        start <- fit$theta
      }
    }
    mean_val[K] <- mean(v)
  }
  K <- which.max(mean_val)
  list(K = K, fit = fit_bernreg(X_tr, m, x_tr, K, nstart = 12), mean_val = mean_val)
}

# ----------------------------------------------------------------------------
# 5. 测试集评价（测试集完全不参与拟合与 K 选择）
#    test logLik = sum_i log P(x_i | 训练拟合模型, z_i)
#    test MAE    = mean | x_i/m - 模型预测均值比例 |
#    test MSE    = mean (x_i/m - 模型预测均值比例)^2
# ----------------------------------------------------------------------------
test_loglik <- function(model, fit, X_te, m, x_te, K = NULL) {
  switch(model,
    logistic = {
      p <- fit$pred(X_te)
      sum(dbinom(x_te, m, p, log = TRUE))
    },
    bbreg = {
      pr <- fit$pred(X_te)
      sum(log(bbreg_pmf(pr$a, pr$b, m, x_te)))
    },
    bernreg = {
      if (is.null(K)) K <- fit$K
      A_tab <- bern_coef_table(m, K)
      pr <- fit$pred(X_te)
      sum(log(sapply(seq_along(x_te), function(i)
        bernreg_pmf(pr$lambda[i, ], m, x_te[i], A_tab))))
    }
  )
}

test_mae <- function(model, fit, X_te, m, x_te, K = NULL) {
  switch(model,
    logistic = mean(abs(x_te / m - fit$pred(X_te))),
    bbreg = { pr <- fit$pred(X_te); mean(abs(x_te / m - pr$mu)) },
    bernreg = {
      if (is.null(K)) K <- fit$K
      pr <- fit$pred(X_te)
      mu1 <- as.numeric(pr$lambda %*% ((0:K + 1) / (K + 2)))   # E(p_i | z_i)
      mean(abs(x_te / m - mu1))
    }
  )
}

test_mse <- function(model, fit, X_te, m, x_te, K = NULL) {
  switch(model,
    logistic = mean((x_te / m - fit$pred(X_te))^2),
    bbreg = { pr <- fit$pred(X_te); mean((x_te / m - pr$mu)^2) },
    bernreg = {
      if (is.null(K)) K <- fit$K
      pr <- fit$pred(X_te)
      mu1 <- as.numeric(pr$lambda %*% ((0:K + 1) / (K + 2)))
      mean((x_te / m - mu1)^2)
    }
  )
}
