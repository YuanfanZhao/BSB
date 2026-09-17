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
# 3. Bern-Bino 回归（论文第 2.5 节，固定 m 版本，已修正模型设定）
# ----------------------------------------------------------------------------
# 模型：
#   X_i | p_i ~ Binomial(m, p_i),   p_i ~ Bernstein(K, lambda_i),
#   lambda_ik = exp(eta_ik) / sum_l exp(eta_il),
#   eta_ik    = phi_k + gamma_k * ( z_i^T beta ),      k = 0, 1, ..., K,
# 其中
#   * z_i 为【不含截距】的协变量向量（长度 p'）；
#   * phi_k 为第 k 个 Bernstein 成分的【专属截距】，phi_0 = 0 作为基线；
#   * beta  为所有成分【共享】的回归系数，只作用于非截距协变量；
#   * gamma_k 为第 k 个成分在共享线性指标上的【缩放尺度】，gamma_0 = 1 作为基线。
# 参数个数 = 2K + p'（phi_0 与 gamma_0 固定，不计入）。
#
# 【关键】截距项必须全部由 phi_k 承担，不能包含在 z_i 中；否则
#   phi_k + gamma_k * beta_0 中的 phi_k 与 gamma_k*beta_0 互为冗余，
#   会引入额外的不可识别方向。因此本函数要求 Z 为“纯协变量”矩阵。
# ----------------------------------------------------------------------------
bernreg_ll_grad <- function(theta, Z, m, x, A_tab) {
  K <- nrow(A_tab) - 1
  n <- length(x); p <- ncol(Z)
  phi <- c(0, theta[1:K]); gam <- c(1, theta[(K + 1):(2 * K)])
  beta <- theta[(2 * K + 1):(2 * K + p)]
  Zb <- as.numeric(Z %*% beta)                               # 共享线性指标（无截距）
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
  g_beta <- colSums(as.numeric(diff %*% gam) * Z)            # p 维
  list(ll = ll, grad = c(g_phi, g_gam, g_beta))
}

# 边缘（不含协变量）Bern-Bino 的 MM 拟合，用于为 phi_k 提供初始值：
#   pmf(x) = C(m,x) * sum_k lambda_k a_k(x),  lambda_k^{new} = sum_x c_x gamma_{kx}/n
bern_marginal_mm <- function(m, x, K, maxit = 100000, eps = 1e-12) {
  counts <- tabulate(x + 1, nbins = m + 1); n <- sum(counts)
  A <- bern_coef_table(m, K)                  # (K+1) x (m+1)
  lc <- lchoose_mx(m, 0:m)
  lam <- rep(1 / (K + 1), K + 1)
  ll <- sum(counts * (lc + log(as.numeric(lam %*% A))))
  for (it in 1:maxit) {
    den <- as.numeric(lam %*% A)
    prop <- (lam * A) / rep(den, each = K + 1)
    lam_new <- as.numeric(prop %*% counts) / n
    ll_new <- sum(counts * (lc + log(as.numeric(lam_new %*% A))))
    if (abs(ll_new - ll) < eps * (abs(ll) + 1e-12)) { lam <- lam_new; ll <- ll_new; break }
    lam <- lam_new; ll <- ll_new
  }
  list(lambda = lam, ll = ll, iter = it)
}

# Z_tr：不含截距的协变量矩阵（n1 x p'）
fit_bernreg <- function(Z_tr, m, x_tr, K, start = NULL, maxit = 800, nstart = 12,
                        nindep = 12, seed = 123) {
  p <- ncol(Z_tr)
  A_tab <- bern_coef_table(m, K)
  fn <- function(th) -bernreg_ll_grad(th, Z_tr, m, x_tr, A_tab)$ll
  gr <- function(th) -bernreg_ll_grad(th, Z_tr, m, x_tr, A_tab)$grad
  extra_starts <- list()
  glm0 <- glm(cbind(x_tr, m - x_tr) ~ Z_tr, family = binomial)
  slopes <- coef(glm0)[-1]
  if (is.null(start)) {
    # 共享回归系数取逻辑回归的斜率（截距不进入 beta，全部由 phi_k 承担）。
    # 起点 1（主）：用【边缘 Bern-Bino 的 MM 解】给出 phi_k，gamma_k = 1。
    mg <- bern_marginal_mm(m, x_tr, K)
    ph <- log(pmax(mg$lambda, 1e-12)); ph <- ph - ph[1]      # 使 phi_0 = 0
    start <- c(ph[-1], rep(1, K), slopes)
    # 起点 2、3：各成分同处基线（phi_k = 0）或取逻辑回归截距水平。
    extra_starts <- list(c(rep(0, K), rep(1, K), slopes),
                         c(rep(coef(glm0)[1], K), rep(1, K), slopes))
  }
  # 多起点（multi-start）：基础起点恒在首位，再补若干随机扰动起点
  # （尺度 0.6 / 1.5 / 2.5 交替），取最优者；固定种子保证可复现。
  anchors <- c(list(start), extra_starts)
  starts <- anchors
  # (i) 锚定扰动起点：围绕三个确定性起点做不同尺度的扰动
  set.seed(seed)
  scales <- c(0.6, 1.5, 2.5, 4)
  while (length(starts) < nstart) {
    i <- length(starts) - length(anchors)
    sc <- scales[i %% length(scales) + 1]
    a  <- anchors[[i %% length(anchors) + 1]]
    starts[[length(starts) + 1]] <- a + c(rnorm(2 * K, 0, sc), rnorm(p, 0, 0.3 * sc))
  }
  # (ii) 独立随机起点：不锚定在边缘解上，用于探索“所有成分都被激活”的解。
  # 似然曲面高度多峰，独立随机起点是找到全局最优的关键，nindep 可调大。
  set.seed(seed + 1)
  for (s in seq_len(nindep)) starts[[length(starts) + 1]] <-
    c(rnorm(K, 0, 2), 1 + rnorm(K, 0, 3), slopes + rnorm(p, 0, 0.5))
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
  list(npar = 2 * K + p, pz = p, K = K, theta = theta, ll = best$ll, conv = best$conv,
       pred = function(Z) {
         p <- ncol(Z)
         phi <- c(0, theta[1:K]); gam <- c(1, theta[(K + 1):(2 * K)])
         beta <- theta[(2 * K + 1):(2 * K + p)]
         Zb <- as.numeric(Z %*% beta)
         n <- nrow(Z)
         eta <- outer(Zb, gam) + matrix(rep(phi, each = n), n, K + 1)
         maxe <- apply(eta, 1, max)
         lam <- exp(eta - maxe); lam <- lam / rowSums(lam)
         list(lambda = lam, phi = phi, gam = gam, beta = beta, index = Zb)
       })
}

# Z：不含截距的协变量矩阵；lambda 为 n x (K+1) 权重矩阵
bernreg_pmf <- function(lambda, m, x, A_tab) {
  exp(lchoose_mx(m, x)) * as.numeric(lambda %*% A_tab[, x + 1])
}


# ----------------------------------------------------------------------------
# 3b. Bern-Bino 回归（变体 A：每个成分各自的回归系数）
# ----------------------------------------------------------------------------
# 模型：
#   lambda_ik = exp(eta_ik) / sum_l exp(eta_il),
#   eta_ik    = z_i^T beta_k,          k = 0, 1, ..., K,
# 其中 z_i 为含截距的协变量向量（长度 p）。
# 可识别性：softmax 对"同一条观测的所有 eta_k 同时加上 z_i^T delta"不变，
#   故 beta_0, ..., beta_K 之间存在 p 维不可识别方向。取 beta_0 = 0（参照成分）
#   即可完全识别，此时自由参数为 K * p 个（theta = vec(beta_1, ..., beta_K)）。
# ----------------------------------------------------------------------------
bernA_ll_grad <- function(theta, Z, m, x, A_tab) {
  K <- nrow(A_tab) - 1
  n <- length(x); p <- ncol(Z)
  B <- matrix(theta, nrow = p, ncol = K)          # p x K，第 k 列 = beta_k
  eta <- cbind(0, Z %*% B)                        # n x (K+1)，成分 0 为参照
  maxe <- apply(eta, 1, max)
  lam <- exp(eta - maxe); lam <- lam / rowSums(lam)
  A <- A_tab[, x + 1, drop = FALSE]
  P <- t(lam) * A
  dens <- colSums(P)
  diff <- t(P) / dens - lam                       # n x (K+1)
  ll <- sum(lchoose_mx(m, x) + log(dens))
  G <- t(Z) %*% diff[, -1, drop = FALSE]          # p x K
  list(ll = ll, grad = as.numeric(G))
}

fit_bernregA <- function(Z_tr, m, x_tr, K, start = NULL, maxit = 2000,
                         nstart = 1, nindep = 6, seed = 99) {
  p <- ncol(Z_tr)
  A_tab <- bern_coef_table(m, K)
  fn <- function(th) -bernA_ll_grad(th, Z_tr, m, x_tr, A_tab)$ll
  gr <- function(th) -bernA_ll_grad(th, Z_tr, m, x_tr, A_tab)$grad
  starts <- list(rep(0, K * p))                    # 均匀权重起点
  set.seed(seed)
  for (s in seq_len(nindep))
    starts[[length(starts) + 1]] <- rnorm(K * p, 0, 0.5)
  best <- NULL
  for (st in starts) {
    o <- tryCatch(optim(st, fn, gr = gr, method = "BFGS",
                        control = list(maxit = maxit, reltol = 1e-10)),
                  error = function(e) NULL)
    if (!is.null(o) && is.finite(o$value) && (is.null(best) || -o$value > best$ll))
      best <- list(theta = o$par, ll = -o$value, conv = o$convergence)
  }
  theta <- best$theta
  list(npar = K * p, K = K, pz = p, theta = theta, ll = best$ll, conv = best$conv,
       pred = function(Z) {
         p <- ncol(Z)
         B <- matrix(theta, nrow = p, ncol = K)
         eta <- cbind(0, Z %*% B)
         mx <- apply(eta, 1, max)
         lam <- exp(eta - mx); lam <- lam / rowSums(lam)
         list(lambda = lam, B = B)
       })
}

# K 选择（变体 A）：训练集 BIC 曲线 与 训练集内部 5 折交叉验证
bernA_fit_bic <- function(Z_tr, m, x_tr, Kmax = 10, nindep = 6) {
  n1 <- length(x_tr)
  curve <- data.frame(K = 1:Kmax, train_ll = NA, npar = NA, BIC = NA)
  best <- NULL
  for (K in 1:Kmax) {
    fit <- fit_bernregA(Z_tr, m, x_tr, K, nindep = nindep)
    bic <- fit$npar * log(n1) - 2 * fit$ll
    curve[K, c("train_ll", "npar", "BIC")] <- list(fit$ll, fit$npar, bic)
    if (is.null(best) || bic < best$bic) best <- list(K = K, fit = fit, bic = bic)
  }
  best$curve <- curve
  best
}

# 注意：以下 K 选择函数与 test_* 评价函数的第一个矩阵参数，对 bernreg 而言
# 必须是【不含截距】的协变量矩阵 Z；对 logistic / bbreg 则仍是含截距的设计矩阵 X。
# ----------------------------------------------------------------------------
# 4. K 选择（只在训练集上进行）
#    (a) 训练集 BIC：K = 1..Kmax，取 BIC 最小者
#    (b) 训练集内部 5 折交叉验证：验证准则为对数似然，取均值最大者
# ----------------------------------------------------------------------------
bernreg_fit_bic <- function(Z_tr, m, x_tr, Kmax = 10, nstart = 12, nindep = 12) {
  n1 <- length(x_tr)
  best <- NULL
  curve <- data.frame(K = 1:Kmax, train_ll = NA, npar = NA, BIC = NA)
  for (K in 1:Kmax) {
    # 每个 K 都用逻辑回归系数构造的基础起点（fit_bernreg 内部会多起点），
    # 不使用"上一 K 的解"作锚点，避免劣质局部解沿 K 链传播。
    fit <- fit_bernreg(Z_tr, m, x_tr, K, start = NULL, nstart = nstart, nindep = nindep)
    bic <- fit$npar * log(n1) - 2 * fit$ll
    curve$train_ll[K] <- fit$ll; curve$npar[K] <- fit$npar; curve$BIC[K] <- bic
    if (is.null(best) || bic < best$bic) best <- list(K = K, fit = fit, bic = bic)
  }
  best$curve <- curve
  best
}

bernreg_fit_cv <- function(Z_tr, m, x_tr, Kmax = 10, folds = 5, reps = 1,
                           seed = 6789, nstart = 8, nindep = 8) {
  n1 <- length(x_tr)
  set.seed(seed)
  mean_val <- numeric(Kmax)
  for (K in 1:Kmax) {
    v <- c(); start <- NULL
    for (r in 1:reps) {
      foldid <- sample(rep(1:folds, length.out = n1))
      for (f in 1:folds) {
        va <- which(foldid == f); tr <- which(foldid != f)
        fit <- fit_bernreg(Z_tr[tr, , drop = FALSE], m, x_tr[tr], K,
                           start = start, nstart = nstart, nindep = nindep)
        v <- c(v, test_loglik("bernreg", fit, Z_tr[va, , drop = FALSE], m, x_tr[va], K))
        start <- fit$theta
      }
    }
    mean_val[K] <- mean(v)
  }
  K <- which.max(mean_val)
  list(K = K, fit = fit_bernreg(Z_tr, m, x_tr, K, nstart = 12, nindep = nindep),
       mean_val = mean_val)
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
