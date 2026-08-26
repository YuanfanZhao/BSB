# ============================================================================
# 00_模型函数库.R
# ----------------------------------------------------------------------------
# 作用：定义本数据集分析所需的全部模型函数（似然、拟合、拟合优度、交叉验证
#       辅助函数）。本文件只定义函数、不执行分析，由 01/02 脚本调用。
#
# 适用数据形态：固定试验次数 m 的二项计数数据。
#   假设观测为 X_i ~ Binomial(m, p_i)，i = 1..n，p_i 为每次试验的成功概率
#   （随个体变化），X_i 取值 0..m。数据以"得分计数"给出：
#   counts = (c_0, c_1, ..., c_m)，其中 c_x 表示 X = x 的个体数，n = sum(counts)。
#
# 本文件包含的五个模型：
#   1) Binomial(m, p)              ：所有个体共用同一个成功概率 p（无个体异质性）
#   2) Beta-Binomial(m, a, b)      ：p ~ Beta(a, b)
#   3) Logit-Normal-Binomial(m)    ：p ~ Logit-Normal(mu, sigma)
#   4) Kumaraswamy-Binomial(m)     ：p ~ Kumaraswamy(a, b)
#   5) Bern-Bino(m, K, lambda)     ：p ~ Bernstein(K, lambda)（本文提出的模型）
#
# 运行方式：在脚本所在目录执行
#   source("00_模型函数库.R")
# 或由 01/02 脚本自动调用。
# ============================================================================

# ----------------------------------------------------------------------------
# 0. 辅助函数
# ----------------------------------------------------------------------------

# 组合数对数 log C(m, x) = log( m! / (x! (m-x)!) )
lchoose_mx <- function(m, x) {
  lgamma(m + 1) - lgamma(x + 1) - lgamma(m - x + 1)
}

# ----------------------------------------------------------------------------
# 1. Bern-Bino 分布（本文提出的模型）
#    模型：X | p ~ Binomial(m, p)，p ~ Bernstein(K, lambda)
#          Bernstein(K, lambda) 的密度为 sum_{k=0}^K lambda_k * Beta(k+1, K-k+1)
# ----------------------------------------------------------------------------

# 系数矩阵 a_k(x) = B(k+1+x, K-k+1+m-x) / B(k+1, K-k+1)
# 这是"X=x 且先验取第 k 个 Beta 成分"时、去掉组合数后的 Beta-Binomial 概率。
# 返回 (m+1) x (K+1) 矩阵，第 (x+1) 行第 (k+1) 列 = a_k(x)。
bern_coef <- function(m, K) {
  x <- 0:m
  sapply(0:K, function(k)
    exp(lbeta(k + 1 + x, K - k + 1 + m - x) - lbeta(k + 1, K - k + 1)))
}

# Bern-Bino 的概率质量函数 P(X=x) = C(m,x) * sum_k lambda_k * a_k(x)
bern_pmf <- function(m, counts, lambda) {
  K <- length(lambda) - 1
  A <- bern_coef(m, K)
  exp(lchoose_mx(m, 0:m)) * as.numeric(A %*% lambda)
}

# Bern-Bino 对数似然：log L(lambda) = sum_x c_x * log P(X=x)
bern_ll <- function(m, counts, lambda) {
  sum(counts * log(bern_pmf(m, counts, lambda)))
}

# MM 算法估计权重 lambda（论文第 2.2 节，迭代式 2.8）
#   E 步（隐变量分配）：gamma_{kx} = lambda_k * a_k(x) / sum_l lambda_l * a_l(x)
#   M 步（更新权重）  ：lambda_k_new = sum_x c_x * gamma_{kx} / n
# 从均匀权重出发迭代，直到对数似然相对变化小于 eps。
bern_mm <- function(m, counts, K, lambda_init = NULL, maxit = 2000, eps = 1e-10) {
  n <- sum(counts)
  A <- bern_coef(m, K)
  if (is.null(lambda_init)) lam <- rep(1 / (K + 1), K + 1) else lam <- lambda_init
  ll <- bern_ll(m, counts, lam)
  for (it in 1:maxit) {
    # 注意：必须用 matrix(lam, byrow=TRUE) 把权重铺成与 A 同形的矩阵再相乘，
    # 直接用 A * lam 在行数(m+1)与权重数(K+1)不等时会发生错误的向量循环对齐。
    lam_mat <- matrix(lam, nrow = m + 1, ncol = K + 1, byrow = TRUE)
    prop <- A * lam_mat                   # 第 k 列乘 lambda_k：lambda_k * a_k(x)
    gamma <- prop / rowSums(prop)         # 后验分配概率 gamma_{kx}
    lam_new <- colSums(counts * gamma) / n  # MM 更新（保持 sum(lambda)=1）
    ll_new <- bern_ll(m, counts, lam_new)
    if (abs(ll_new - ll) < eps * (abs(ll) + 1e-12)) break
    lam <- lam_new
    ll <- ll_new
  }
  list(lambda = lam, ll = ll, iter = it)   # lambda: 权重向量；ll: 对数似然；iter: 迭代次数
}

# 在 K = 1..Kmax 网格上拟合 Bern-Bino，并用 BIC 选择最优 K。
# BIC(K) = (K+1) * log(n) - 2 * logLik(K)，其中 K+1 为权重参数个数。
# 返回：最优 K、对应拟合对象、最小 BIC 值。
bern_fit_bic <- function(m, counts, Kmax = 30) {
  n <- sum(counts)
  best <- NULL
  for (K in 1:Kmax) {
    f <- bern_mm(m, counts, K)
    bic <- (K + 1) * log(n) - 2 * f$ll
    if (is.null(best) || bic < best$bic) {
      best <- list(K = K, fit = f, bic = bic)
    }
  }
  best
}

# ----------------------------------------------------------------------------
# 2. 对比模型（均通过"先验密度 + 数值积分"构造似然）
# ----------------------------------------------------------------------------

# 2.1 Binomial(m, p)：固定成功概率 p
#     对数似然 = sum_x c_x * log C(m,x) * p^x (1-p)^(m-x)
binom_ll <- function(m, counts, p) {
  sum(counts * dbinom(0:m, m, p, log = TRUE))
}

# 2.2 Beta-Binomial(m, a, b)：p ~ Beta(a, b)
#     对数似然 = sum_x c_x * [ log C(m,x) + log B(a+x, b+m-x) - log B(a,b) ]
#     参数用 (log a, log b) 表示以保证 a, b > 0，optim 在其上优化。
bb_ll <- function(params, m, counts) {
  a <- exp(params[1])
  b <- exp(params[2])
  x <- 0:m
  sum(counts * (lchoose_mx(m, x) + lbeta(a + x, b + m - x) - lbeta(a, b)))
}

# 2.3 Logit-Normal-Binomial：p ~ Logit-Normal(mu, sigma)
#     Logit-Normal 密度：f(p) = phi( (logit(p)-mu)/sigma ) / (sigma p (1-p))
#     为避免 p=0/1 处奇异，用换元 u = logit(p)：
#       P(X=x) = C(m,x) * ∫_{-20}^{20} logit^{-1}(u)^x (1-logit^{-1}(u))^{m-x} * phi(u; mu, sigma) du
f_ln <- function(p, mu, sigma) {
  out <- rep(0, length(p))
  ok <- p > 0 & p < 1
  q <- p[ok]
  out[ok] <- exp(-(log(q / (1 - q)) - mu)^2 / (2 * sigma^2)) /
    (sqrt(2 * pi) * sigma * q * (1 - q))
  out
}
int_ln <- function(x, m, mu, sigma, ppow = 0) {
  # ppow：对 logit^{-1}(u) 再乘的幂，用于计算 E(p^ppow)（见 moments_model）
  f <- function(u) {
    q <- plogis(u)                                   # logit 逆变换
    q^(x + ppow) * (1 - q)^(m - x) * dnorm(u, mu, sigma)
  }
  v <- tryCatch(integrate(f, -20, 20, rel.tol = 1e-8, subdivisions = 400,
                             stop.on.error = FALSE)$value, error = function(e) NA)
  if (!is.finite(v) || v <= 0) 1e-300 else v   # 数值保护：非有限或过小值时给极小正数
}
ln_pmf <- function(m, x, mu, sigma) exp(lchoose_mx(m, x)) * int_ln(x, m, mu, sigma)
# 对数似然（参数 (mu, log sigma)，sigma 取对数保证为正）
ln_ll <- function(params, m, counts) {
  mu <- params[1]
  sigma <- exp(params[2])
  sum(counts * log(sapply(0:m, function(xx) ln_pmf(m, xx, mu, sigma))))
}

# 2.4 Kumaraswamy-Binomial：p ~ Kumaraswamy(a, b)
#     Kumaraswamy 密度：f(p) = a b p^{a-1} (1-p^a)^{b-1}
#     用换元 u = p^a 消除 p=0 处奇异：
#       P(X=x) = C(m,x) * ∫_0^1 u^{x/a} (1-u^{1/a})^{m-x} * b (1-u)^{b-1} du
int_km <- function(x, m, a, b) {
  f <- function(u) u^(x / a) * (1 - u^(1 / a))^(m - x) * b * (1 - u)^(b - 1)
  v <- tryCatch(integrate(f, 0, 1, rel.tol = 1e-8, subdivisions = 400,
                             stop.on.error = FALSE)$value, error = function(e) NA)
  if (!is.finite(v) || v <= 0) 1e-300 else v   # 数值保护
}
km_pmf <- function(m, x, a, b) exp(lchoose_mx(m, x)) * int_km(x, m, a, b)
# 对数似然（参数 (log a, log b)）
km_ll <- function(params, m, counts) {
  a <- exp(params[1])
  b <- exp(params[2])
  sum(counts * log(sapply(0:m, function(xx) km_pmf(m, xx, a, b))))
}

# ----------------------------------------------------------------------------
# 3. 统一拟合接口 fit_model
#    对给定模型估计参数（最大似然），返回：
#       $ll   ：对数似然最大值
#       $npar ：参数个数
#       $pmf  ：拟合的 pmf（长度 m+1，即 P(X=0),...,P(X=m)）
#       $extra：模型参数（供矩、绘图等使用）
# ----------------------------------------------------------------------------
fit_model <- function(model, m, counts, K = NULL) {
  n <- sum(counts)
  switch(model,
    bin = {
      p_hat <- sum((0:m) * counts) / (n * m)         # 固定 p 的 MLE = 总成功次数/总试验次数
      list(ll = binom_ll(m, counts, p_hat), npar = 1,
           pmf = dbinom(0:m, m, p_hat), extra = list(p = p_hat))
    },
    bb = {
      opt <- optim(c(0, 0), function(par) -bb_ll(par, m, counts),
                   method = "L-BFGS-B", lower = c(-8, -8), upper = c(8, 8))
      a <- exp(opt$par[1]); b <- exp(opt$par[2])
      list(ll = -opt$value, npar = 2,
           pmf = exp(lchoose_mx(m, 0:m) + lbeta(a + 0:m, b + m - 0:m) - lbeta(a, b)),
           extra = list(a = a, b = b))
    },
    ln = {
      opt <- optim(c(0, log(0.5)), function(par) -ln_ll(par, m, counts),
                   method = "L-BFGS-B",
                   lower = c(-5, log(0.05)), upper = c(5, log(3)))
      mu <- opt$par[1]; sigma <- exp(opt$par[2])
      list(ll = -opt$value, npar = 2,
           pmf = sapply(0:m, function(xx) ln_pmf(m, xx, mu, sigma)),
           extra = list(mu = mu, sigma = sigma))
    },
    km = {
      opt <- optim(c(0, 0), function(par) -km_ll(par, m, counts),
                   method = "L-BFGS-B", lower = c(-8, -8), upper = c(8, 8))
      a <- exp(opt$par[1]); b <- exp(opt$par[2])
      list(ll = -opt$value, npar = 2,
           pmf = sapply(0:m, function(xx) km_pmf(m, xx, a, b)),
           extra = list(a = a, b = b))
    },
    bern = {
      f <- bern_mm(m, counts, K)
      list(ll = f$ll, npar = K + 1,
           pmf = bern_pmf(m, counts, f$lambda),
           extra = list(lambda = f$lambda, K = K))
    }
  )
}

# ----------------------------------------------------------------------------
# 4. 模型矩与 Pearson 拟合优度
# ----------------------------------------------------------------------------

# 计算混合分布的前两阶矩 mu1 = E(p)，mu2 = E(p^2)。
# 对二项混合模型，X 的方差为
#   Var(X) = m*mu1 - m*mu2 + m*(m-1)*(mu2 - mu1^2)
moments_model <- function(model, fit, m) {
  ex <- fit$extra
  switch(model,
    bin  = c(mu1 = ex$p, mu2 = ex$p^2),
    bb   = c(mu1 = ex$a / (ex$a + ex$b),
             mu2 = ex$a * (ex$a + 1) / ((ex$a + ex$b) * (ex$a + ex$b + 1))),
    ln   = c(mu1 = int_ln(0, 0, ex$mu, ex$sigma, 1),
             mu2 = int_ln(0, 0, ex$mu, ex$sigma, 2)),
    km   = c(mu1 = ex$b * beta(1 + 1 / ex$a, ex$b),
             mu2 = ex$b * beta(1 + 2 / ex$a, ex$b)),
    bern = {
      K <- ex$K; k <- 0:K
      c(mu1 = sum(ex$lambda * (k + 1) / (K + 2)),
        mu2 = sum(ex$lambda * (k + 1) * (k + 2) / ((K + 2) * (K + 3))))
    }
  )
}

# Pearson 统计量：X2 = sum_x c_x * (x - m*mu1)^2 / Var(X)
pearson_stat <- function(m, counts, pmf) {
  n <- sum(counts)
  x <- 0:m
  pmf <- pmf / sum(pmf)                       # 归一化
  mu1 <- sum(x * pmf) / m                     # E(p) = E(X)/m
  mu2 <- sum(x^2 * pmf) / m^2                 # E(p^2) = E(X^2)/m^2
  Var <- m * mu1 - m * mu2 + m * (m - 1) * (mu2 - mu1^2)
  sum(counts * (x - m * mu1)^2 / Var)
}

# 参数自助拟合优度 p 值：
#   在拟合模型下模拟 B 组同规模数据，重新计算 Pearson 统计量，
#   p 值 = (1 + {模拟 X2 >= 观测 X2}) / (B + 1)。
# 说明：这是"样本内"的拟合优度检验（论文第 4.2 节同款）。
gof_p <- function(m, counts, pmf, B = 999, seed = 12345) {
  n <- sum(counts)
  x <- 0:m
  pmf <- pmf / sum(pmf)
  X2_obs <- pearson_stat(m, counts, pmf)
  set.seed(seed)
  cnt <- 0
  for (b in 1:B) {
    cc <- as.numeric(rmultinom(1, n, pmf))     # 从拟合分布模拟一组计数
    mu1 <- sum(x * cc) / (n * m)
    mu2 <- sum(x^2 * cc) / (n * m^2)
    Var <- m * mu1 - m * mu2 + m * (m - 1) * (mu2 - mu1^2)
    X2 <- sum(cc * (x - m * mu1)^2 / Var)
    if (X2 >= X2_obs) cnt <- cnt + 1
  }
  (1 + cnt) / (B + 1)
}

# L1 距离（全变差距离的一半）与 MSE：衡量拟合 pmf 与经验分布的接近程度
dist_emp <- function(m, counts, pmf) {
  emp <- counts / sum(counts)
  pmf <- pmf / sum(pmf)
  c(L1 = 0.5 * sum(abs(emp - pmf)), MSE = mean((emp - pmf)^2))
}

# ----------------------------------------------------------------------------
# 5. 交叉验证辅助函数
# ----------------------------------------------------------------------------

# 把 n 个样本随机划分为 folds 折（每折人数尽量相等）。
# 返回长度为 n 的向量，第 i 个元素为第 i 个样本所在的折号（1..folds）。
make_folds <- function(n, folds, seed) {
  set.seed(seed)
  sample(rep(1:folds, length.out = n))
}
