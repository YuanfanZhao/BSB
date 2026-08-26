# ============================================================================
# 01_模型拟合与样本内比较.R —— Women's Mobility 数据
# ----------------------------------------------------------------------------
# 数据：孟加拉国生育力调查 1989（Bangladesh Fertility Survey, 1989）农村女性
#       子样本，n = 8445 人；每位女性回答 8 个二值条目（能否独自进行某项活动）。
#       把 8 个条目看作 8 次伯努利试验，X = 允许的活动数，取值 0..8，m = 8。
#
# 本脚本做什么：
#   1) 读取数据（得分计数表），给出描述统计与过度离散检验；
#   2) 用最大似然拟合五个模型：Binomial、Beta-Binomial、Logit-Normal、
#      Kumaraswamy、Bern-Bino（Bern-Bino 的 K 在 1..30 网格上用 BIC 选定）；
#   3) 输出"样本内"比较表：logLik / AIC / BIC / L1 / MSE / 拟合优度 p 值；
#   4) 把结果保存到 results/ 目录。
#
# 运行方式：在脚本所在目录执行
#   Rscript 01_模型拟合与样本内比较.R
# ============================================================================

# 0) 加载函数库（定义全部模型函数）与数据 -------------------------------
source("00_functions.R")

# 读取数据：每行一名女性，score 为其允许的活动数（0..8）
scores <- read.csv("mobility_data.csv")$score
m <- max(scores)                   # 试验次数 m = 8
counts <- tabulate(scores + 1, nbins = m + 1)  # counts[x+1] = 得分为 x 的人数
n <- length(scores)                # 总人数 n = 8445
cat("n =", n, " m =", m, "\n")
cat("得分分布：", paste(counts, collapse = ", "), "\n")

# 1) 数据描述：固定 p 的过度离散统计量 -----------------------------------
# 若所有个体共用同一成功概率 p（MLE p.hat = 总成功次数/总试验次数），
# 则 phi = sum_i (X_i - m*p.hat)^2 / (m*p.hat*(1-p.hat)) / (n-1) 近似服从
# 卡方分布（自由度为 n-1）。phi >> 1 说明个体间存在异质性（过度离散）。
p_hat <- sum((0:m) * counts) / (n * m)
phi   <- sum(counts * ((0:m) - m * p_hat)^2 / (m * p_hat * (1 - p_hat))) / (n - 1)
p_phi <- pchisq((n - 1) * phi, df = n - 1, lower.tail = FALSE)
cat(sprintf("p_hat = %.4f,  过度离散统计量 phi = %.3f (p = %.3g)\n", p_hat, phi, p_phi))

# 2) 拟合五个模型 ---------------------------------------------------------
models <- c("bin", "bb", "ln", "km", "bern")
fits <- list()

# 2.1 四个参数模型：fit_model 内部用 optim 做最大似然
for (md in models[1:4]) {
  fits[[md]] <- fit_model(md, m, counts)
}

# 2.2 Bern-Bino：在 K=1..30 上做 MM 算法，按 BIC 选 K
bern_best <- bern_fit_bic(m, counts, Kmax = 30)
fits[["bern"]] <- list(ll = bern_best$fit$ll, npar = bern_best$K + 1,
                       pmf = bern_pmf(m, counts, bern_best$fit$lambda),
                       extra = list(lambda = bern_best$fit$lambda, K = bern_best$K))
cat("Bern-Bino: BIC 选定的 K =", bern_best$K, "\n")
cat("权重 lambda =", paste(round(bern_best$fit$lambda, 4), collapse = ", "), "\n")

# 3) 样本内比较表 ----------------------------------------------------------
# 准则说明：
#   logLik : 对数似然（越大越好）
#   AIC    : 2*npar - 2*logLik（越小越好，对参数个数有惩罚）
#   BIC    : npar*log(n) - 2*logLik（越小越好，惩罚比 AIC 更重）
#   L1     : 拟合 pmf 与经验分布的全变差距离的一半（越小越好）
#   MSE    : 拟合 pmf 与经验分布的均方误差（越小越好）
#   GOF p  : Pearson 拟合优度的参数自助 p 值（越大越好，>0.05 表示不被拒绝）
dist_vals <- sapply(models, function(md) dist_emp(m, counts, fits[[md]]$pmf))
tab <- data.frame(
  Model  = c("Binomial", "Beta-Binomial", "Logit-Normal", "Kumaraswamy",
             sprintf("Bern-Bino (K=%d)", bern_best$K)),
  K      = c(NA, NA, NA, NA, bern_best$K),
  npar   = sapply(fits, function(f) f$npar),
  logLik = sapply(fits, function(f) round(f$ll, 2)),
  AIC    = sapply(fits, function(f) round(2 * f$npar - 2 * f$ll, 2)),
  BIC    = sapply(fits, function(f) round(f$npar * log(n) - 2 * f$ll, 2)),
  L1     = round(dist_vals["L1", ], 4),
  MSE    = round(dist_vals["MSE", ], 6),
  GOFp   = sapply(fits, function(f) round(gof_p(m, counts, f$pmf), 3)),
  row.names = NULL
)
print(tab)
write.csv(tab, "results/in_sample_comparison.csv", row.names = FALSE)

# 4) 保存 Bern-Bino 的 K 选择曲线（BIC/AIC vs K）与拟合参数 ---------------
kg <- data.frame(K = 1:30, logLik = NA, AIC = NA, BIC = NA)
for (K in 1:30) {
  f <- bern_mm(m, counts, K)
  kg$logLik[K] <- f$ll
  kg$AIC[K]    <- 2 * (K + 1) - 2 * f$ll
  kg$BIC[K]    <- (K + 1) * log(n) - 2 * f$ll
}
write.csv(round(kg, 3), "results/k_selection.csv", row.names = FALSE)

params <- data.frame(
  Model = c("Binomial", "Beta-Binomial", "Logit-Normal", "Kumaraswamy", "Bern-Bino"),
  Estimated = c(
    sprintf("p = %.4f", fits$bin$extra$p),
    sprintf("a = %.4f, b = %.4f", fits$bb$extra$a, fits$bb$extra$b),
    sprintf("mu = %.4f, sigma = %.4f", fits$ln$extra$mu, fits$ln$extra$sigma),
    sprintf("a = %.4f, b = %.4f", fits$km$extra$a, fits$km$extra$b),
    sprintf("K = %d, lambda = (%s)", bern_best$K,
            paste(round(bern_best$fit$lambda, 4), collapse = ", "))))
write.csv(params, "results/fitted_params.csv", row.names = FALSE)

cat("结果已保存到 results/ 目录：in_sample_comparison.csv, k_selection.csv, fitted_params.csv\n")
