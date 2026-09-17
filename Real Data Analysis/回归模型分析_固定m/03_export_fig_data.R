# NOTE: superseded by 05_report_figures.R, which uses the corrected
# Bern-Bino specification (intercept carried by phi_k, not by beta).
# ============================================================================
# 03_export_fig_data.R --- 导出作图数据（测试集拟合 pmf、先验密度、K 选择曲线）
# ----------------------------------------------------------------------------
# 读取 02 的结果（results/regression_comparison.csv），对每个数据集：
#   1) 按 70/30 划分（seed=12345），只在训练集拟合；
#   2) 取 Bern-Bino 在 BIC 与 CV 两个候选 K 中"测试集 logLik 更优"的那个 K 作为最终 K；
#   3) 计算测试集上各模型（Logistic / Beta-Binomial / Bern-Bino）的平均预测 pmf
#      与测试集经验频率，输出 results/fig_<label>.csv；
#   4) 计算 Bern-Bino 平均先验密度 pi(p)（对测试集协变量取平均），
#      输出 results/prior_<label>.csv；
#   5) K 选择曲线（BIC / CV）已由 02 保存为 results/kgrid_*_<label>.csv。
# 之后由 make_figures.py 作图（风格与数值模拟实验一致）。
# ============================================================================
source("00_regression_lib.R")
suppressMessages(library(psychotools))

std <- function(v) as.numeric(scale(v))
make_design_me <- function(d) cbind(1,
  tests = std(d$tests), genderM = as.numeric(d$gender == "male"),
  study571 = as.numeric(d$study == "571"), semester = std(d$semester))
make_design_pisa00 <- function(d) cbind(1,
  female = d$female, hisei = std(d$hisei), age = std(d$age))

cfg <- list(
  list(file = "data/mathexam_data.csv", label = "MathExam14W", design = make_design_me),
  list(file = "data/pisa00_data.csv",  label = "PISA2000Read", design = make_design_pisa00))

res <- read.csv("results/regression_comparison.csv")
dir.create("results", showWarnings = FALSE)

for (c in cfg) {
  dat <- read.csv(c$file, stringsAsFactors = TRUE)
  X <- c$design(dat); m <- dat$m[1]; x <- dat$x; n <- length(x)
  set.seed(12345)
  test_idx <- sample(1:n, size = floor(0.3 * n)); train_idx <- setdiff(1:n, test_idx)
  X_tr <- X[train_idx, , drop = FALSE]; x_tr <- x[train_idx]
  X_te <- X[test_idx,  , drop = FALSE]; x_te <- x[test_idx]
  n2 <- length(x_te)

  # 最终 K：02 结果中 Bern-Bino(BIC) 与 Bern-Bino(CV) 测试 logLik 更优者
  sub <- res[res$label == c$label & res$Model %in% c("BernBino_BIC", "BernBino_CV"), ]
  Kfin <- sub$K[which.max(sub$test_logLik)]
  cat(sprintf("[%s] final K = %d\n", c$label, Kfin))

  # 训练集拟合
  fit_log <- fit_logistic(X_tr, m, x_tr)
  fit_bb  <- fit_bbreg(X_tr, m, x_tr)
  fit_ber <- fit_bernreg(X_tr, m, x_tr, Kfin, start = NULL, nstart = 12)

  # 测试集平均预测 pmf（对每个 x=0..m 求测试集平均 P(X=x)）
  xs <- 0:m
  pred_log  <- numeric(m + 1)
  pred_bb   <- numeric(m + 1)
  pred_bern <- numeric(m + 1)
  A_tab <- bern_coef_table(m, Kfin)
  pr_ber <- fit_ber$pred(X_te)
  pr_bb  <- fit_bb$pred(X_te)
  pr_log <- fit_log$pred(X_te)
  for (j in seq_along(xs)) {
    v <- xs[j]
    pred_log[j]  <- mean(dbinom(v, m, pr_log))
    pred_bb[j]   <- mean(bbreg_pmf(pr_bb$a, pr_bb$b, m, rep(v, n2)))
    pred_bern[j] <- mean(sapply(seq_len(n2), function(i)
      bernreg_pmf(pr_ber$lambda[i, ], m, v, A_tab)))
  }
  obs <- tabulate(x_te + 1, nbins = m + 1) / n2
  out <- data.frame(x = xs, observed = obs, logistic = pred_log,
                    bbreg = pred_bb, bern = pred_bern)
  write.csv(out, sprintf("results/fig_%s.csv", c$label), row.names = FALSE)

  # 平均先验密度 pi(p) = mean_i sum_k lambda_ik * Beta(k+1, K-k+1)
  pg <- seq(0, 1, length.out = 201)
  prior <- numeric(length(pg))
  for (j in seq_along(pg)) {
    bdens <- dbeta(pg[j], (0:Kfin) + 1, Kfin - (0:Kfin) + 1)   # (K+1) 个 Beta 密度
    prior[j] <- mean(pr_ber$lambda %*% bdens)
  }
  pr <- data.frame(p = pg, mean_prior = prior)
  write.csv(pr, sprintf("results/prior_%s.csv", c$label), row.names = FALSE)
  cat(sprintf("[%s] exported fig_%s.csv, prior_%s.csv\n", c$label, c$label, c$label))
}
cat("Done.\n")
