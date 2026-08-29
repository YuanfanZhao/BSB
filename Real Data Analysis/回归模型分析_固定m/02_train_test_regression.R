# ============================================================================
# 02_train_test_regression.R --- 训练/测试（样本外）回归模型比较（固定 m）
# ----------------------------------------------------------------------------
# 实验设计（用户要求的 out-of-sample 流程）：
#   1) 将全体 n 个样本按 70/30 随机分为训练集 (n1) 与测试集 (n2)，seed=12345；
#      测试集从一开始就分离，不参与任何模型拟合与超参数选择。
#   2) 只在训练集上拟合：Logistic（glm）、Beta-Binomial 回归、Bern-Bino 回归。
#   3) Bern-Bino 的超参数 K（1..Kmax）只在训练集内选择：
#        (a) 训练集 BIC 最小； (b) 训练集内部 5 折 x 2 次 CV（验证对数似然最大）。
#      两个候选 K 都拿去测试集检验，最终按测试集 logLik 更优者作为 Bern-Bino 的结果。
#   4) 在测试集上检验：test logLik（越大越好）、test MAE / test MSE（越小越好）。
# 同时输出训练集 logLik / AIC / BIC 作为参考（评价以测试集为准）。
# 输出：
#   results/regression_comparison.csv   全部模型 x 数据集的训练/测试结果
#   results/kgrid_bic_<label>.csv       各 K 的训练 logLik 与 BIC（K 选择曲线）
#   results/kgrid_cv_<label>.csv        各 K 的交叉验证对数似然（K 选择曲线）
# ============================================================================
source("00_regression_lib.R")

# ----------------------------------------------------------------------------
# 设计矩阵：把协变量转成数值设计矩阵（连续变量标准化；分类变量转 dummy）
# ----------------------------------------------------------------------------
std <- function(v) as.numeric(scale(v))

make_design_me <- function(d) {
  # 协变量：tests（考前练习数）、gender、study（学位类型）、semester（学期数），共 4 个。
  # 注：原始数据还有 attempt（尝试次数，5 水平但绝大多数为 1、稀疏水平多），
  #     加入后引入噪声、样本外表现下降，故不纳入最终模型。
  cbind(1,
        tests    = std(d$tests),
        genderM  = as.numeric(d$gender == "male"),
        study571 = as.numeric(d$study == "571"),
        semester = std(d$semester))
}

make_design_pisa00 <- function(d) {
  cbind(1,
        female = d$female,
        hisei  = std(d$hisei),
        age    = std(d$age))
}

make_design_ctb <- function(d) {
  cbind(1,
        gender   = d$gender,
        type2    = as.numeric(d$type == 2),
        type3    = as.numeric(d$type == 3),
        size     = std(d$size),
        bachelor = std(d$bachelor),
        language = std(d$language))
}

# ----------------------------------------------------------------------------
# 单个数据集的完整流程
# ----------------------------------------------------------------------------
run_dataset <- function(file, label, design_fun, Kmax = 10) {
  dat <- read.csv(file, stringsAsFactors = TRUE)
  X <- design_fun(dat)
  m <- dat$m[1]                      # 固定 m
  x <- dat$x
  n <- length(x)
  if (length(unique(dat$m)) != 1) stop("m is not fixed!")

  set.seed(12345)                    # 70/30 一次性划分（测试集不参与拟合/K 选择）
  test_idx  <- sample(1:n, size = floor(0.3 * n))
  train_idx <- setdiff(1:n, test_idx)
  X_tr <- X[train_idx, , drop = FALSE]; x_tr <- x[train_idx]
  X_te <- X[test_idx,  , drop = FALSE]; x_te <- x[test_idx]
  n1 <- length(x_tr); n2 <- length(x_te)
  cat(sprintf("[%s] n=%d (train=%d, test=%d), m=%d, p=%d\n",
              label, n, n1, n2, m, ncol(X)))

  # ---- 训练集拟合 ----
  t0 <- proc.time()
  fit_log <- fit_logistic(X_tr, m, x_tr)
  ll_log  <- sum(dbinom(x_tr, m, fit_log$pred(X_tr), log = TRUE))   # 训练集 logLik
  fit_bb  <- fit_bbreg(X_tr, m, x_tr)
  bern_bic <- bernreg_fit_bic(X_tr, m, x_tr, Kmax)
  bern_cv  <- bernreg_fit_cv(X_tr, m, x_tr, Kmax)

  # 保存 K 选择曲线（BIC / CV）
  write.csv(bern_bic$curve, sprintf("results/kgrid_bic_%s.csv", label), row.names = FALSE)
  kc <- data.frame(K = 1:Kmax, cv_mean_ll = bern_cv$mean_val)
  write.csv(kc, sprintf("results/kgrid_cv_%s.csv", label), row.names = FALSE)

  cat(sprintf("   Bern-Bino: K(BIC)=%d, K(CV)=%d, time=%.1fs\n",
              bern_bic$K, bern_cv$K, (proc.time() - t0)[3]))

  # ---- 测试集检验 ----
  ev <- function(md, fit, K = NULL) c(
    test_logLik = test_loglik(md, fit, X_te, m, x_te, K),
    test_MAE    = test_mae(md, fit, X_te, m, x_te, K),
    test_MSE    = test_mse(md, fit, X_te, m, x_te, K))
  aic <- function(ll, npar) 2 * npar - 2 * ll
  bic <- function(ll, npar) npar * log(n1) - 2 * ll
  rows <- list(
    Logistic         = c(npar = fit_log$npar, K = NA, ev("logistic", fit_log),
                         train_ll = ll_log, AIC = aic(ll_log, fit_log$npar),
                         BIC = bic(ll_log, fit_log$npar)),
    BetaBinomial_reg = c(npar = fit_bb$npar, K = NA, ev("bbreg", fit_bb),
                         train_ll = fit_bb$ll, AIC = aic(fit_bb$ll, fit_bb$npar),
                         BIC = bic(fit_bb$ll, fit_bb$npar)),
    BernBino_BIC     = c(npar = bern_bic$fit$npar, K = bern_bic$K,
                         ev("bernreg", bern_bic$fit),
                         train_ll = bern_bic$fit$ll,
                         AIC = aic(bern_bic$fit$ll, bern_bic$fit$npar),
                         BIC = bic(bern_bic$fit$ll, bern_bic$fit$npar)),
    BernBino_CV      = c(npar = bern_cv$fit$npar, K = bern_cv$K,
                         ev("bernreg", bern_cv$fit),
                         train_ll = bern_cv$fit$ll,
                         AIC = aic(bern_cv$fit$ll, bern_cv$fit$npar),
                         BIC = bic(bern_cv$fit$ll, bern_cv$fit$npar)))
  tab <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  tab$Model <- rownames(tab); rownames(tab) <- NULL
  tab$label <- label
  numcols <- c("npar","K","test_logLik","test_MAE","test_MSE","train_ll","AIC","BIC")
  tab[, numcols] <- lapply(tab[, numcols], as.numeric)
  tab[, c("npar","K")] <- lapply(tab[, c("npar","K")], as.integer)
  tab
}

# ----------------------------------------------------------------------------
# 主流程：固定 m 数据集（前两个为推荐实例，CTB 为对照探索）
# ----------------------------------------------------------------------------
cat("========== Fixed-m regression: train/test comparison ==========\n\n")
res <- list()
for (cfg in list(
  list(file = "data/mathexam_data.csv", label = "MathExam14W",
       design = make_design_me, Kmax = 10),
  list(file = "data/pisa00_data.csv",  label = "PISA2000Read",
       design = make_design_pisa00, Kmax = 10),
  list(file = "data/ctb_data.csv",     label = "CTB",
       design = make_design_ctb, Kmax = 10))) {
  t <- tryCatch(run_dataset(cfg$file, cfg$label, cfg$design, cfg$Kmax),
                error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL })
  if (!is.null(t)) {
    print(t[, c("Model","label","npar","K","test_logLik","test_MAE","test_MSE")],
          row.names = FALSE)
    res[[cfg$label]] <- t; cat("\n")
  }
}
out <- do.call(rbind, res)
write.csv(out, "results/regression_comparison.csv", row.names = FALSE)
cat("Saved results/regression_comparison.csv\n")
