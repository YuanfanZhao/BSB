# ============================================================================
# 02_交叉验证_样本外比较.R —— deAyala 数学测验数据
# ----------------------------------------------------------------------------
# 目的：在"训练集上拟合、测试集上验证"的框架下比较五个模型，防止过拟合。
# 背景：01 号脚本是在全部数据上拟合并在同一数据上评价（样本内比较）。
#       由于 Bern-Bino 参数较多（K+1 个权重），样本内比较可能高估其表现，
#       因此这里用 K 折交叉验证做样本外（out-of-sample）比较。
#
# 流程（每一步都对应一个编号小节）：
#   0) 加载函数库与数据
#   1) 把全部样本随机划分为 K 折（5 折，重复 3 次，共 15 个"训练/测试"划分）
#   2) 【在训练集上拟合】对每个训练折拟合全部五个模型；
#      Bern-Bino 的 K 也在训练折上用 BIC 重新选择（不使用任何测试信息）
#   3) 【在测试集上验证】用训练折得到的模型，在对应测试折上计算
#      测试对数似然（CV logLik）与测试 L1 / MSE
#   4) 汇总 15 个测试折的结果，得到每个模型的 CV logLik、CV L1、CV MSE
#   5) 保存结果
#
# 运行方式：在脚本所在目录执行
#   Rscript 02_cv_outsample.R
# ============================================================================

# 0) 加载函数库与数据 ------------------------------------------------------
source("00_functions.R")
scores <- read.csv("deayala_data.csv")$score
m <- max(scores)
n <- length(scores)

# 1) 生成训练/测试划分：5 折 x 3 次重复 ------------------------------------
# 每次重复把 n 个样本随机分成 5 折（每折约 n/5 人）。
# 第 r 次重复的第 f 折作为测试集，其余 4 折合并为训练集。
folds <- 5
reps  <- 3
set.seed(12345)                       # 固定随机种子，保证结果可复现
fold_ids <- list()                    # fold_ids[[r]] 是第 r 次重复的折号向量
for (r in 1:reps) {
  fold_ids[[r]] <- sample(rep(1:folds, length.out = n))
}

models <- c("bin", "bb", "ln", "km", "bern")

# 2) & 3) 逐折：训练集拟合 + 测试集验证 ------------------------------------
# 说明：
#   - 训练集拟合：fit_model(model, m, counts_tr) 用训练折计数估计参数；
#     Bern-Bino 用 bern_fit_bic(m, counts_tr) 在训练折上按 BIC 选 K 并拟合。
#   - 测试集验证：用训练得到的 pmf，计算测试折的对数似然
#     test_ll = sum_x c_x^test * log( pmf_train(x) )，
#     以及测试 L1 / MSE（比较测试折经验分布与训练拟合分布）。
per_fold <- list()
for (r in 1:reps) {
  foldid <- fold_ids[[r]]                       # 第 r 次重复的折号
  for (f in 1:folds) {
    test_idx  <- which(foldid == f)             # 本折样本 = 测试集
    train_idx <- which(foldid != f)             # 其余样本 = 训练集

    # 训练/测试折的得分计数
    counts_tr <- tabulate(scores[train_idx] + 1, nbins = m + 1)
    counts_te <- tabulate(scores[test_idx]  + 1, nbins = m + 1)

    # ---- 第 2 步：在训练集上拟合五个模型 ----
    fits_tr <- lapply(models[1:4], function(md) fit_model(md, m, counts_tr))
    names(fits_tr) <- models[1:4]
    bb <- bern_fit_bic(m, counts_tr, Kmax = 30)      # K 在训练折上用 BIC 重选
    fits_tr[["bern"]] <- list(ll = bb$fit$ll, npar = bb$K + 1,
                              pmf = bern_pmf(m, counts_tr, bb$fit$lambda),
                              extra = list(lambda = bb$fit$lambda, K = bb$K))

    # ---- 第 3 步：在测试集上验证 ----
    row <- data.frame(rep = r, fold = f, Kb = bb$K)
    for (md in models) {
      pmf <- fits_tr[[md]]$pmf
      # 测试对数似然：只对有观测的得分求和，避免 0 * log(0)
      test_ll <- sum(ifelse(counts_te > 0, counts_te * log(pmf), 0))
      d <- dist_emp(m, counts_te, pmf)               # 测试 L1 / MSE
      row[[paste0("ll_", md)]]  <- test_ll
      row[[paste0("L1_", md)]]  <- d["L1"]
      row[[paste0("MSE_", md)]] <- d["MSE"]
    }
    per_fold[[length(per_fold) + 1]] <- row
  }
}
res <- do.call(rbind, per_fold)

# 4) 汇总：CV logLik = 15 个测试折对数似然之和（等价于"留一折"式总样本外似然）
#          CV L1 / MSE = 15 个测试折的平均值
summary_tab <- data.frame(
  Model = c("Binomial", "Beta-Binomial", "Logit-Normal", "Kumaraswamy", "Bern-Bino"),
  CV_logLik = sapply(models, function(md) sum(res[[paste0("ll_", md)]])),
  CV_L1     = sapply(models, function(md) mean(res[[paste0("L1_", md)]])),
  CV_MSE    = sapply(models, function(md) mean(res[[paste0("MSE_", md)]]))
)
print(summary_tab, row.names = FALSE)
cat("各训练折 BIC 选定的 K：", paste(res$Kb, collapse = ", "), "\n")

# 5) 保存结果 ---------------------------------------------------------------
write.csv(res,          "results/cv_per_fold.csv",      row.names = FALSE)
write.csv(summary_tab,  "results/cv_out_of_sample.csv", row.names = FALSE)
cat("结果已保存到 results/：cv_per_fold.csv（15 个折的明细），cv_out_of_sample.csv（汇总）\n")
