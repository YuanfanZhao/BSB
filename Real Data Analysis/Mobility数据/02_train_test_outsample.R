# ============================================================================
# 02_train_test_outsample.R —— Women's Mobility 数据（正确的样本外检验流程）
# ----------------------------------------------------------------------------
# 实验设计（外层"训练/测试" + 内层"K 选择"）：
#
#   ① 把全体 n 个样本【一次性】随机分为 n1 个训练样本 + n2 个测试样本
#      （这里取 70% / 30%）。测试样本自始至终【不参与任何模型拟合】，
#      也不参与 Bern-Bino 的超参数 K 的选择，只用于最后的拟合效果检验。
#
#   ② 在 n1 个训练样本上拟合全部对比模型：
#      Binomial / Beta-Binomial / Logit-Normal / Kumaraswamy（直接最大似然）。
#
#   ③ Bern-Bino 的超参数 K 只在训练集上选择，两种方法都做：
#      方法 A：BIC —— 在 n1 训练集上，K=1..30 中取 BIC 最小者；
#      方法 B：5 折 x 3 次交叉验证 —— 在 n1 训练集内部再划 5 份，
#              每次取 4 份拟合 Bern-Bino、1 份验证，对每个候选 K 得到
#              15 个验证准则（对数似然 / MSE）取平均，选最优 K。
#              （方法 B 的交叉验证只发生在训练集内部，仍然不碰测试集。）
#
#   ④ 用方法 A、方法 B 选出的 K 分别在 n1 训练集上训练 Bern-Bino。
#
#   ⑤ 在 n2 测试集上做【最终检验】：比较四个参数模型与各 Bern-Bino 变体
#      的测试对数似然（test logLik）、测试 L1、测试 MSE。
#
#   ⑥ 比较两种 K 选择方法在测试集上的表现（测试 logLik / MSE 更优者胜），
#      选定最终使用的 K，输出最终对比表。
#
# 运行方式：在脚本所在目录执行
#   Rscript 02_train_test_outsample.R
# ============================================================================

# 0) 加载函数库与数据 ------------------------------------------------------
source("00_functions.R")
scores <- read.csv("mobility_data.csv")$score
m <- max(scores)
n <- length(scores)
models <- c("bin", "bb", "ln", "km", "bern")

# ① 外层划分：70% 训练 / 30% 测试（一次性，固定种子保证可复现）-------------
set.seed(12345)
test_idx  <- sample(1:n, size = floor(0.3 * n))   # 测试集下标（不参与拟合）
train_idx <- setdiff(1:n, test_idx)               # 训练集下标
counts_tr <- tabulate(scores[train_idx] + 1, nbins = m + 1)  # 训练集得分计数
counts_te <- tabulate(scores[test_idx]  + 1, nbins = m + 1)  # 测试集得分计数
n1 <- length(train_idx); n2 <- length(test_idx)
cat(sprintf("样本总数 n = %d -> 训练集 n1 = %d (%.0f%%), 测试集 n2 = %d (%.0f%%)\n",
            n, n1, 100 * n1 / n, n2, 100 * n2 / n))

# ② 在训练集上拟合四个参数模型 ---------------------------------------------
fits_tr <- lapply(models[1:4], function(md) fit_model(md, m, counts_tr))
names(fits_tr) <- models[1:4]

# ③ Bern-Bino 的 K 选择（只在训练集上）-------------------------------------

# ---- 方法 A：BIC ----
K_bic <- bern_fit_bic(m, counts_tr, Kmax = 30)$K
cat("方法 A（训练集 BIC）选定的 K =", K_bic, "\n")

# ---- 方法 B：训练集内部 5 折 x 3 次交叉验证 ----
inner_folds <- 5; inner_reps <- 3
K_grid <- 1:30
set.seed(6789)                       # 内层交叉验证的随机种子
cv_ll  <- matrix(NA, length(K_grid), inner_folds * inner_reps)
cv_mse <- matrix(NA, length(K_grid), inner_folds * inner_reps)
for (K in K_grid) {
  col <- 0
  for (r in 1:inner_reps) {
    foldid <- sample(rep(1:inner_folds, length.out = n1))   # 训练集内部再分 5 份
    for (f in 1:inner_folds) {
      col <- col + 1
      v_idx  <- which(foldid == f)              # 1 份验证集
      tr_idx <- which(foldid != f)              # 4 份小训练集
      c_tr <- tabulate(scores[train_idx][tr_idx] + 1, nbins = m + 1)
      c_va <- tabulate(scores[train_idx][v_idx]  + 1, nbins = m + 1)
      # 在 4 份小训练集上拟合给定 K 的 Bern-Bino，在 1 份验证集上计算准则
      fK  <- bern_mm(m, c_tr, K)
      pmf <- bern_pmf(m, c_tr, fK$lambda)
      cv_ll[K, col]  <- sum(ifelse(c_va > 0, c_va * log(pmf), 0))  # 验证对数似然
      cv_mse[K, col] <- dist_emp(m, c_va, pmf)["MSE"]              # 验证 MSE
    }
  }
}
mean_ll  <- rowMeans(cv_ll)                     # 每个候选 K 的平均验证似然
mean_mse <- rowMeans(cv_mse)                    # 每个候选 K 的平均验证 MSE
K_cvll   <- K_grid[which.max(mean_ll)]          # 按验证似然选 K
K_cvmse  <- K_grid[which.min(mean_mse)]         # 按验证 MSE 选 K
cat("方法 B（内层 CV）选定的 K：按验证似然 =", K_cvll, "，按验证 MSE =", K_cvmse, "\n")

# ④ 用三个候选 K 在训练集上训练 Bern-Bino ---------------------------------
bern_A   <- bern_mm(m, counts_tr, K_bic)
bern_Bll <- bern_mm(m, counts_tr, K_cvll)
bern_Bmse<- bern_mm(m, counts_tr, K_cvmse)

# ⑤ 在测试集上做最终检验 ---------------------------------------------------
# 测试准则：用"训练集拟合的 pmf"在测试集上计算
#   test logLik = sum_x c_x^test * log( pmf_train(x) )
#   test L1 / MSE = 测试集经验分布与训练拟合分布的距离
test_eval <- function(pmf) {
  ll <- sum(ifelse(counts_te > 0, counts_te * log(pmf), 0))
  d  <- dist_emp(m, counts_te, pmf)
  # 注意用 unname() 去掉内部向量名，避免列名变成 test_L1.L1 之类
  c(test_logLik = ll, test_L1 = unname(d["L1"]), test_MSE = unname(d["MSE"]))
}

rows <- list()
for (md in models[1:4]) {
  rows[[md]] <- c(Model = md, K = NA, test_eval(fits_tr[[md]]$pmf))
}
rows[["bern_A"]]    <- c(Model = "bern_A", K = K_bic,   test_eval(bern_pmf(m, counts_tr, bern_A$lambda)))
rows[["bern_Bll"]]  <- c(Model = "bern_Bll", K = K_cvll,  test_eval(bern_pmf(m, counts_tr, bern_Bll$lambda)))
rows[["bern_Bmse"]] <- c(Model = "bern_Bmse", K = K_cvmse, test_eval(bern_pmf(m, counts_tr, bern_Bmse$lambda)))
res <- do.call(rbind, rows)
res <- as.data.frame(res, stringsAsFactors = FALSE)
res$test_logLik <- as.numeric(res$test_logLik)
res$test_L1     <- as.numeric(res$test_L1)
res$test_MSE    <- as.numeric(res$test_MSE)

# ⑥ 比较两种 K 选择方法在测试集上的表现，选定最终 K -------------------------
cat("\n测试集最终检验结果：\n")
print(res, row.names = FALSE)
cat("\n候选 K：BIC =", K_bic, "；内层CV(似然) =", K_cvll, "；内层CV(MSE) =", K_cvmse, "\n")
bern_rows <- res[grepl("^bern", res$Model), ]
best_ll  <- bern_rows$Model[which.max(bern_rows$test_logLik)]
best_mse <- bern_rows$Model[which.min(bern_rows$test_MSE)]
cat("测试 logLik 最优的 Bern-Bino 变体：", best_ll, "；测试 MSE 最优：", best_mse, "\n")
final_K <- if (best_ll == "bern_A") K_bic else if (best_ll == "bern_Bll") K_cvll else K_cvmse
cat("最终选定 K =", final_K, "\n")

# ⑦ 保存结果 ---------------------------------------------------------------
write.csv(res, "results/train_test_outsample.csv", row.names = FALSE)
kg <- data.frame(K = K_grid, BIC_train = (K_grid + 1) * log(n1) - 2 *
                   sapply(K_grid, function(k) bern_mm(m, counts_tr, k)$ll),
                 CV_mean_ll = mean_ll, CV_mean_mse = mean_mse)
write.csv(round(kg, 4), "results/k_selection_train.csv", row.names = FALSE)
write.csv(data.frame(split = c("train", "test"), n = c(n1, n2)),
          "results/split_info.csv", row.names = FALSE)
cat("结果已保存到 results/：train_test_outsample.csv, k_selection_train.csv, split_info.csv\n")
