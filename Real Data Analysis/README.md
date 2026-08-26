# Real Data Analysis 目录说明

本目录存放论文第 5 节"实例数据分析"相关的全部数据、代码与结果。

## 目录结构

| 文件夹 | 内容与状态 |
|---|---|
| **Mobility数据/** | 【推荐入口】Women's Mobility 数据的最终分析：数据、函数库（00_functions.R）、样本内比较（01）、**正确的样本外检验 70/30 训练/测试（02_train_test_outsample.R）**、结果与 README。 |
| **deAyala数据/** | 【推荐入口】deAyala 数学测验数据的最终分析（结构同上）。 |
| **dataset_search/** | 【搜索/筛选阶段】寻找合适数据集的探索过程：通用比较流水线 compare_models.R、26 个候选数据集的筛选脚本与结果、论文作图脚本（export_fig_data.R + make_section5_figures.py）、两份报告。 |
| **大鼠致畸数据/** | 早期尝试（58 窝低铁致畸数据）：样本量小、结构简单，Bern-Bino 优势不显著，**已被 Mobility/deAyala 取代**。 |
| **学生压力数据/**、**自行车租赁数据/** | 更早的不适用尝试（非真正"重复伯努利试验"的二项数据），已弃用。 |

## 五个对比模型

1. `Binomial(m, p)` —— 固定成功概率
2. `Beta-Binomial(m, α, β)` —— p ~ Beta(α, β)
3. `Logit-Normal-Binomial(m)` —— p ~ Logit-Normal(μ, σ)
4. `Kumaraswamy-Binomial(m)` —— p ~ Kumaraswamy(a, b)
5. `Bern-Bino(m, K, λ)` —— p ~ Bernstein(K, λ)（论文提出的模型）

## 两类分析

- **样本内比较**（`01_fit_and_insample.R`）：在全部数据上拟合、同一数据上评价
  （logLik / AIC / BIC / L1 / MSE / 拟合优度自助 p 值），作为参考。
- **样本外检验**（`02_train_test_outsample.R`，正确的流程）：
  1. 把全体 n 个样本**一次性**随机分为训练集 n1（70%）与**完全保留的测试集** n2（30%）；
  2. **只在训练集上**拟合全部模型；
  3. Bern-Bino 的 K **只在训练集上**选择：方法 A 用 BIC；方法 B 用训练集内部的
     5 折 × 3 次交叉验证（4 份拟合、1 份验证，准则为验证似然或 MSE）；
  4. 用候选 K 在训练集上训练 Bern-Bino；
  5. 在**从未参与拟合的测试集**上做最终检验（test logLik / test L1 / test MSE）；
  6. 比较 BIC 与内层 CV 两种 K 选择在测试集上的表现，选定最终 K。

## 快速开始

```r
# 例如 Mobility 数据（进入 Mobility数据 文件夹后）
Rscript 01_fit_and_insample.R        # 样本内比较（参考）
Rscript 02_train_test_outsample.R    # 正确的样本外检验
```

详见各数据文件夹内的 README.md。
