# Real Data Analysis 目录说明

本目录存放论文第 5 节"实例数据分析"相关的全部数据、代码与结果。

## 目录结构

| 文件夹 | 内容与状态 |
|---|---|
| **Mobility数据/** | 【推荐入口】Women's Mobility 数据的最终分析：数据、函数库（00_functions.R）、样本内比较（01）、交叉验证样本外比较（02）、结果与 README。 |
| **deAyala数据/** | 【推荐入口】deAyala 数学测验数据的最终分析：数据、函数库（00_functions.R）、样本内比较（01）、交叉验证样本外比较（02）、结果与 README。 |
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
  （logLik / AIC / BIC / L1 / MSE / 拟合优度自助 p 值）。
- **样本外比较**（`02_cv_outsample.R`）：5 折 × 3 次交叉验证——训练集拟合
  （Bern-Bino 的 K 在训练集上用 BIC 重新选择），测试集验证（CV logLik / CV L1 / CV MSE），
  用于防止过拟合。

## 快速开始

```r
# 例如 Mobility 数据（进入 Mobility数据 文件夹后）
Rscript 01_fit_and_insample.R   # 样本内比较
Rscript 02_cv_outsample.R       # 样本外（交叉验证）比较
```

详见各数据文件夹内的 README.md。
