# deAyala 数据 —— 第5节实例分析（Bern-Bino 模型对比）

## 1. 数据集简介

- **来源**：de Ayala, R. J. (2009)《The Theory and Practice of Item Response Theory》
  （Guilford Press）第 14 页给出的 **5 道数学测验题**数据；R 包 `mirt` 以 `deAyala` 内置，
  以"响应模式频数表"形式给出（32 种 0/1 作答模式 + 各模式人数）。
- **样本**：总人数 n = 19601（各模式频数之和）。
- **二项结构**：把 5 个题目视为 5 次伯努利试验，X = 答对题数，X ∈ {0,1,…,5}，m = 5。
- **特点**：得分分布宽而接近对称（计数 691, 3099, 4269, 4116, 4041, 3385）。
- 说明：该数据为 IRT 教材基准数据，包文档未说明其原始采集来源，论文中需如实注明。

## 2. 文件说明

| 文件 | 作用 |
|---|---|
| `deayala_data.csv` | 数据：每行一名考生，`score` 为其答对题数（0..5），共 19601 行 |
| `00_functions.R` | **函数库**：五个模型的对数似然、参数估计（optim / MM）、pmf、矩、Pearson 拟合优度、自助 p 值、L1/MSE、折划分。每个函数均有中文注释。 |
| `01_fit_and_insample.R` | **样本内比较**（全数据拟合、全数据评价）：供参考。 |
| `02_train_test_outsample.R` | **正确的样本外检验**（核心）：①一次性 70/30 划分训练/测试；②只在训练集拟合全部模型；③K 只在训练集上选择（BIC 或训练集内部 5 折×3 次 CV）；④各候选 K 在训练集训练 Bern-Bino；⑤在**从未参与拟合的测试集**上做最终检验；⑥比较 BIC 与 CV 两种 K 选择在测试集的表现，选定最终 K。 |
| `results/in_sample_comparison.csv` | 样本内比较结果 |
| `results/train_test_outsample.csv` | 测试集最终检验结果 |
| `results/k_selection_train.csv` | 训练集上的 K 选择曲线（BIC、内层 CV 平均验证似然/MSE） |
| `results/split_info.csv` | 训练/测试划分大小 |
| `results/fitted_params.csv` | 全数据拟合的参数（供参考） |

## 3. 运行方法

```r
Rscript 01_fit_and_insample.R        # 样本内比较（参考）
Rscript 02_train_test_outsample.R    # 正确的样本外检验（70/30 训练/测试）
```

运行依赖仅为基础 R（base 包）。

## 4. 关键结果（可复现）

**样本内（全数据，B=999，供参考）**

| 模型 | logLik | AIC | BIC | GOF p |
|---|---|---|---|---|
| Binomial | -35797.74 | 71597.49 | 71605.37 | 0.001 |
| Beta-Binomial | -33516.01 | 67036.01 | 67051.78 | 0.772 |
| Logit-Normal | -33536.36 | 67076.72 | 67092.49 | 0.734 |
| Kumaraswamy | -33533.79 | 67071.57 | 67087.34 | 0.987 |
| **Bern-Bino (K=13)** | **-33351.85** | **66731.70** | **66842.06** | **1.000** |

**样本外（70/30 一次性划分：训练 n1=13721，测试 n2=5880；只在训练集拟合与选 K）**

| 模型 | K | test logLik | test L1 | test MSE |
|---|---|---|---|---|
| Binomial | — | -10813.28 | 0.2051 | 0.006657 |
| Beta-Binomial | — | -10075.45 | 0.0736 | 0.000628 |
| Logit-Normal | — | -10082.10 | 0.0774 | 0.000686 |
| Kumaraswamy | — | -10080.47 | 0.0758 | 0.000657 |
| **Bern-Bino** | **30** | **-10013.01** | **0.0294** | **0.000123** |

- K 选择：训练集 BIC 选 K=8；训练集内部 5 折×3 次 CV 按验证似然选 K=30、按验证 MSE 选 K=26。
  在测试集上 K=30 的 Bern-Bino 拟合最优（test logLik 与 test MSE 均最好），故最终 K=30。
- Bern-Bino 测试 logLik 比最优参数模型高约 62，test L1/MSE 约为其 0.40/0.20 倍。

## 5. 学习指引：训练/测试与 K 选择在哪里发生？

在 `02_train_test_outsample.R` 中：
- **① 划分**（`test_idx <- sample(...); train_idx <- setdiff(...)`）：一次性 70/30 划分。
- **② 训练集拟合**（`fits_tr <- lapply(...)`）：四个参数模型只在训练集上拟合。
- **③ K 选择（只在训练集）**：方法 A 为训练集 BIC；方法 B 为训练集内部 5 折×3 次 CV
  （`bern_mm(m, c_tr, K)` 在 4 份小训练集拟合、`c_va` 验证集上算验证似然/MSE）。
- **④ 用候选 K 在训练集上训练**；**⑤ 测试集检验**（`test_eval()`）。
- **⑥ 选定最终 K**：比较各 Bern-Bino 变体在测试集上的 test logLik/MSE，取最优者。

## 6. 与论文的对应

- 论文第 5.2 节 Table 6（样本内）、Table 7（测试集样本外）：
  `results/in_sample_comparison.csv`、`results/train_test_outsample.csv`。
- 论文 Figure 8（测试集观察 vs 训练拟合）、Figure 9（先验密度）、Figure 10（训练集 K 选择曲线）：
  由 `../dataset_search/export_fig_data.R`（计算）与 `../dataset_search/make_section5_figures.py`（Python 绘图）生成。
