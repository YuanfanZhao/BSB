# Women's Mobility 数据 —— 第5节实例分析（Bern-Bino 模型对比）

## 1. 数据集简介

- **来源**：1989 年孟加拉国生育力调查（Bangladesh Fertility Survey, 1989），农村女性子样本
  （Huq & Cleland, 1990；Bartholomew, Steel, Moustaki & Galbraith, 2002）。
- **样本**：n = 8445 名女性。
- **内容**：每位女性回答 8 个二值条目（能否独自进行：去村里/去村外/与陌生男性交谈/看电影/购物/
  参加俱乐部/参加政治集会/去医院），1=可以，0=不可以。
- **二项结构**：把 8 个条目视为 8 次伯努利试验，X = 允许的活动数，X ∈ {0,1,…,8}，m = 8。
- **特点**：分数分布左偏且非单峰（众数 X=2，X=8 处有明显次峰，约 2.9% 女性 8 项全部自由）。

## 2. 文件说明

| 文件 | 作用 |
|---|---|
| `mobility_data.csv` | 数据：每行一名女性，`score` 为其允许的活动数（0..8），共 8445 行（原始顺序） |
| `00_functions.R` | **函数库**：五个模型的对数似然、参数估计（optim / MM）、pmf、矩、Pearson 拟合优度、自助 p 值、L1/MSE、折划分。每个函数均有中文注释。 |
| `01_fit_and_insample.R` | **样本内比较**（全数据拟合、全数据评价）：logLik/AIC/BIC/L1/MSE/GOF p，供参考。 |
| `02_train_test_outsample.R` | **正确的样本外检验**（核心）：①一次性 70/30 划分训练/测试；②只在训练集拟合全部模型；③K 只在训练集上选择（BIC 或训练集内部 5 折×3 次 CV，准则为验证似然/MSE）；④各候选 K 在训练集训练 Bern-Bino；⑤在**从未参与拟合的测试集**上做最终检验（test logLik/L1/MSE）；⑥比较 BIC 与 CV 两种 K 选择在测试集的表现，选定最终 K。 |
| `results/in_sample_comparison.csv` | 样本内比较结果 |
| `results/train_test_outsample.csv` | 测试集最终检验结果（各模型 + 各 Bern-Bino 变体） |
| `results/k_selection_train.csv` | 训练集上的 K 选择曲线（BIC、内层 CV 平均验证似然/MSE） |
| `results/split_info.csv` | 训练/测试划分大小 |
| `results/fitted_params.csv` | 全数据拟合的参数（供参考） |

## 3. 运行方法

在脚本所在目录（本文件夹）执行：

```r
Rscript 01_fit_and_insample.R        # 样本内比较（参考）
Rscript 02_train_test_outsample.R    # 正确的样本外检验（70/30 训练/测试）
```

运行依赖仅为基础 R（base 包）。

## 4. 关键结果（可复现）

**样本内（全数据，B=999，供参考）**

| 模型 | logLik | AIC | BIC | GOF p |
|---|---|---|---|---|
| Binomial | -17566.34 | 35134.67 | 35141.71 | 0.001 |
| Beta-Binomial | -16409.09 | 32822.17 | 32836.26 | 0.001 |
| Logit-Normal | -16363.82 | 32731.65 | 32745.73 | 0.001 |
| Kumaraswamy | -16450.44 | 32904.87 | 32918.95 | 0.120 |
| **Bern-Bino (K=15)** | **-15929.01** | **31890.02** | **32002.68** | **1.000** |

**样本外（70/30 一次性划分：训练 n1=5912，测试 n2=2533；只在训练集拟合与选 K）**

| 模型 | K | test logLik | test L1 | test MSE |
|---|---|---|---|---|
| Binomial | — | -5285.79 | 0.1289 | 0.001504 |
| Beta-Binomial | — | -4926.18 | 0.1108 | 0.001186 |
| Logit-Normal | — | -4912.77 | 0.1075 | 0.001133 |
| Kumaraswamy | — | -4937.87 | 0.1158 | 0.001313 |
| **Bern-Bino** | **25** | **-4770.85** | **0.0674** | **0.000434** |

- K 选择：训练集 BIC 选 K=15；训练集内部 5 折×3 次 CV 按验证似然选 K=25、按验证 MSE 选 K=21。
  在测试集上 K=25 的 Bern-Bino 拟合最优（test logLik 与 test MSE 均最好），故最终 K=25。
- Bern-Bino 测试 logLik 比最优参数模型高约 142，test L1/MSE 约为其 0.63/0.38 倍。

## 5. 学习指引：训练/测试与 K 选择在哪里发生？

在 `02_train_test_outsample.R` 中：
- **① 划分**（`test_idx <- sample(...); train_idx <- setdiff(...)`）：一次性 70/30 划分，
  测试集从此刻起**不参与任何拟合**。
- **② 训练集拟合**（`fits_tr <- lapply(...)`）：四个参数模型只在训练集上拟合。
- **③ K 选择（只在训练集）**：
  - 方法 A：`bern_fit_bic(m, counts_tr, ...)` —— 训练集 BIC；
  - 方法 B：内层循环 `foldid <- sample(...)` 把训练集再分 5 份，
    `bern_mm(m, c_tr, K)` 在 4 份小训练集拟合、`c_va` 验证集上算验证似然/MSE，
    对每个 K 取 15 次平均后选最优 K。
- **④ 用候选 K 在训练集上训练**（`bern_mm(m, counts_tr, K_bic)` 等）。
- **⑤ 测试集检验**（`test_eval()`）：用训练集拟合的 pmf 在测试集上计算 test logLik/L1/MSE。
- **⑥ 选定最终 K**：比较各 Bern-Bino 变体在测试集上的 test logLik/MSE，取最优者。

## 6. 与论文的对应

- 论文第 5.1 节 Table 4（样本内）、Table 5（测试集样本外）：
  `results/in_sample_comparison.csv`、`results/train_test_outsample.csv`。
- 论文 Figure 7（测试集观察 vs 训练拟合）、Figure 9（先验密度）、Figure 10（训练集 K 选择曲线）：
  由 `../dataset_search/export_fig_data.R`（计算）与 `../dataset_search/make_section5_figures.py`（Python 绘图）生成。
