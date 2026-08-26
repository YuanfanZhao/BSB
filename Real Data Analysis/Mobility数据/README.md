# Women's Mobility 数据 —— 第5节实例分析（Bern-Bino 模型对比）

## 1. 数据集简介

- **来源**：1989 年孟加拉国生育力调查（Bangladesh Fertility Survey, 1989），农村女性子样本
  （Huq & Cleland, 1990；Bartholomew, Steel, Moustaki & Galbraith, 2002）。
- **样本**：n = 8445 名女性。
- **内容**：每位女性回答 8 个二值条目（能否独自进行：去村里/去村外/与陌生男性交谈/看电影/购物/
  参加俱乐部/参加政治集会/去医院），1=可以，0=不可以。
- **二项结构**：把 8 个条目视为 8 次伯努利试验，X = 允许的活动数，X ∈ {0,1,…,8}，m = 8。
- **特点**：分数分布左偏且非单峰（众数 X=2，X=8 处有明显次峰，约 2.9% 女性 8 项全部自由），
  存在明显个体异质性（过度离散），适合检验 Bernstein 先验的灵活性。

## 2. 文件说明

| 文件 | 作用 |
|---|---|
| `mobility_data.csv` | 数据：每行一名女性，`score` 为其允许的活动数（0..8），共 8445 行（原始顺序） |
| `00_functions.R` | **函数库**：五个模型（Binomial/Beta-Binomial/Logit-Normal/Kumaraswamy/Bern-Bino）的对数似然、参数估计（optim / MM 算法）、pmf、模型矩、Pearson 拟合优度、参数自助 p 值、L1/MSE、交叉验证折划分。每个函数均有中文注释。 |
| `01_fit_and_insample.R` | **样本内比较**：读取数据→描述统计与过度离散检验→拟合五个模型（Bern-Bino 在 K=1..30 按 BIC 选 K）→输出样本内比较表（logLik/AIC/BIC/L1/MSE/GOF p） |
| `02_cv_outsample.R` | **样本外（交叉验证）比较**：5 折 × 3 次的训练/测试划分；每个训练折上拟合全部模型（K 在训练折上 BIC 重选），每个测试折上计算测试 logLik / L1 / MSE |
| `results/in_sample_comparison.csv` | 样本内比较结果 |
| `results/cv_out_of_sample.csv` | 交叉验证汇总结果（CV logLik / CV L1 / CV MSE） |
| `results/cv_per_fold.csv` | 15 个训练/测试折的逐折明细 |
| `results/k_selection.csv` | Bern-Bino 的 K 选择曲线（K, logLik, AIC, BIC） |
| `results/fitted_params.csv` | 各模型估计出的参数 |

## 3. 运行方法

在脚本所在目录（本文件夹）执行：

```r
Rscript 01_fit_and_insample.R    # 样本内比较
Rscript 02_cv_outsample.R        # 样本外（交叉验证）比较
```

两个脚本都会把结果写入 `results/` 目录。运行依赖仅为基础 R（base 包），无需额外安装包。

## 4. 关键结果（可复现）

**样本内比较（B=999）**

| 模型 | 参数个数 | logLik | AIC | BIC | GOF p |
|---|---|---|---|---|---|
| Binomial | 1 | -17566.34 | 35134.67 | 35141.71 | 0.001 |
| Beta-Binomial | 2 | -16409.09 | 32822.17 | 32836.26 | 0.001 |
| Logit-Normal | 2 | -16363.82 | 32731.65 | 32745.73 | 0.001 |
| Kumaraswamy | 2 | -16450.44 | 32904.87 | 32918.95 | 0.120 |
| **Bern-Bino (K=15)** | 16 | **-15929.01** | **31890.02** | **32002.68** | **1.000** |

**样本外（5 折 × 3 次）**

| 模型 | CV logLik | CV L1 | CV MSE |
|---|---|---|---|
| Binomial | -52704.25 | 0.1246 | 0.001370 |
| Beta-Binomial | -49237.36 | 0.1022 | 0.000966 |
| Logit-Normal | -49101.33 | 0.1007 | 0.000904 |
| Kumaraswamy | -49361.93 | 0.1094 | 0.001102 |
| **Bern-Bino** | **-47795.48** | **0.0564** | **0.000311** |

- Bern-Bino 在样本内与样本外**所有准则上均最优**；15 个训练折全部以 BIC 选定 K=15。
- 结论：优势来自数据真实结构（X=8 次峰），不是过拟合。

## 5. 学习指引：训练/测试在哪里发生？

在 `02_cv_outsample.R` 中：
- **第 1 步**（`folds <- 5; reps <- 3; set.seed(12345); fold_ids <- ...`）：
  生成 15 组"训练/测试"划分（随机 5 折，每折作为一次测试集，其余为训练集）。
- **第 2 步**（循环内 `counts_tr` 之后，`fits_tr <- lapply(...)` 与 `bern_fit_bic(...)`）：
  在**训练集**上拟合五个模型；注意 Bern-Bino 的 K 也在训练集上用 BIC 重新选择，
  **绝不使用测试集信息**。
- **第 3 步**（`test_ll <- sum(ifelse(counts_te > 0, counts_te * log(pmf), 0))`）：
  在**测试集**上验证：用训练出的 pmf 计算测试集的对数似然（以及测试 L1/MSE）。
- **第 4 步**（`summary_tab`）：汇总 15 个测试折，得到每个模型的 CV logLik 等。

## 6. 与论文的对应

- 论文第 5.1 节 Table 4（样本内）、Table 5（CV）：即 `results/in_sample_comparison.csv`、
  `results/cv_out_of_sample.csv` 中的数字。
- 论文 Figure 7（拟合+残差）、Figure 9（先验密度）、Figure 10（BIC vs K）由
  `../dataset_search/export_fig_data.R`（计算）与 `../dataset_search/make_section5_figures.py`（Python 绘图）生成。
