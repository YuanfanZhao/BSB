# deAyala 数据 —— 第5节实例分析（Bern-Bino 模型对比）

## 1. 数据集简介

- **来源**：de Ayala, R. J. (2009)《The Theory and Practice of Item Response Theory》
  （Guilford Press）第 14 页给出的 **5 道数学测验题**数据；R 包 `mirt` 以 `deAyala` 内置，
  以"响应模式频数表"形式给出（32 种 0/1 作答模式 + 各模式人数）。
- **样本**：总人数 n = 19601（各模式频数之和）。
- **二项结构**：把 5 个题目视为 5 次伯努利试验，X = 答对题数，X ∈ {0,1,…,5}，m = 5。
- **特点**：得分分布宽而接近对称（计数 691, 3099, 4269, 4116, 4041, 3385），
  形态偏离简单 Beta/Logit-Normal 单峰先验，适合展示 Bernstein 先验的灵活性。
- 说明：该数据为 IRT 教材基准数据，包文档未说明其原始采集来源（真实测验或模拟响应），论文中需如实注明。

## 2. 文件说明

| 文件 | 作用 |
|---|---|
| `deayala_data.csv` | 数据：每行一名考生，`score` 为其答对题数（0..5），共 19601 行 |
| `00_functions.R` | **函数库**：五个模型（Binomial/Beta-Binomial/Logit-Normal/Kumaraswamy/Bern-Bino）的对数似然、参数估计（optim / MM 算法）、pmf、模型矩、Pearson 拟合优度、参数自助 p 值、L1/MSE、交叉验证折划分。每个函数均有中文注释。 |
| `01_fit_and_insample.R` | **样本内比较**：读取数据→描述统计与过度离散检验→拟合五个模型（Bern-Bino 在 K=1..30 按 BIC 选 K）→输出样本内比较表 |
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

运行依赖仅为基础 R（base 包）。

## 4. 关键结果（可复现）

**样本内比较（B=999）**

| 模型 | 参数个数 | logLik | AIC | BIC | GOF p |
|---|---|---|---|---|---|
| Binomial | 1 | -35797.74 | 71597.49 | 71605.37 | 0.001 |
| Beta-Binomial | 2 | -33516.01 | 67036.01 | 67051.78 | 0.772 |
| Logit-Normal | 2 | -33536.36 | 67076.72 | 67092.49 | 0.734 |
| Kumaraswamy | 2 | -33533.79 | 67071.57 | 67087.34 | 0.987 |
| **Bern-Bino (K=13)** | 14 | **-33351.85** | **66731.70** | **66842.06** | **1.000** |

**样本外（5 折 × 3 次）**

| 模型 | CV logLik | CV L1 | CV MSE |
|---|---|---|---|
| Binomial | -107398.5 | 0.1978 | 0.006133 |
| Beta-Binomial | -100553.5 | 0.0671 | 0.000548 |
| Logit-Normal | -100614.6 | 0.0711 | 0.000602 |
| Kumaraswamy | -100606.8 | 0.0692 | 0.000580 |
| **Bern-Bino** | **-100139.2** | **0.0383** | **0.000224** |

- Bern-Bino 在样本内与样本外**所有准则上均最优**；15 个训练折全部以 BIC 选定 K=8。
- 结论：优势来自数据真实结构（宽峰、多成分能力分布），不是过拟合。

## 5. 学习指引：训练/测试在哪里发生？

在 `02_cv_outsample.R` 中：
- **第 1 步**（`folds <- 5; reps <- 3; set.seed(12345); fold_ids <- ...`）：
  生成 15 组"训练/测试"划分。
- **第 2 步**（循环内 `counts_tr` 之后，`fits_tr <- lapply(...)` 与 `bern_fit_bic(...)`）：
  在**训练集**上拟合五个模型；Bern-Bino 的 K 也在训练集上用 BIC 重新选择。
- **第 3 步**（`test_ll <- sum(ifelse(counts_te > 0, counts_te * log(pmf), 0))`）：
  在**测试集**上验证：用训练出的 pmf 计算测试集的对数似然（以及测试 L1/MSE）。
- **第 4 步**（`summary_tab`）：汇总 15 个测试折，得到每个模型的 CV logLik 等。

## 6. 与论文的对应

- 论文第 5.2 节 Table 6（样本内）、Table 7（CV）：即 `results/in_sample_comparison.csv`、
  `results/cv_out_of_sample.csv` 中的数字。
- 论文 Figure 8（拟合+残差）、Figure 9（先验密度）、Figure 10（BIC vs K）由
  `../dataset_search/export_fig_data.R`（计算）与 `../dataset_search/make_section5_figures.py`（Python 绘图）生成。
