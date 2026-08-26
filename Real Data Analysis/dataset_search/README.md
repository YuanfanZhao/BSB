# dataset_search 目录说明（数据集搜索与筛选阶段）

本文件夹是**寻找合适实例数据集**的探索阶段的全部产出，也是论文图的生成工具所在。
**正式的、可复现的每数据集分析请使用 `../Mobility数据/` 与 `../deAyala数据/` 两个文件夹。**

## 文件分类

### 1. 通用五模型比较流水线（核心）
- `compare_models.R`：统一的模型比较函数库（Binomial/Beta-Binomial/Logit-Normal/
  Kumaraswamy/Bern-Bino 的拟合、logLik/AIC/BIC/L1/MSE/GOF、K 选择、交叉验证），
  可被其他脚本 `source()` 使用。`Mobility数据/00_functions.R` 与 `deAyala数据/00_functions.R`
  是其"固定 m、教学注释版"。

### 2. 数据集搜索与筛选（探索阶段）
- `extract_datasets.R`：从 R 包中提取候选二项数据集（N, R）并计算过度离散度。
- `run_screening.R`：对全部候选数据集运行 `compare_one` 比较，输出各数据集"Bern-Bino 是否胜出"。
- `extra_tests.R` / `recompute_fastm.R`：补充候选的尝试与修正后的重算。
- `datasets_cache.rds` / `screening_results.rds`：搜索阶段的缓存与结果。

### 3. Mobility / deAyala 的早期探索脚本（已被正式文件夹取代）
- `mobility_analysis.R`、`mobility_figures.R`、`mobility_Kgrid.R`、`dataset_info.R`。
- `mobility_scores.csv`、`mobility_items.csv`、`mobility_full_comparison.csv`、`mobility_K_grid.csv`。
- `mobility_fit.png`、`mobility_prior_density.pdf`、`deayala_fit.png`、`deayala_prior_density.pdf`。

### 4. 交叉验证（早期版本，已被取代）
- `cv_validation.R`：早期"反复 5 折、每折轮流当测试"的样本外比较，结果见
  `cv_results.rds`、`cv_out_of_sample_results.csv`。
- **注意**：最终采用更严格的"一次性 70/30 训练/测试 + 训练集内选 K"流程，
  见 `Mobility数据/02_train_test_outsample.R` 与 `deAyala数据/02_train_test_outsample.R`。

### 5. 论文作图（Figure 7–10，与最终协议一致）
- `export_fig_data.R`：使用与 02 相同的 70/30 划分（seed=12345），**只在训练集拟合**、
  用最终 K（Mobility=25、deAyala=30）拟合 Bern-Bino，导出：
  `fig_*.csv`（测试集经验频率 + 训练集拟合 pmf）、`prior_*.csv`（训练集拟合的先验密度）；
  `kgrid_*.csv` 为训练集 K 选择曲线（由各数据文件夹的 results/k_selection_train.csv 复制）。
- `make_section5_figures.py`：Python（matplotlib，模拟实验作图风格）读取上述 CSV，
  生成 `Thesis/figures/plot16–19.png`（论文 Figure 7–10）。

### 6. 文档
- `mobility_doc.txt`、`deayala_doc.txt`：R 包原始帮助文档摘录。
- `数据集搜索报告.md`：26 个候选数据集的筛选结果汇总与推荐结论。
- `Mobility与deAyala数据集分析报告.md`：两个入选数据集的详细介绍与完整对比结果。

## 结论

- 最终入选数据集：**Mobility**（真实调查，n=8445，m=8）与 **deAyala**（教材基准，n=19601，m=5）。
- 按"一次性 70/30 训练/测试、只在训练集拟合与选 K、测试集最终检验"的协议，
  两者上 Bern-Bino 的测试集 logLik/L1/MSE 均最优（详见各数据文件夹 README 与论文第 5 节）。
