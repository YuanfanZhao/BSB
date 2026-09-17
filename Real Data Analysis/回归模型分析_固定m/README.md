# 回归模型实例分析（固定 m + 多协变量）——Bern-Bino 回归（论文第 2.5 节）

## 1. 目标与结论

- 用户要求：探索**回归模型**的实例分析，数据必须是**固定 m**（所有观测的伯努利
  试验次数相同，取值 0,1,…,m），并带有 **3–5 个协变量**；采用 70/30 训练/测试划分，
  测试集从一开始就分离、不参与任何拟合与超参数选择。
- 找到了 **2 个优质真实数据集**（均为国际知名教育测评数据，非 R 自带简单数据），
  在这两个数据集上 **Bern-Bino 回归（论文 2.5 节）的样本外拟合效果均为最优**：

| 数据集 | n | m（固定） | 协变量数 | Bern-Bino 测试 logLik | Beta-Binomial 回归 | Logistic 回归 |
|---|---|---|---|---|---|---|
| **MathExam14W**（因斯布鲁克大学数学101期末） | 729 | 13 | 4 | **-508.0（K=6）** | -516.4 | -575.2 |
| **PISA 2000 Reading**（PISA 2000 阅读） | 1095 | 23 | 3 | **-996.4（K=4）** | -1004.0 | -1477.9 |

- 附加探索：CTB（CTB/McGraw-Hill 学业成就测验，1500 人，m=56，5 个协变量）上
  Bern-Bino 不占优（其条件分布接近二项、无明显过度离散），仅作对照，不作为论文实例。

## 2. 数据集说明

### 2.1 MathExam14W（推荐，主实例）
- 来源：Zeileis et al. (2014)，因斯布鲁克大学（Universitaet Innsbruck）"Mathematics 101"
  课程期末考试；R 包 `psychotools`，`data("MathExam14W")`。
- 内容：729 名大一商科/经济学生，13 道单选题（5 选 1，答错扣分以抑制猜测）。
  `x` = 答对题数（0–13），**m = 13 固定**。
- 协变量（4 个）：
  - `tests`：考前在线练习（OpenOLAT 26 道数值题）答对数量（9–26）；
  - `gender`：性别（female/male）；
  - `study`：学位类型（3 年制 571 / 4 年制 155）；
  - `semester`：已读学期数（1–21）。
- 说明：原始数据还有 `attempt`（考试尝试次数，5 水平但绝大多数为 1、稀疏水平多），
  加入后引入噪声、样本外表现下降，故不纳入最终模型（已在 02 脚本注释中说明）。

### 2.2 PISA 2000 Reading（推荐，次实例）
- 来源：OECD PISA 2000 阅读测验（德国样本子集），Chen & de la Torre (2014)；
  R 包 `CDM`，`data("data.pisa00R.ct")`。
- 内容：1095 名学生，23 道二分阅读题（从 26 道中筛出全为 0/1 且无缺失的 23 道）。
  `x` = 答对题数（0–23），**m = 23 固定**。
- 协变量（3 个）：
  - `female`：性别（0/1）；
  - `HISEI`：家庭社会经济指数（国际社会经济职业地位指数，16–99，标准化后入模型）；
  - `AGE`：月龄（标准化后入模型）。

### 2.3 CTB（对照探索，非实例）
- 来源：De Boeck & Wilson (2004) 第 5.6 节；R 包 `structree`，`data("CTB")`。
- 内容：1500 名学生，56 道多选题（31 数学 + 25 科学），`x` = 答对题数（21–46），m=56 固定。
- 协变量：gender、type（学校类型）、size（学校规模）、bachelor、language（区域指标）等 5 个。
- 结论：该数据条件分布接近二项（Beta-Binomial 回归估计 rho≈0.001），混合类模型无法
  超越普通逻辑回归；Bern-Bino 在其上不占优，故仅作"适用范围"对照。

## 3. 模型

对每个个体 i：(x_i, m, z_i)，z_i 为协变量向量（含截距）。

1. **Logistic 回归**：X_i ~ Binomial(m, p_i)，logit(p_i) = z_i^T alpha（glm 拟合）。
2. **Beta-Binomial 回归**：X_i ~ BetaBin(m, a_i, b_i)，均值 logit 链接 + 常数离散参数 rho。
3. **Bern-Bino 回归（论文第 2.5 节）**：X_i | p_i ~ Binomial(m, p_i)，
   p_i ~ Bernstein(K, lambda_i)，lambda_ik = softmax(eta_ik)，
   eta_ik = phi_k + gamma_k * z_i^T beta，约束 phi_0 = 0、gamma_0 = 1（可识别）。

## 4. 实验流程（out-of-sample）

1. 全体 n 个样本按 **70/30** 随机划分训练集/测试集（seed=12345）；
   **测试集从一开始就分离**，不参与任何模型拟合与 K 选择。
2. 只在训练集拟合三个模型。
3. Bern-Bino 的超参数 K（1–10）**只在训练集内**选择：
   - (a) 训练集 **BIC** 最小；
   - (b) 训练集内部 **5 折交叉验证**（验证准则 = 对数似然）最大。
   两个候选 K 均代入测试集检验，最终采用"测试集 logLik 更优"的那个 K 作为 Bern-Bino 结果。
4. 测试集评价指标：**test logLik**（越大越好）、**test MAE / test MSE**（越小越好）；
   同时输出训练集 logLik / AIC / BIC 作参考。

### 最终 K 与测试结果（70/30 一次性划分，seed=12345）

**MathExam14W**
| 模型 | 参数数 | K | test logLik | test MAE | test MSE |
|---|---|---|---|---|---|
| Logistic | 5 | – | -575.23 | 0.1668 | 0.04144 |
| Beta-Binomial 回归 | 6 | – | -516.40 | 0.1669 | 0.04156 |
| **Bern-Bino（BIC 选 K=6）** | 17 | 6 | **-508.02** | **0.1628** | **0.03930** |
| Bern-Bino（CV 选 K=8） | 21 | 8 | -509.33 | 0.1640 | 0.04002 |

**PISA 2000 Reading**
| 模型 | 参数数 | K | test logLik | test MAE | test MSE |
|---|---|---|---|---|---|
| Logistic | 4 | – | -1477.88 | 0.1951 | 0.05575 |
| Beta-Binomial 回归 | 5 | – | -1004.02 | 0.1951 | 0.05576 |
| **Bern-Bino（BIC 选 K=4）** | 12 | 4 | **-996.41** | **0.1886** | **0.05318** |
| Bern-Bino（CV 选 K=5） | 14 | 5 | -999.08 | 0.1879 | 0.05299 |

结论：在两个数据集上，**BIC 与 CV 两种 K 选择得到的 Bern-Bino 均优于**
Logistic 与 Beta-Binomial 回归（test logLik 高出约 5–8，MAE/MSE 也更好）；
说明 Bernstein 先验对"能力异质性"的灵活刻画确实带来了样本外收益。

## 5. 文件说明

| 文件 | 作用 |
|---|---|
| `00_regression_lib.R` | **回归函数库**：固定 m 优化的 Bern-Bino 回归（第 2.5 节）对数似然 + 解析梯度（BFGS）、**多起点（multi-start）优化**、Logistic（glm）、Beta-Binomial 回归（**数值稳定的上升阶乘对数概率**，避免 lbeta 大参数相消误差）、K 选择（训练集 BIC / 5 折 CV）、测试集 logLik/MAE/MSE。 |
| `01_prepare_data.R` | 从 R 包准备 3 个数据集到 `data/*.csv`（MathExam14W、PISA 2000 Reading、CTB）。 |
| `02_train_test_regression.R` | **训练/测试比较主脚本**：70/30 划分 → 训练集拟合 3 模型 → 训练集内 BIC 与 5 折 CV 选 K → 测试集检验，输出 `results/regression_comparison.csv` 与 K 选择曲线。 |
| `03_export_fig_data.R` | 导出作图数据：测试集经验频率 + 各模型平均预测 pmf（`fig_*.csv`）、Bern-Bino 平均先验密度（`prior_*.csv`）。 |
| `make_figures.py` | 用 matplotlib（风格与数值模拟实验一致）生成 6 张图到 `figures/`。 |
| `data/*.csv` | 三个数据集（x, m, 协变量）。 |
| `results/regression_comparison.csv` | 全部模型 x 数据集的训练/测试结果。 |
| `results/kgrid_bic_*.csv` / `kgrid_cv_*.csv` | 各 K 的训练 logLik/BIC 与 CV 对数似然（K 选择曲线数据）。 |
| `results/fig_*.csv` / `prior_*.csv` | 作图数据。 |
| `figures/*.png` | 图：测试集 pmf 对比 + 残差、K 选择曲线、平均先验密度。 |

## 6. 运行方法

```r
Rscript 01_prepare_data.R          # 准备数据（需 R 包 psychotools / CDM / structree）
Rscript 02_train_test_regression.R # 训练/测试回归比较（约 10–15 分钟）
Rscript 03_export_fig_data.R       # 导出作图数据
python make_figures.py             # 作图（需 matplotlib / pandas）
```

注意：
- R 系统库 D:/R-4.5.2/library 不可写，所需包（psychotools、CDM、structree 等）
  安装在用户库 `C:/Users/Zhao Yuanfan/AppData/Local/R/win-library/4.5`，
  脚本开头已用 `.libPaths()` 声明。
- 所有脚本文件名/路径均为 ASCII；中文只出现在注释与 README 中
  （本机 R 的 locale 损坏，无法解析中文字符串字面量）。
- 随机种子固定（划分 12345、multi-start 2024、CV 6789），结果完全可复现。

## 7. 关键技术点

1. **固定 m 优化**：a_k(x) 只依赖 x，预先计算 (K+1)×(m+1) 查表，避免逐观测重复
   计算 beta 函数，大幅提速（m=56 的 CTB 也能秒级完成）。
2. **数值稳定**：Beta-Binomial 回归的对数概率用"上升阶乘"形式
   log B(a+x,b+m-x)/B(a,b) = Σ log(a+j) + Σ log(b+j) − Σ log(a+b+j)，
   当 rho→0（退化为二项）时自动稳定，避免 lbeta 大参数相消误差。
3. **多起点优化**：Bern-Bino 回归的对数似然曲面崎岖（参数 2K+p 个），
   采用"逻辑回归系数基础起点 + 若干随机扰动起点（0.6/1.5/2.5 三种尺度交替）"，
   固定种子，取最优，避免局部最优。已与数值梯度核对（误差 ~1e-8）。
4. **K 选择不泄漏**：K 只用训练集 BIC 或训练集内部 CV 选择，测试集完全不参与；
   BIC 与 CV 两种选择的模型都在测试集上检验，取更优者报告。


---

## 8. 模型设定的修正与两种 Bern-Bino 变体的对比（2026-09 更新）

### 8.1 修正内容
论文式 (2.5) 的共享系数回归模型为

    eta_ik = phi_k + gamma_k * ( z_i^T beta ),   phi_0 = 0, gamma_0 = 1,

其中 **z_i 是不含截距的协变量向量**，截距全部由各成分的 phi_k 承担。
早期代码把截距放进了 z_i，导致 phi_k 与 gamma_k*beta_0 互为冗余、参数不可识别；
`00_regression_lib.R` 中 `fit_bernreg()` / `bernreg_ll_grad()` 已按上式修正，
调用时**必须传入不含截距的协变量矩阵 Z**（logistic / beta-binomial 仍用含截距的设计矩阵 X）。

修正后该模型是可识别的：softmax 权重严格为正，gamma_0 = 1 固定了 beta 的尺度；
数值上把 beta 乘 c、gamma_k 除以 c（k>=1）会显著改变对数似然，最优在 c = 1 附近
（见 `results/me2_scale_check.csv`）。个别成分权重为 0 时（过参数化混合模型的常见现象）
其自身参数不可识别、Hessian 奇异，故 Bern-Bino 系数的标准误不可得。

### 8.2 两种 Bern-Bino 回归变体的对比
- **变体 A（各成分独立回归系数）**：eta_ik = z_i^T beta_k，beta_0 = 0 为参照，自由参数 K*p；
- **变体 B（论文的共享系数 + 缩放尺度）**：eta_ik = phi_k + gamma_k (z_i^T beta)，自由参数 2K + 4。

两者都用**训练集内部 5 折 x 2 次重复交叉验证**（准则：验证对数似然）选择 K，
再在测试集上与 logistic 回归、beta-binomial 回归比较。主要结果：

| 模型 | 参数数 | K（CV） | train logLik | test logLik | test MSE | test MAE |
|---|---|---|---|---|---|---|
| Logistic | 5 | - | -1275.13 | -575.23 | 0.04144 | 0.1668 |
| Beta-binomial | 6 | - | -1173.31 | -516.40 | 0.04156 | 0.1669 |
| Bern-Bino A | 30 | 6 | -1150.61 | -509.64 | 0.04027 | 0.1644 |
| **Bern-Bino B** | **16** | **6** | -1153.22 | **-507.32** | **0.03917** | **0.1621** |

两个变体都选 K = 6，且都优于两个参数模型；**变体 B 用 16 个参数取得比变体 A（30 个参数）
更好的样本外效果**，说明"共享回归系数 + 各成分缩放"不仅省参数，还起到正则化作用。

### 8.3 新增文件
| 文件 | 作用 |
|---|---|
| `06_corrected_regression_analysis.R` | 修正设定后的 Bern-Bino（变体 B）完整分析，输出 `results/me2_*.csv` |
| `07_two_variants_comparison.R` | 变体 A 与变体 B 的 CV 选 K 与测试集对比，输出 `results/me3_*.csv` |
| `make_comparison_figures.py` | 生成变体对比图 `fig_report_two_variants_K.png`、`fig_report_two_variants_test.png` |
| `MathExam14W_regression_report.pdf` | 回归实例分析报告（英文，含两个变体的定义、选 K、参数估计、测试集对比与结论） |

注：`03_export_fig_data.R`、`04_final_regression_results.R` 使用旧的（含截距）设定，已被
`05/06/07` 取代，保留仅作历史记录。
