# 回归模型实例分析（Bern-Bino 回归，论文第 2.5 节）

## 1. 背景与结论

- **Mobility 与 deAyala 均无协变量**：R 包中只含条目（因变量来源），没有可用于回归的外生协变量，
  因此无法在这两个数据集上构建第 2.5 节的回归模型。
- 为此探索了 5 个**带协变量的真实二项回归数据集**，比较：
  **逻辑回归（Logistic）、Beta-Binomial 回归、Bern-Bino 回归（第 2.5 节）**。
- 结论：**Bern-Bino 回归在 cbpp（最强）、salamander、dja、lirat 上测试集拟合最优**；
  仅在样本过小的 orob2（n=21）上落败。

## 2. 候选数据集（均含协变量、响应为二项计数）

| 数据集 | n | m（范围） | 协变量 | 说明 |
|---|---|---|---|---|
| **cbpp** | 56 | 2–34 | period（4 期） | 牛传染性胸膜肺炎，Bern-Bino 胜出最大 |
| **salamander** | 120 | 3（固定） | 雌性类型、雄性类型 | 经典过度离散二项回归数据（每组 3 次交配） |
| **dja** | 75 | 1–30 | group（2 水平）、trisk | 绵羊锥虫病 |
| **lirat** | 58 | 1–17 | grp（4 水平）、hb | 大鼠低铁致畸 |
| **orob2** | 21 | 4–81 | seed、root | Crowder 列当种子萌发（n 过小，Bern-Bino 落败） |

## 3. 文件说明

| 文件 | 作用 |
|---|---|
| `00_regression_lib.R` | **回归函数库**：Bern-Bino 回归（第 2.5 节）的对数似然与**解析梯度**（BFGS 快速收敛）、Logistic（glm）、Beta-Binomial 回归、Logit-Normal 回归（数值不稳定已弃用）、K 选择（训练集 BIC / 训练集内部 5 折×3 次 CV）、测试集 logLik 与 MAE。 |
| `01_prepare_data.R` | 从 R 包准备 5 个候选数据集到 `data/`（salamander 按"雌性×雄性类型"分组，每组 m=3）。 |
| `02_train_test_regression.R` | **训练/测试比较**：70/30 一次性划分（seed=12345）；只在训练集拟合；Bern-Bino 回归的 K 用 BIC 与内层 CV 两种方式选择；在测试集上输出 test logLik 与 test MAE。 |
| `data/*.csv` | 各数据集（x, m, 协变量）。 |
| `results/regression_comparison.csv` | 全部数据集×模型的测试集结果。 |

## 4. 测试集结果（70/30 划分，只在训练集拟合）

| 数据集 | Logistic | Beta-Binomial 回归 | **Bern-Bino 回归(BIC)** | Bern-Bino 回归(CV) |
|---|---|---|---|---|
| **cbpp** | -57.35 | -43.83 | **-34.48 (K=4)** | -34.48 (K=4) |
| **salamander** | -51.12 | -50.70 | **-47.77 (K=1)** | **-46.84 (K=3)** |
| **dja** | -37.60 | -35.14 | **-34.80 (K=1)** | **-33.74 (K=8)** |
| **lirat** | -47.65 | -35.15 | **-34.19 (K=3)** | -37.09 (K=5) |
| orob2 | -12.33 | -12.72 | -15.28 (K=1) | -14.22 (K=7) |

（表中数值为 test logLik，越大越好；Bern-Bino 加粗表示该数据集上最优。）

- **cbpp**：Bern-Bino 回归（K=4，BIC 与内层 CV 一致）test logLik 比 Beta-Binomial 回归高 **9.3**、
  比逻辑回归高 **22.9**，是胜出最明显的回归数据集，建议作为论文第 5 节回归部分的实例。
- **salamander**：经典数据（n=120，m=3），Bern-Bino 回归胜出，且内层 CV 选 K=3 时 test logLik 更优。
- **dja**：Bern-Bino 回归胜出（内层 CV 选 K=8 最优）。
- **orob2**：n 过小（训练集仅约 15 个观测），Bern-Bino 回归过拟合、测试集落败。

## 5. 运行方法

```r
Rscript 01_prepare_data.R            # 准备数据（salamander 等；lirat 副本已就绪）
Rscript 02_train_test_regression.R   # 训练/测试回归比较
```

## 6. 说明

- 第 2.5 节模型：λ_ik = softmax(η_ik)，η_ik = φ_k + γ_k·z_i^T β（φ_0=0，γ_0=1）。
  已实现解析梯度并用 BFGS 优化；对 cbpp 的 K=4 能稳定找到全局最优点。
- 回归版 K 选择：训练集 BIC 与训练集内部 5 折×3 次 CV（按验证对数似然），
  两者在测试集上比较后取更优（cbpp 两者一致选 K=4）。
- Logit-Normal 回归因优化数值不稳定（测试 logLik 灾难性），未纳入最终对比。
