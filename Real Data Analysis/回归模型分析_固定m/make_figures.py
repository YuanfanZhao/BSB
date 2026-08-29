###############################################################################
# make_figures.py --- 回归实例分析作图（风格与数值模拟实验一致）
# 数据来自 02 / 03 脚本：
#   results/fig_<label>.csv       测试集经验频率 + 各模型平均预测 pmf
#   results/prior_<label>.csv     Bern-Bino 平均先验密度 pi(p)
#   results/kgrid_bic_<label>.csv 训练集 BIC 随 K 变化
#   results/kgrid_cv_<label>.csv  训练集 5 折 CV 对数似然随 K 变化
# 输出 figures/ 下的 PNG（300 dpi）。
###############################################################################
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

BASE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)

def set_academic_style():
    plt.style.use("default")
    rcParams["font.family"] = "serif"
    rcParams["font.serif"] = ["Times New Roman"]
    rcParams["mathtext.fontset"] = "stix"
    rcParams["font.size"] = 12
    rcParams["axes.labelsize"] = 14
    rcParams["axes.titlesize"] = 15
    rcParams["xtick.labelsize"] = 11
    rcParams["ytick.labelsize"] = 11
    rcParams["legend.fontsize"] = 11
    rcParams["figure.figsize"] = (10, 6)
    rcParams["figure.dpi"] = 100
    rcParams["savefig.dpi"] = 300
    rcParams["axes.linewidth"] = 1.2
    rcParams["lines.linewidth"] = 2
    rcParams["grid.linewidth"] = 0.8
    rcParams["grid.alpha"] = 0.25
set_academic_style()

MODELS = ["logistic", "bbreg", "bern"]
LABELS = ["Logistic", "Beta-Binomial", "Bern-Bino"]
COLORS = ["#E41A1C", "#4DAF4A", "#FF7F00"]
LSTYLES = ["-", "--", "-."]
LINEW = 2.5

DATASETS = [
    ("MathExam14W", "Mathematics 101 Exam (MathExam14W)", "Number of solved items $x$", 13),
    ("PISA2000Read", "PISA 2000 Reading (PISA2000Read)", "Number of correct reading items $x$", 23),
]

# ---------------- (1) 测试集经验频率 vs 各模型平均预测 pmf + 残差 ----------------
def plot_fit(label, title, xlabel, m):
    fig = pd.read_csv(os.path.join(RES, "fig_%s.csv" % label))
    x = fig["x"].values
    obs = fig["observed"].values
    fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True)
    ax1.bar(x, obs, alpha=0.45, color="lightsteelblue", edgecolor="steelblue",
            linewidth=0.6, label="Observed (test set)")
    for j, md in enumerate(MODELS):
        y = fig[md].values
        ax1.plot(x, y, color=COLORS[j], linestyle=LSTYLES[j], linewidth=LINEW,
                 label=LABELS[j])
        ax2.plot(x, y - obs, color=COLORS[j], linestyle=LSTYLES[j], linewidth=LINEW,
                 label=LABELS[j])
    ax1.set_ylabel("Probability")
    ax1.set_title(title + ": fitted pmf vs test-set frequencies")
    ax1.legend(loc="upper right", ncol=2)
    ax1.grid(True, alpha=0.3)
    ax2.axhline(0, color="black", alpha=0.4)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("Residual (fitted $-$ observed)")
    ax2.legend(loc="upper right", ncol=3, fontsize=10)
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(os.path.join(FIG, "fig_pmf_%s.png" % label))
    plt.close(fig2)
    print("saved fig_pmf_%s.png" % label)

# ---------------- (2) K 选择曲线（BIC 与 CV） ----------------
def plot_kgrid(label, title):
    bic = pd.read_csv(os.path.join(RES, "kgrid_bic_%s.csv" % label))
    cv = pd.read_csv(os.path.join(RES, "kgrid_cv_%s.csv" % label))
    Kb = int(bic.loc[bic["BIC"].idxmin(), "K"])
    Kc = int(cv.loc[cv["cv_mean_ll"].idxmax(), "K"])
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
    ax1.plot(bic["K"], bic["BIC"], "o-", color="#377EB8", linewidth=2)
    ax1.axvline(Kb, color="#E41A1C", linestyle="--", linewidth=1.5,
                label="BIC-optimal K = %d" % Kb)
    ax1.set_xlabel("$K$"); ax1.set_ylabel("Training-set BIC")
    ax1.set_title(title + ": BIC")
    ax1.legend(); ax1.grid(True, alpha=0.3)
    ax2.plot(cv["K"], cv["cv_mean_ll"], "s-", color="#4DAF4A", linewidth=2)
    ax2.axvline(Kc, color="#E41A1C", linestyle="--", linewidth=1.5,
                label="CV-optimal K = %d" % Kc)
    ax2.set_xlabel("$K$"); ax2.set_ylabel("Mean validation log-likelihood")
    ax2.set_title(title + ": 5-fold CV")
    ax2.legend(); ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(os.path.join(FIG, "fig_kgrid_%s.png" % label))
    plt.close(fig2)
    print("saved fig_kgrid_%s.png" % label)

# ---------------- (3) Bern-Bino 平均先验密度 pi(p) ----------------
def plot_prior(label, title):
    pr = pd.read_csv(os.path.join(RES, "prior_%s.csv" % label))
    fig2, ax = plt.subplots(figsize=(8, 5))
    ax.plot(pr["p"], pr["mean_prior"], color="#FF7F00", linewidth=2.5)
    ax.fill_between(pr["p"], 0, pr["mean_prior"], color="#FF7F00", alpha=0.2)
    ax.set_xlabel("$p$"); ax.set_ylabel("Density")
    ax.set_title(title + ": average implied prior density $\\pi(p)$")
    ax.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(os.path.join(FIG, "fig_prior_%s.png" % label))
    plt.close(fig2)
    print("saved fig_prior_%s.png" % label)

for label, title, xlabel, m in DATASETS:
    plot_fit(label, title, xlabel, m)
    plot_kgrid(label, title)
    plot_prior(label, title)
print("all figures done.")
