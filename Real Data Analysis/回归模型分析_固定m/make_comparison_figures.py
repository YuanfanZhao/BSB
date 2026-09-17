# -*- coding: utf-8 -*-
"""Figures for the two-variant Bern-Bino regression comparison."""
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

BASE = os.path.dirname(os.path.abspath(__file__))
RES, FIG = os.path.join(BASE, "results"), os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)
plt.style.use("default")
rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman"],
                 "mathtext.fontset": "stix", "font.size": 11, "axes.labelsize": 12,
                 "axes.titlesize": 12, "xtick.labelsize": 10, "ytick.labelsize": 10,
                 "legend.fontsize": 10, "figure.dpi": 110, "savefig.dpi": 300,
                 "axes.linewidth": 1.1, "lines.linewidth": 2, "grid.alpha": 0.25})

ks = pd.read_csv(os.path.join(RES, "me3_k_selection.csv"))
KA = int(ks.loc[ks["A_CV_logLik"].idxmax(), "K"])
KB = int(ks.loc[ks["B_CV_logLik"].idxmax(), "K"])

fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.9))
ax = axes[0]
ax.plot(ks["K"], ks["A_CV_logLik"], "o-", color="#E41A1C", label="Variant A (K=%d)" % KA)
ax.plot(ks["K"], ks["B_CV_logLik"], "s-", color="#377EB8", label="Variant B (K=%d)" % KB)
ax.axvline(KA, color="#E41A1C", ls=":", lw=1.3)
ax.axvline(KB, color="#377EB8", ls="--", lw=1.3)
ax.set_xlabel("Bernstein degree $K$"); ax.set_ylabel("Mean validation log-likelihood")
ax.set_title("(a) CV with validation log-likelihood"); ax.legend(); ax.grid(True, alpha=0.3)
ax = axes[1]
ax.plot(ks["K"], ks["A_CV_MSE"], "o-", color="#E41A1C", label="Variant A")
ax.plot(ks["K"], ks["B_CV_MSE"], "s-", color="#377EB8", label="Variant B")
ax.axvline(KA, color="#E41A1C", ls=":", lw=1.3); ax.axvline(KB, color="#377EB8", ls="--", lw=1.3)
ax.set_xlabel("Bernstein degree $K$"); ax.set_ylabel("Mean validation MSE")
ax.set_title("(b) CV with validation MSE"); ax.legend(); ax.grid(True, alpha=0.3)
ax = axes[2]
ax.plot(ks["K"], ks["A_BIC"], "o-", color="#E41A1C", label="Variant A")
ax.plot(ks["K"], ks["B_BIC"], "s-", color="#377EB8", label="Variant B")
ax.set_xlabel("Bernstein degree $K$"); ax.set_ylabel("BIC (training set)")
ax.set_title("(c) BIC"); ax.legend(); ax.grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_report_two_variants_K.png")); plt.close(fig)
print("saved fig_report_two_variants_K.png (K_A=%d, K_B=%d)" % (KA, KB))

te = pd.read_csv(os.path.join(RES, "me3_test_metrics.csv"))
labels = {"Logistic": "Logistic", "BetaBinom": "Beta-\nbinomial",
          "BernA": "Bern-Bino\nA (own $\\beta_k$)", "BernB": "Bern-Bino\nB (shared $\\beta$)"}
x = list(range(len(te))); names = [labels[m] for m in te["Model"]]
cols = ["#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00"]
fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
axes[0].bar(x, -te["test_logLik"], color=cols); axes[0].set_ylabel("$-$\u2009test log-likelihood")
axes[0].set_title("(a) Test log-likelihood (lower is better)")
axes[1].bar(x, te["test_MAE"], color=cols); axes[1].set_ylabel("test MAE")
axes[1].set_title("(b) Test MAE (lower is better)")
axes[2].bar(x, te["test_MSE"], color=cols); axes[2].set_ylabel("test MSE")
axes[2].set_title("(c) Test MSE (lower is better)")
for ax in axes:
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=8.5); ax.grid(True, axis="y", alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_report_two_variants_test.png")); plt.close(fig)
print("saved fig_report_two_variants_test.png")
