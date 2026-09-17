# -*- coding: utf-8 -*-
"""Figures for the MathExam14W regression report (corrected model)."""
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

BASE = os.path.dirname(os.path.abspath(__file__))
RES, FIG = os.path.join(BASE, "results"), os.path.join(BASE, "figures")
os.makedirs(FIG, exist_ok=True)
plt.style.use("default")
rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman"],
                 "mathtext.fontset": "stix", "font.size": 11,
                 "axes.labelsize": 12, "axes.titlesize": 12,
                 "xtick.labelsize": 10, "ytick.labelsize": 10,
                 "legend.fontsize": 10, "figure.dpi": 110, "savefig.dpi": 300,
                 "axes.linewidth": 1.1, "lines.linewidth": 2, "grid.alpha": 0.25})

ks = pd.read_csv(os.path.join(RES, "me2_k_selection.csv"))
K_bic = int(ks.loc[ks["BIC"].idxmin(), "K"])
K_cvll = int(ks.loc[ks["CV_mean_logLik"].idxmax(), "K"])
K_cvmse = int(ks.loc[ks["CV_mean_MSE"].idxmin(), "K"])

fig, axes = plt.subplots(1, 3, figsize=(13.2, 3.9))
axes[0].plot(ks["K"], ks["BIC"], "o-", color="#377EB8")
axes[0].axvline(K_bic, color="#E41A1C", ls="--", lw=1.4, label="BIC: K=%d" % K_bic)
axes[0].set_xlabel("Bernstein degree $K$"); axes[0].set_ylabel("BIC (training set)")
axes[0].set_title("(a) BIC"); axes[0].legend(); axes[0].grid(True, alpha=0.3)
axes[1].plot(ks["K"], ks["CV_mean_logLik"], "s-", color="#4DAF4A")
axes[1].axvline(K_cvll, color="#E41A1C", ls="--", lw=1.4, label="CV logLik: K=%d" % K_cvll)
axes[1].set_xlabel("Bernstein degree $K$")
axes[1].set_ylabel("Mean validation log-likelihood")
axes[1].set_title("(b) 5-fold repeated CV, log-likelihood")
axes[1].legend(); axes[1].grid(True, alpha=0.3)
axes[2].plot(ks["K"], ks["CV_mean_MSE"], "^-", color="#984EA3")
axes[2].axvline(K_cvmse, color="#E41A1C", ls="--", lw=1.4, label="CV MSE: K=%d" % K_cvmse)
axes[2].set_xlabel("Bernstein degree $K$"); axes[2].set_ylabel("Mean validation MSE")
axes[2].set_title("(c) 5-fold repeated CV, MSE")
axes[2].legend(); axes[2].grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_report_kselect.png")); plt.close(fig)
print("saved fig_report_kselect.png (BIC=%d, CV-ll=%d, CV-mse=%d)" % (K_bic, K_cvll, K_cvmse))

pr = pd.read_csv(os.path.join(RES, "me2_prior_profiles.csv"))
cal = pd.read_csv(os.path.join(RES, "me2_calibration.csv"))
fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.3))
axes[0].plot(pr["p"], pr["q10"], color="#E41A1C", label="10th pct. of the index")
axes[0].plot(pr["p"], pr["q50"], color="#377EB8", label="median index")
axes[0].plot(pr["p"], pr["q90"], color="#4DAF4A", label="90th pct. of the index")
axes[0].set_xlabel("$p$"); axes[0].set_ylabel("Density")
axes[0].set_title("(a) Fitted Bernstein prior $\\pi(p\\,|\\,z)$")
axes[0].legend(); axes[0].grid(True, alpha=0.3)
axes[1].plot([0, 1], [0, 1], color="grey", ls=":", lw=1.5, label="45 degrees")
axes[1].plot(cal["bern"], cal["obs"], "o-", color="#FF7F00", label="Bern-Bino")
axes[1].plot(cal["bbreg"], cal["obs"], "s--", color="#4DAF4A", label="Beta-Binomial")
axes[1].plot(cal["logistic"], cal["obs"], "^-.", color="#E41A1C", label="Logistic")
axes[1].set_xlabel("Mean predicted proportion in decile")
axes[1].set_ylabel("Observed proportion in decile")
axes[1].set_title("(b) Calibration on the test set")
axes[1].legend(); axes[1].grid(True, alpha=0.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_report_prior_cal.png")); plt.close(fig)
print("saved fig_report_prior_cal.png")
