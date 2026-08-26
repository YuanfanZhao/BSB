###############################################################################
# make_section5_figures.py  (plotting only; fits/data computed by export_fig_data.R)
# New protocol: models are fit on the 70% training set; the figures show
#   - plot16/17: held-out TEST-set observed frequencies vs training-set fitted
#     pmfs of the five models (plus residuals)
#   - plot18: estimated prior densities pi(p) from training-set fits
#   - plot19: K-selection curves on the training set (BIC, inner-CV mean
#     validation log-likelihood / MSE), marking the final K
# Style follows the paper's simulation experiments (matplotlib).
###############################################################################
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

BASE = r'E:\Thesis\BSB\Real Data Analysis\dataset_search'
FIG = r'E:\Thesis\BSB\Thesis\figures'

def set_academic_style():
    plt.style.use('default')
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Times New Roman']
    rcParams['mathtext.fontset'] = 'stix'
    rcParams['font.size'] = 12
    rcParams['axes.labelsize'] = 14
    rcParams['axes.titlesize'] = 16
    rcParams['xtick.labelsize'] = 12
    rcParams['ytick.labelsize'] = 12
    rcParams['legend.fontsize'] = 12
    rcParams['figure.figsize'] = (10, 6)
    rcParams['figure.dpi'] = 100
    rcParams['savefig.dpi'] = 300
    rcParams['axes.linewidth'] = 1.2
    rcParams['lines.linewidth'] = 2
    rcParams['grid.linewidth'] = 0.8
    rcParams['grid.alpha'] = 0.25
set_academic_style()

MODELS = ['bin', 'bb', 'ln', 'km', 'bern']
LABELS = ['Binomial', 'Beta-Binomial', 'Logit-Normal', 'Kumaraswamy', 'Bern-Bino']
COLORS = ['#E41A1C', '#FF00F2', '#4DAF4A', '#984EA3', '#FF7F00']
LSTYLES = ['-', '--', '-.', ':', (0, (3, 1, 1, 1, 1, 1))]
FINAL_K = {'mobility': 25, 'deayala': 30}

def load(name):
    fig = pd.read_csv(os.path.join(BASE, 'fig_%s.csv' % name))
    prior = pd.read_csv(os.path.join(BASE, 'prior_%s.csv' % name))
    kg = pd.read_csv(os.path.join(BASE, 'kgrid_%s.csv' % name))
    return fig, prior, kg

# ---------------- plot16 / plot17 : test observed vs train-fitted pmfs ----
def plot_fit(name, fname, title, xlabel, m):
    fig, _, _ = load(name)
    x = fig['x'].values; obs = fig['observed'].values
    fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    ax1.bar(x, obs, alpha=0.5, color='lightblue', label='Observed (test set)')
    for j, md in enumerate(MODELS):
        y = fig[md].values
        ax1.plot(x, y, color=COLORS[j], linestyle=LSTYLES[j], linewidth=2.5, label=LABELS[j])
        ax2.plot(x, y - obs, color=COLORS[j], linestyle=LSTYLES[j], linewidth=2.5, label=LABELS[j])
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel('Probability')
    ax1.set_title(title + ' (held-out test set)')
    ax1.legend(loc='upper right', ncol=2)
    ax1.grid(True, alpha=0.3)
    ax2.axhline(0, color='black', alpha=0.4)
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel('Residual (Fitted - Observed)')
    ax2.legend(loc='upper right', ncol=2, fontsize=10)
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(os.path.join(FIG, fname))
    plt.close(fig2)

# ---------------- plot18 : prior densities (train fits) -------------------
def plot_prior(fname):
    fig2, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, name, title in zip(axes, ['mobility', 'deayala'], ['Mobility', 'deAyala']):
        _, prior, _ = load(name)
        p = prior['p'].values
        ax.plot(p, prior['bern'].values, color=COLORS[4], linestyle='-', linewidth=2.5, label='Bern-Bino')
        ax.plot(p, prior['bb'].values, color=COLORS[1], linestyle='--', linewidth=2, label='Beta')
        ax.plot(p, prior['ln'].values, color=COLORS[2], linestyle='-.', linewidth=2, label='Logit-Normal')
        ax.plot(p, prior['km'].values, color=COLORS[3], linestyle=':', linewidth=2, label='Kumaraswamy')
        ax.set_xlabel('p')
        ax.set_ylabel('Density')
        ax.set_title(title)
        ax.legend(loc='upper center', fontsize=10)
        ax.grid(True, alpha=0.3)
    fig2.suptitle('Estimated prior densities of the success probability p (training-set fits)', fontsize=15)
    fig2.tight_layout(rect=[0, 0, 1, 0.95])
    fig2.savefig(os.path.join(FIG, fname))
    plt.close(fig2)

# ---------------- plot19 : K selection on the training set -----------------
def plot_bic(fname):
    fig2, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, name, title in zip(axes, ['mobility', 'deayala'], ['Mobility', 'deAyala']):
        _, _, kg = load(name)
        K = kg['K'].values
        kb = int(FINAL_K[name])
        # 左图：训练集 BIC vs K
        ax.plot(K, kg['BIC_train'].values, color='#1f77b4', linewidth=2)
        ax.plot(K, kg['BIC_train'].values, 'o', color='#1f77b4', markersize=4)
        ax.plot(kb, kg['BIC_train'].values[kg['K'] == kb][0], 'o', color='#d62728', markersize=9)
        ax.set_xlabel('K'); ax.set_ylabel('BIC (training set)')
        ax.set_title(title + '  (final K = %d)' % kb)
        ax.grid(True, alpha=0.3)
        # 右图：内层 CV 的平均验证似然与 MSE（双纵轴）
        ax2 = ax.twinx()
        l1, = ax.plot(K, kg['CV_mean_ll'].values, color='#2ca02c', linewidth=2, label='CV mean logLik')
        ax.plot(K, kg['CV_mean_ll'].values, 'o', color='#2ca02c', markersize=4)
        l2, = ax2.plot(K, kg['CV_mean_mse'].values, color='#984EA3', linestyle='--', linewidth=2, label='CV mean MSE')
        ax2.plot(kb, kg['CV_mean_mse'].values[kg['K'] == kb][0], 'o', color='#d62728', markersize=9)
        ax.set_xlabel('K'); ax.set_ylabel('CV mean logLik')
        ax2.set_ylabel('CV mean MSE')
        ax.set_title(title + '  (final K = %d)' % kb)
        ax.grid(True, alpha=0.3)
        ax.legend(handles=[l1, l2], loc='center right', fontsize=9)
    fig2.suptitle('Bernstein degree K selection on the training set (BIC and inner 5-fold CV)', fontsize=15)
    fig2.tight_layout(rect=[0, 0, 1, 0.95])
    fig2.savefig(os.path.join(FIG, fname))
    plt.close(fig2)

plot_fit('mobility', 'plot16.png', 'Women\'s Mobility: fitted score distributions', 'Score (number of permitted activities)', 8)
plot_fit('deayala', 'plot17.png', 'deAyala: fitted score distributions', 'Score (number of correct answers)', 5)
plot_prior('plot18.png')
plot_bic('plot19.png')
print('figures written to', FIG)
