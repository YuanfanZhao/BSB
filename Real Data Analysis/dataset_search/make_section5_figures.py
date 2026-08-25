###############################################################################
# make_section5_figures.py  (plotting only; fits computed by export_fig_data.R)
# Produces Section-5 figures in the style of the paper's simulation experiments:
#   plot16.png  Mobility  : observed vs fitted pmf (5 models) + residuals
#   plot17.png  deAyala   : observed vs fitted pmf (5 models) + residuals
#   plot18.png  estimated prior densities pi(p) (Mobility & deAyala)
#   plot19.png  BIC vs K (Mobility & deAyala)
###############################################################################
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

BASE = r'E:\Thesis\BSB\Real Data Analysis\dataset_search'
FIG = r'E:\Thesis\BSB\Thesis\figures'

# ---------------- academic style (from simulation notebooks) ----------------
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

# ---------------- common style elements ----------------
MODELS = ['bin', 'bb', 'ln', 'km', 'bern']
LABELS = ['Binomial', 'Beta-Binomial', 'Logit-Normal', 'Kumaraswamy', 'Bern-Bino']
COLORS = ['#E41A1C', '#FF00F2', '#4DAF4A', '#984EA3', '#FF7F00']
LSTYLES = ['-', '--', '-.', ':', (0, (3, 1, 1, 1, 1, 1))]

def load(name):
    fig = pd.read_csv(os.path.join(BASE, 'fig_%s.csv' % name))
    prior = pd.read_csv(os.path.join(BASE, 'prior_%s.csv' % name))
    kg = pd.read_csv(os.path.join(BASE, 'kgrid_%s.csv' % name))
    return fig, prior, kg

# ---------------- plot16 / plot17 : fit + residuals ----------------
def plot_fit(name, fname, title, xlabel, m):
    fig, kg, _ = load(name)
    x = fig['x'].values; obs = fig['observed'].values
    fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    ax1.bar(x, obs, alpha=0.5, color='lightblue', label='Observed')
    for j, md in enumerate(MODELS):
        y = fig[md].values
        ax1.plot(x, y, color=COLORS[j], linestyle=LSTYLES[j], linewidth=2.5, label=LABELS[j])
        ax2.plot(x, y - obs, color=COLORS[j], linestyle=LSTYLES[j], linewidth=2.5, label=LABELS[j])
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel('Probability')
    ax1.set_title(title + ' (n=%d, m=%d)' % (int(fig['observed'].sum()), m))
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

# ---------------- plot18 : prior densities ----------------
def plot_prior(fname):
    fig2, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, name, title in zip(axes, ['mobility', 'deayala'],
                               ['Mobility', 'deAyala']):
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
    fig2.suptitle('Estimated prior densities of the success probability p', fontsize=15)
    fig2.tight_layout(rect=[0, 0, 1, 0.95])
    fig2.savefig(os.path.join(FIG, fname))
    plt.close(fig2)

# ---------------- plot19 : BIC vs K ----------------
def plot_bic(fname):
    fig2, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, name, title in zip(axes, ['mobility', 'deayala'], ['Mobility', 'deAyala']):
        _, _, kg = load(name)
        K = kg['K'].values; BIC = kg['BIC'].values
        ax.plot(K, BIC, color='#1f77b4', linewidth=2)
        ax.plot(K, BIC, 'o', color='#1f77b4', markersize=4)
        kb = int(K[np.argmin(BIC)])
        ax.plot(kb, BIC.min(), 'o', color='#d62728', markersize=9)
        ax.set_xlabel('K')
        ax.set_ylabel('BIC')
        ax.set_title(title + '  (best K = %d)' % kb)
        ax.grid(True, alpha=0.3)
    fig2.suptitle('BIC for selecting the Bernstein degree K', fontsize=15)
    fig2.tight_layout(rect=[0, 0, 1, 0.95])
    fig2.savefig(os.path.join(FIG, fname))
    plt.close(fig2)

plot_fit('mobility', 'plot16.png', 'Women\'s Mobility: fitted score distributions', 'Score (number of permitted activities)', 8)
plot_fit('deayala', 'plot17.png', 'deAyala: fitted score distributions', 'Score (number of correct answers)', 5)
plot_prior('plot18.png')
plot_bic('plot19.png')
print('figures written to', FIG)
