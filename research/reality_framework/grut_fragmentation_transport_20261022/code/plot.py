#!/usr/bin/env python3
"""Standalone scientific figure; plotted clumps are not identified objects."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
fig, ax = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
for n, color in ((64, '#2c7fb8'), (96, '#f28e2b'), (128, '#5945a5')):
    with (ROOT/f'data/n{n}_seed17_dt0.02.csv').open() as f:
        rows = list(csv.DictReader(f))
    t = np.array([float(r['tau']) for r in rows])
    rho = np.array([float(r['peak_density']) for r in rows])
    ax[0].semilogy(t, rho, label=f'{n}³', color=color)
    ax[0].plot(t[-1], rho[-1], 'x', color=color)
ax[0].set(xlabel=r'Time $\tau=\Gamma_0 t$', ylabel=r'Peak charge density / $n_0$',
          title='Collapse reaches the grid limit')
ax[0].legend(frameon=False)
array = np.load(ROOT/'data/n128_seed17_dt0.02.npz')['tau_6_density_midplane']
im = ax[1].imshow(array.T, origin='lower', extent=(0, 4*np.pi, 0, 4*np.pi), cmap='viridis', vmin=float(array.min()), vmax=float(array.max()))
ax[1].set(xlabel=r'$k_0 x$', ylabel=r'$k_0 y$', title=r'Resolved midplane at $\tau=6$')
fig.colorbar(im, ax=ax[1], shrink=.8, label=r'$n/n_0$')
curves = np.load(ROOT/'TRANSPORT_RESULTS.npz')
for gamma, color in ((1, '#2c7fb8'), (10, '#f28e2b'), (1000, '#5945a5')):
    z, _, weighted = curves[f'gamma_{gamma:g}']
    cumulative = cumulative_trapezoid(weighted, -np.log1p(z), initial=0.)
    mask = (z >= 1000) & (z <= 1e7)
    ax[2].loglog(z[mask], np.maximum(cumulative[mask], 1e-12), label=rf'$\gamma_i={gamma}$', color=color)
ax[2].axhline(6e-5, color='#333333', linestyle=':', label='Historical heat screen')
ax[2].set(xlim=(1e7, 1e3), ylim=(1e-9, .1), xlabel='Redshift',
          ylabel=r'Cumulative visible heat / photon energy',
          title='Conditional scalar debris, 20 MeV source')
ax[2].legend(frameon=False, fontsize=8, loc='upper left')
fig.savefig(ROOT/'FIGURE.png', dpi=170)
print(ROOT/'FIGURE.png')
