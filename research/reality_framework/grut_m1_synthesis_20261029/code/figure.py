"""Draw the locally tested phase holder and conditional transport scales."""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': 'white', 'savefig.facecolor': 'white'})
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6), constrained_layout=True)
phase = json.loads((ROOT/'phase_holder.json').read_text())['states']
axes[0].plot([v['xi']/1e8 for v in phase], [v['delta_sigma_minus_6x'] for v in phase], 'o-', color='#29527a', lw=2)
axes[0].axhline(-np.pi/2, color='#c85d3f', ls='--', label='Requested phase')
axes[0].scatter([4], [next(v for v in phase if v['xi']==4e8)['delta_sigma_minus_6x']], s=90, facecolor='#e1a546', edgecolor='black', zorder=5)
axes[0].set(xlabel=r'Holder coefficient $\xi_6$  [$10^8$]', ylabel=r'Relative phase $\delta$ [rad]', title='Local de Sitter minimum')
axes[0].legend(frameon=False, loc='lower right')
ref=np.load(ROOT/'reference_long.npz')['original_prepared_no_holder'];hold=np.load(ROOT/'holder_long.npz')['local_equilibrium_holder_4e8']
axes[1].plot(ref[:,0],-ref[:,4]/1000,lw=1.5,color='#53798c',label='Original preparation')
axes[1].plot(hold[:,0],-hold[:,4]/1000,lw=1.5,color='#ce7e36',label=r'Holder, $\xi_6=4\times10^8$')
axes[1].set(xlabel=r'Reheating time $H_i t$',ylabel=r'Comoving $-N_X$  [$10^3$]',title='Reheating charge (spectator)')
axes[1].set_xlim(0,800);axes[1].legend(frameon=False,loc='lower right',fontsize=9)
f=np.load(ROOT/'free_stream.npz')
for n,c in [(64,'#9caab3'),(96,'#2d836f')]:
 x=f[f'lambda_Mpc_{n}']*1e6;w=f[f'weights_{n}'];idx=np.argsort(x);axes[2].plot(x[idx],np.cumsum(w[idx])/np.sum(w),color=c,label=f'{n} source bins')
axes[2].axvline(35.5063,color='#2d836f',ls=':',lw=1)
axes[2].set(xlabel='Collisionless distance to equality [comoving pc]',ylabel='Fraction of free $X$',title='Conditional free flight')
axes[2].set_xlim(0,120);axes[2].set_ylim(0,1.02);axes[2].legend(frameon=False,loc='lower right')
for ax in axes:ax.grid(alpha=.17)
fig.suptitle('M1 synthesis: phase selection, reheating, and transport',fontsize=15,fontweight='bold')
fig.text(.5,-.035,'Local scalar model and homogeneous gas assumptions; these are not a final isocurvature or matter-power prediction.',ha='center',fontsize=10,color='#505050')
fig.savefig(ROOT/'synthesis.png',dpi=180,bbox_inches='tight')
