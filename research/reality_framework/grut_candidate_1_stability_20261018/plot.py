#!/usr/bin/env python3
"""Scientific figure from branch and nonlinear trajectory evidence."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
branches=json.loads((ROOT/'branch_results.json').read_text())
states=[x for x in branches['states'] if x['label']=='continuation']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
ax=axes[0,0]
low=[x for x in states if x['Omega']<=.96];high=[x for x in states if x['Omega']>=.95]
for seq,color,label in [(low,'#007c70','Lower-frequency branch'),(high,'#ca5c20','Higher-frequency branch')]:
    ax.plot([x['Qbar'] for x in seq],[x['E_over_Q'] for x in seq],'o-',color=color,label=label,ms=4)
ax.axhline(1,color='#777777',ls='--',lw=1,label='Free carriers')
ax.scatter([399.4384624],[.9110865197],s=65,marker='*',color='#1e4f8a',zorder=4,label='Candidate 1 preparation')
ax.set(xlabel=r'Charge $\bar Q$',ylabel=r'Energy per charge $\bar E/\bar Q$',title='One law has distinct stationary branches')
ax.legend(frameon=False,fontsize=8)
colors={'control':'#333333','small':'#1f78b4','stress':'#008675','time':'#9b5b97','coarse':'#999999','space':'#dc7c25','domain':'#736daa'}
runs={name:json.loads((ROOT/f'evolution_{name}.json').read_text()) for name in colors}
ax=axes[0,1]
for name in ['control','small','stress']:
    rows=runs[name]['trajectory']
    ax.plot([r['tau'] for r in rows],[r['core_Q_fraction']['8'] for r in rows],label=name,color=colors[name],lw=2)
ax.set(xlabel=r'Time $g_X v t$',ylabel='Charge retained within radius 8',title='Follow the core, not only conserved totals')
ax.legend(frameon=False)
ax=axes[1,0]
for name in ['control','small','stress','time']:
    rows=runs[name]['trajectory']
    ax.plot([r['tau'] for r in rows],[r['phase_aligned_profile_error'] for r in rows],label=name,color=colors[name],lw=1.6,ls='--' if name=='time' else '-')
ax.set(xlabel=r'Time $g_X v t$',ylabel='Phase-aligned profile deviation',title='A strained state oscillates during the test')
ax.legend(frameon=False)
ax=axes[1,1]
for name in ['stress','time','space','domain']:
    rows=runs[name]['trajectory'];E0=rows[0]['Ebar']
    ax.plot([r['tau'] for r in rows],[(r['Ebar']/E0-1)*1e5 for r in rows],label=name,color=colors[name],lw=1.5)
ax.set(xlabel=r'Time $g_X v t$',ylabel=r'Relative energy change $\times10^5$',title='Time, grid and box checks remain visible')
ax.legend(frameon=False,ncol=2)
for ax in axes.flat:ax.grid(alpha=.16)
fig.suptitle('GRUT Candidate 1 · finite-charge and nonlinear stability tests',fontsize=15,fontweight='bold')
fig.savefig(ROOT/'STABILITY_FIGURE.png',dpi=170)
fig.savefig(ROOT/'STABILITY_FIGURE.pdf')
print('Created STABILITY_FIGURE.png and STABILITY_FIGURE.pdf')
