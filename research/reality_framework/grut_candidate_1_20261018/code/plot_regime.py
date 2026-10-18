#!/usr/bin/env python3
"""Reproduce the scientific profile figure from saved calculation outputs."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
d=np.genfromtxt(ROOT/'results_profile.csv',delimiter=',',names=True)
reporter=json.loads((ROOT/'reporter_results.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False,'axes.titlesize':12})
fig,axes=plt.subplots(2,2,figsize=(11,7.8),layout='constrained')
r=d['rho']; mask=r<=14
ax=axes[0,0]
ax.plot(r[mask],d['F'][mask],label=r'Charge field $F=f/v$',color='#225ea8',lw=2.2)
ax.plot(r[mask],d['S'][mask],label=r'Organization $S=\phi/v$',color='#008b78',lw=2.2)
ax.set(xlabel=r'Radius $g_X v r$',ylabel='Field amplitude',title='A populated core depletes the organization field')
ax.legend(frameon=False,fontsize=9)
ax=axes[0,1]
ax.plot(r[mask],d['energy_over_g2v4'][mask],label='Energy density',color='#225ea8',lw=2)
ax.plot(r[mask],d['pr_over_g2v4'][mask],label='Radial pressure',color='#d95f0e',lw=1.8)
ax.plot(r[mask],d['pt_over_g2v4'][mask],label='Tangential pressure',color='#8266aa',lw=1.8)
ax.axhline(0,color='#bbbbbb',lw=.6)
ax.set(xlabel=r'Radius $g_X v r$',ylabel=r'Stress in units $g_X^2 v^4$',title='The gravitational source retains the pressures')
ax.legend(frameon=False,fontsize=9)
ax=axes[1,0]
ax.plot(r[mask],-d['Phi'][mask]*1e10,label=r'$-\Phi$ (clock potential)',color='#d95f0e',lw=2)
ax.plot(r[mask],-d['Psi'][mask]*1e10,label=r'$-\Psi$ (spatial potential)',color='#225ea8',lw=2)
ax.set(xlabel=r'Radius $g_X v r$',ylabel=r'Potential magnitude $\times 10^{10}$',title='Full-stress Einstein response in the weak field')
ax.legend(frameon=False,fontsize=9)
ax=axes[1,1]
k=np.array([x['k_over_gv'] for x in reporter['phase_results']])
delta=np.array([x['phase_delta_s_wave'] for x in reporter['phase_results']])
born=np.array([x['born_delta'] for x in reporter['phase_results']])
ax.plot(k,2*delta*1000,'o-',color='#008b78',label='Calculated elastic S phase',lw=2)
ax.plot(k,2*born*1000,'x',color='#333333',label='First Born approximation',ms=7)
ax.set(xlabel=r'Reporter momentum $p/(g_X v)$',ylabel='Relative phase (milliradians)',title='One profile predicts three reporter phases')
ax.legend(frameon=False,fontsize=9)
for ax in axes.flat:
    ax.grid(alpha=.15)
fig.suptitle('GRUT Candidate 1 · calculated leading regime',fontsize=16,fontweight='bold')
fig.savefig(ROOT.parent/'REGIME_FIGURE.png',dpi=170)
fig.savefig(ROOT.parent/'REGIME_FIGURE.pdf')
print('Created REGIME_FIGURE.png and REGIME_FIGURE.pdf')
