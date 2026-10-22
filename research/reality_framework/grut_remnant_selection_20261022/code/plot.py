#!/usr/bin/env python3
"""Static scientific figure from the recorded calculation."""
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
R=json.loads((ROOT/'RESULTS.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(2,2,figsize=(12,8.5),constrained_layout=True)
colors=['#2563a6','#b64b37']
for T,col in zip([1e9,1e12],colors):
    rows=[x for x in R['histories'] if x['T_initial_GeV']==T and x['initial_du_dN']==0]
    xx=[x['initial_u']/math.log(10) for x in rows]
    yy=[x['final_u']/math.log(10) for x in rows]
    ax[0,0].plot(xx,yy,'o-',ms=4,lw=1,color=col,label=f'$T_i=10^{{{int(math.log10(T))}}}$ GeV')
ax[0,0].plot([-44,43],[-44,43],':',color='.55',label='Unchanged position')
ax[0,0].set(xlabel=r'Initial $\log_{10}(\sigma/\sigma_0)$',ylabel=r'At 160 GeV: $\log_{10}(\sigma/\sigma_0)$',
            title='Initial position survives across a wide range')
ax[0,0].legend(loc='upper left',fontsize=9)
ax[0,0].text(.39,.08,r'$u_f/u_i\simeq1-5.14\times10^{-60}$ near minimum',transform=ax[0,0].transAxes,
             fontsize=9,bbox=dict(facecolor='#f1f5f9',edgecolor='none',pad=5))

for u,col in zip([60,90,95],['#178579','#2563a6','#b64b37']):
    d=np.genfromtxt(ROOT/f'data/history_u{u}.csv',delimiter=',',names=True)
    ax[0,1].plot(d['N'],d['u'],color=col,label=f'$u_i={u}$')
    ax[1,0].semilogy(d['N'],d['Omega_radiation'],color=col,label=f'$u_i={u}$')
ax[0,1].axhline(0,color='.45',ls=':',label='Required minimum')
ax[0,1].set(xlabel='E-folds after initial radiation state',ylabel=r'$u=\ln(\sigma/\sigma_0)$',
            title=r'Large displacements remain far from $u=0$',ylim=(-5,101))
ax[0,1].legend(fontsize=9)
ax[1,0].set(xlabel='E-folds after initial radiation state',ylabel=r'Radiation fraction $\Omega_r$',
            title='Scalar domination is included in the evolution',ylim=(1e-26,1.3))
ax[1,0].legend(fontsize=9)

bb=np.linspace(.25,4,600)
ww=np.where(bb<1,1/3,4/(3*bb*bb)-1)
ax[1,1].axvspan(.25,1,color='#dceaf6',label='Radiation tracking: $b<1$')
ax[1,1].axvspan(math.sqrt(2),4,color='#f8e6df',label=r'Acceleration: $b>\sqrt{2}$')
ax[1,1].plot(bb,ww,color='#304154',lw=2,label='Asymptotic total equation of state')
ax[1,1].axhline(-1/3,color='.5',ls=':')
ax[1,1].plot([3],[4/27-1],'o',color='#b64b37')
ax[1,1].set(xlabel='Fixed pole coefficient $b$',ylabel=r'$w_{\rm total}$',ylim=(-1,.5),xlim=(.25,4),
            title='Tracking and acceleration require different slopes')
ax[1,1].legend(loc='upper right',fontsize=8)
fig.suptitle('Remnant position from reheating to the electroweak epoch',fontsize=16)
fig.savefig(ROOT/'overview.png',dpi=170)
plt.close(fig)
