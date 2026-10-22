from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
p=Path(__file__).resolve().parents[1];i=json.loads((p/'INPUTS.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(1,2,figsize=(11.5,4.4),layout='constrained')
ns=np.linspace(.95,.99,300);V0=i['lam']*i['v_GeV']**4/4
for rr,col in [(0.,'#176a94'),(.036,'#be5033')]:
 rho=128*V0**2/(9*np.pi**2*i['As']*i['M_GeV']**4*(8*(1-ns)-rr)**2)
 ax[0].plot(ns,rho/i['rho_reference_GeV4'],label='r=%g'%rr,color=col,lw=2)
ax[0].axvspan(i['ns_snapshot']-2*i['ns_snapshot_sigma'],i['ns_snapshot']+2*i['ns_snapshot_sigma'],color='gray',alpha=.15,label='Dated tilt sensitivity range')
ax[0].axhline(1,color='#666666',ls='--',lw=1)
ax[0].set(xlabel='Spectral tilt ns',ylabel='Initial remnant / known reference density',ylim=(0,2),title='Conditional same-action slow-roll relation');ax[0].legend(loc='upper left',frameon=False)
d=np.genfromtxt(p/'fp_dx0.05_L80.csv',delimiter=',',names=True);x=d['canonical_over_Hdiv2pi']
for t,col in [(10,'#176a94'),(100,'#be5033')]:ax[1].plot(x,d['prob_N%g'%t]/.05,color=col,lw=2,label='N=%g'%t)
ax[1].axvline(0,color='#666666',lw=.7);ax[1].set(xlim=(-35,20),xlabel='Canonical field / (H / 2pi)',ylabel='Probability density',title='Fixed-inflaton stochastic control, c=10');ax[1].legend(frameon=False)
for a in ax:a.grid(alpha=.15)
fig.savefig(p/'overview.png',dpi=160);plt.close(fig)
