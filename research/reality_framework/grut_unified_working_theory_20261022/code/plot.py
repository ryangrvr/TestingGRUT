#!/usr/bin/env python3
"""Scientific overview of the single-action result, not an empirical fit."""
import json
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binom

ROOT=Path(__file__).resolve().parents[1]
R=json.loads((ROOT/'RESULTS.json').read_text());I=json.loads((ROOT/'INPUTS.json').read_text())
rows={x['Q']:x for x in R['profiles']}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(2,2,figsize=(11.5,8.5),constrained_layout=True)
colors=['#147e76','#2364a6','#b44737']
for Q,col in zip(I['main_charges'],colors):
    p=np.genfromtxt(ROOT/f'data/profile_Q{Q}.csv',delimiter=',',names=True)
    ax[0,0].plot(p['rho'],p['S'],color=col,label=f'Q={Q:,}')
    row=rows[Q];mass=row['mass_over_gv'];exact=-row['scattering_length_over_inverse_gv']
    lower=-row['exact_scattering_length_upper_bound'];upper=-row['exact_scattering_length_lower_bound']
    ax[0,1].errorbar(mass,exact,yerr=[[exact-lower],[upper-exact]],fmt='o',capsize=5,color=col)
ax[0,0].set(xlim=(0,12),ylim=(0,1.05),xlabel=r'Radius $g_Xv_\phi r$',ylabel=r'Organization field $\phi/v_\phi$',
            title='One field profile supports each object')
ax[0,0].legend(fontsize=9)
mass=np.linspace(0,72000,200)
coefficient=I['h_D']**2*I['g_X']**2/(4*math.pi*I['lambda_phi'])
ax[0,1].plot(mass,coefficient*mass,'--',color='.35',label='Exact classical Born sum rule')
ax[0,1].set(xlabel=r'Object mass $E/(g_Xv_\phi)$',ylabel=r'$-g_Xv_\phi a_D$',title='Mass constrains the probe response')
ax[0,1].text(.03,.92,'Points: exact scalar scattering\nBars: controlled intervals',va='top',transform=ax[0,1].transAxes,fontsize=9)
ax[0,1].legend(loc='lower right',fontsize=8)
for Q,col in zip(I['readout_charges'],colors):
    p=np.genfromtxt(ROOT/f'data/phases_Q{Q}.csv',delimiter=',',names=True)
    ax[1,0].plot(p['k_over_gv'],p['delta_l0']-p['delta_l1'],color=col,label=f'Q={Q:,}')
ax[1,0].set(xlabel=r'Probe momentum $k/(g_Xv_\phi)$',ylabel=r'Relative phase $\delta_0-\delta_1$ (rad)',
            title='The existing probe reads a charge label')
ax[1,0].legend(fontsize=9)
info=R['readout'];nn=np.unique(np.geomspace(1,25000,350).astype(int))
optimal=.5*(1+np.sqrt(-np.expm1(nn*math.log(info['overlap_squared']))))
p0,p1=sorted(info['measurement_plus_probabilities'],reverse=True)
cut=np.ceil(nn*np.log((1-p1)/(1-p0))/(np.log(p0/p1)+np.log((1-p1)/(1-p0)))).astype(int)
fixed=1-.5*(binom.cdf(cut-1,nn,p0)+binom.sf(cut-1,nn,p1))
ax[1,1].semilogx(nn,optimal,color='#2364a6',label='Optimal quantum measurement')
ax[1,1].semilogx(nn,fixed,color='#b44737',label='Specified fixed measurement')
ax[1,1].axhline(.95,color='.5',ls=':')
ax[1,1].scatter([info['ideal_independent_probe_count_optimal_95'],info['fixed_measurement_probe_count_95']],[.95,.95],
                 color=['#2364a6','#b44737'],s=25)
ax[1,1].set(xlabel='Number of ideal independent probes',ylabel='Label discrimination success',ylim=(.49,1.01),
            title='Readout assumes an unchanged background')
ax[1,1].legend(loc='upper left',fontsize=8)
fig.suptitle('Candidate 1: one profile links matter, mass and a readable signal',fontsize=15)
fig.savefig(ROOT/'overview.png',dpi=170)
plt.close(fig)
