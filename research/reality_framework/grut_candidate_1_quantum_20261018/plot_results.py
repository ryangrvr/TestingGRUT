#!/usr/bin/env python3
"""Static scientific figure from the recorded scalar-loop diagnostics."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def main():
    profile = np.genfromtxt(HERE/'potential_profile_flags.csv',delimiter=',',names=True)
    vals = np.genfromtxt(HERE/'flow_scale_values.csv',delimiter=',',names=True)
    flow = json.loads((HERE/'flow_results.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig, ax = plt.subplots(1,3,figsize=(12.8,3.6),layout='constrained')
    rho = profile['rho']
    wall = profile['negative_mode_count'] > 0
    ax[0].plot(rho, profile['F'],label=r'$F=f/v$',color='#2363a5')
    ax[0].plot(rho, profile['phi'],label=r'$S=\phi/v$',color='#ca6a17')
    ax[0].fill_between(rho,0,2.05,where=wall,color='#bd3942',alpha=.13)
    ax[0].set(xlim=(0,10),ylim=(0,2.05),xlabel=r'$\rho=g_Xvr$',ylabel='Field amplitude',title='Prepared classical solution')
    ax[0].legend(frameon=False)
    ax[1].plot(rho,profile['mass2_0'],color='#2363a5')
    ax[1].fill_between(rho,profile['mass2_0'],0,where=wall,color='#bd3942',alpha=.28)
    ax[1].axhline(0,color='#333333',linewidth=.8)
    ax[1].set(xlim=(0,10),xlabel=r'$\rho=g_Xvr$',ylabel=r'Lowest static Hessian eigenvalue / $v^2$',title='Constant-field wall limitation')
    ax[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    for i,label,color in [(0,'Vacuum background','#2363a5'),(3,'Homogeneous bulk background','#ca6a17')]:
        rows=vals[vals['background_index']==i]
        reference=rows[rows['mu_over_mu0']==1]['V0_plus_CW_over_v4'][0]
        ax[2].plot(rows['mu_over_mu0'],(rows['V0_plus_CW_over_v4']-reference)*1e9,'o-',label=label,color=color,markersize=4)
    ax[2].set_xscale('log',base=2)
    ax[2].set(xticks=[.25,.5,1,2,4],xticklabels=['1/4','1/2','1','2','4'],xlabel=r'$\mu/\mu_0$',ylabel=r'Change in $(V_0+V_1)$ / $(10^{-9}v^4)$',title='Residual scale dependence')
    ax[2].axhline(0,color='#333333',linewidth=.8)
    ax[2].legend(frameon=False,fontsize=8,loc='lower left')
    for a in ax:
        a.grid(alpha=.14)
    fig.savefig(HERE/'quantum_diagnostics.png',dpi=190)
    plt.close(fig)


if __name__=='__main__':
    main()
