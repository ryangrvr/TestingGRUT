#!/usr/bin/env python3
"""Scientific plots of recorded mode controls and scalar pole thresholds."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def main():
    modes=json.loads((HERE/'mode_results.json').read_text())
    poles=json.loads((HERE/'neutral_pole_results.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,3,figsize=(12.9,3.6),layout='constrained')
    hs=[];nus=[]
    for N in [120,180,400,800]:
        s=next(s for s in modes['spectra'] if s['N']==N and s['L']==30 and s['ell']==1 and s['method'] in ['dense','sparse'])
        ns=[r['s_imag'] for r in s['smallest_roots'] if r['symmetry']=='translation' and r['s_imag']>0]
        hs.append(30/(N+1));nus.append(min(ns))
    ax[0].plot(hs,nus,'o-',color='#2363a5',label='Computed translation mode')
    ax[0].plot([0,max(hs)],[0,np.mean(np.array(nus[-2:])/np.array(hs[-2:]))*max(hs)],'--',color='#777777',label=r'Proportional to $h$')
    ax[0].set(xlabel=r'Radial spacing $h$',ylabel=r'$\nu/(g_Xv)$',title='Translation zero-mode recovery',xlim=(0,.265))
    ax[0].legend(frameon=False,fontsize=8,loc='upper left')
    for L,color in [(30,'#ca6a17'),(40,'#2363a5')]:
        vals=[]
        for ell in range(6):
            s=next(s for s in modes['spectra'] if s['ell']==ell and s['L']==L and s['method'] in ['matched_box','matched_ell5'])
            vals.append(min(r['s_imag'] for r in s['smallest_roots'] if r['s_imag']>1e-5 and r['symmetry'] is None))
        ax[1].plot(range(6),vals,'o-',label=f'L={L}, h=0.05',color=color)
    ax[1].axhline(1-modes['protocol']['Omega'],color='#333333',linestyle='--',label='Charged continuum threshold')
    ax[1].set(xlabel=r'Angular sector $\ell$',ylabel=r'$\nu/(g_Xv)$',title='Box levels move toward continuum',xticks=list(range(6)),ylim=(.197,.253))
    ax[1].legend(frameon=False,fontsize=7.5,loc='upper left')
    rows=[poles['results'][i] for i in [0,1,3]]
    x=np.arange(3);width=.34
    cur=[(r['strict_curvature_squared_over_v2']/r['tree_mass_squared_over_v2']-1)*1e6 for r in rows]
    pol=[(r['strict_pole_squared_over_v2']/r['tree_mass_squared_over_v2']-1)*1e6 for r in rows]
    ax[2].bar(x-width/2,cur,width,label='Zero-momentum curvature',color='#9ab6d4')
    ax[2].bar(x+width/2,pol,width,label='Propagator pole',color='#2363a5')
    ax[2].axhline(0,color='#333333',linewidth=.7)
    ax[2].set(xticks=x,xticklabels=[r'$\phi$',r'$X$',r'$D$'],ylabel=r'Relative squared-mass shift (ppm)',title='Scalar-loop vacuum thresholds')
    ax[2].legend(frameon=False,fontsize=7.5,loc='upper right')
    for a in ax:a.grid(alpha=.12,axis='y');a.set_axisbelow(True)
    fig.savefig(HERE/'fluctuation_diagnostics.png',dpi=190)
    plt.close(fig)


if __name__=='__main__':main()
