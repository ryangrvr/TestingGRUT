#!/usr/bin/env python3
"""Figures expose both construction and unresolved convergence."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]


def main():
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.facecolor':'white','axes.titlesize':12})
    core=json.loads((ROOT/'CORE_RESULTS.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(11.7,4.5))
    for name,color in [('A','#bb683d'),('A_fine','#4077a4'),('A_finer','#367760')]:
        row=next(x for x in core['runs'] if x['name']==name)
        data=np.loadtxt(ROOT/row['history_file'],delimiter=',',skiprows=1)
        label=f"h={row['parameters']['h']:g}, dt={row['parameters']['dt']:g}"
        axes[0].plot(data[:,0],data[:,12]/data[0,2],color=color,lw=1.6,label=label)
        err=np.maximum(abs((data[:,1]+data[:,3])/data[0,1]-1),1e-16)
        axes[1].semilogy(data[:,0],err,color=color,lw=1.2,label=label)
    axes[0].set(xlabel=r'Time $g_Xvt$',ylabel=r'Charge inside $\rho=60$ / initial charge',
                title='Retained charge: unresolved spatial convergence',ylim=(.2,.65))
    axes[1].set(xlabel=r'Time $g_Xvt$',ylabel='Relative energy + receiver error',
                title='Hamiltonian accounting improves with refinement',ylim=(1e-7,1e-3))
    for ax in axes:ax.grid(alpha=.18);ax.legend(fontsize=8)
    fig.suptitle('Same scalar action and prepared Gaussian patch; these are finite spherical prototypes',fontsize=13)
    fig.tight_layout(rect=[0,0,1,.94]);fig.savefig(ROOT/'core_convergence.png',dpi=160)
    d=json.loads((ROOT/'DAUGHTER_RESULTS.json').read_text())
    rows=d['free_streaming']['rows'];fig,ax=plt.subplots(figsize=(8.5,4.8))
    for gamma,color in [(1.,'#367760'),(10.,'#4077a4'),(1000.,'#bb683d')]:
        r=[r for r in rows if r['parent_gamma']==gamma]
        ax.loglog([v['a_formation'] for v in r],[v['mean_free_stream_Mpc_today'] for v in r],
                  '-o',color=color,label=fr'Parent $\gamma_s={gamma:g}$')
    cold=[r for r in rows if r['parent_gamma']==1.]
    ax.annotate('20 MeV: 51 pc',(cold[0]['a_formation'],cold[0]['mean_free_stream_Mpc_today']),
                xytext=(10,12),textcoords='offset points',fontsize=10)
    ax.annotate('1 keV: 0.64 Mpc',(cold[2]['a_formation'],cold[2]['mean_free_stream_Mpc_today']),
                xytext=(12,-23),textcoords='offset points',fontsize=10)
    ax.set(xlabel='Supplied formation scale factor',ylabel='Mean comoving travel to today [Mpc]',
           title='Conditional daughter propagation; no Lyman-alpha fit')
    ax.grid(alpha=.18,which='both');ax.legend(fontsize=10)
    fig.text(.13,.02,'At fixed source energy, boosted parents produce fewer daughters; their matter fraction is recomputed.',fontsize=8.5)
    fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(ROOT/'daughter_history.png',dpi=160)


if __name__=='__main__':main()
