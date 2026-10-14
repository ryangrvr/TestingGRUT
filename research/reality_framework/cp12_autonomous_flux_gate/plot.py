"""Standalone scientific figure from saved trajectories; no fitted dimension."""
from pathlib import Path
import json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent


def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,3,figsize=(13,4.2),layout='constrained')
    colors=['#1769aa','#d66716','#247a50']
    for name,color in zip(['N4_Q4_B256_T20.json','N5_Q5_B128_T10.json','N6_Q6_B128_T10.json'],colors):
        d=json.loads((P/name).read_text());rs=[r for r in d['rows'] if r['T']<=10]
        t=np.array([r['T'] for r in rs]);err=np.array([r['bounded_fraction_error_upper'] for r in rs])
        for ax,key in zip(axes[:2],['removable_fraction','flux_protected_fraction']):
            y=np.array([r[key] for r in rs])
            ax.plot(t,y,label=f"N = Q = {d['N']}",color=color,lw=1.6)
            ax.fill_between(t,np.maximum(0,y-err),np.minimum(1,y+err),color=color,alpha=.16)
    axes[0].axhline(.1,color='#777777',lw=1,ls=':',label='Low-r diagnostic threshold')
    axes[0].set(title='Partial protection, not low removability',ylabel=r'$\langle r\rangle/Q$',ylim=(0,1.03))
    axes[1].set(title='Flux protection develops autonomously',ylabel=r'$\langle k\rangle/Q$',ylim=(0,.3))
    for name,color in zip(['RECORDS_N4_Q4_B256.json','RECORDS_N5_Q5_B128.json'],colors):
        d=json.loads((P/name).read_text());rs=d['rows'];cap=math.log2(d['graph_count'])
        t=[r['T'] for r in rs]
        axes[2].plot(t,[r['graph_entropy_bits']/cap for r in rs],marker='o',ms=3,color=color,
                     label=f"Graph entanglement, N={d['N']}")
        axes[2].plot(t,[r['graph_classical_entropy_bits']/cap for r in rs],ls='--',color=color,alpha=.7,
                     label=f"Graph-basis entropy, N={d['N']}")
    axes[2].set(title='Internal records; broad graph distribution',ylabel='Fraction of maximum graph entropy',ylim=(0,1.05))
    for ax in axes:
        ax.set(xlabel='Time T',xlim=(0,10));ax.grid(alpha=.16);ax.legend(fontsize=8,loc='best',frameon=False)
    fig.suptitle('Unforced frozen CP-12 from a uniform neutral graph state | finite systems only',fontsize=12)
    fig.savefig(P/'flux_gate.png',dpi=180)
    fig.savefig(P/'flux_gate.svg')
    print('flux_gate.png and flux_gate.svg generated')


if __name__=='__main__':main()
