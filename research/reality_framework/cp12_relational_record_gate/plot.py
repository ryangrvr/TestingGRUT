"""Static scientific figure: protected loops are not extended space."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

P=Path(__file__).resolve().parent


def main():
    ns=np.array([4,5,6]);rows=[]
    for n,b in [(4,256),(5,128),(6,128)]:
        d=json.loads((P/f'RELATIONAL_N{n}_Q{n}_B{b}.json').read_text());rows.append(d['rows'][-1])
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','font.family':'DejaVu Sans'})
    fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    colors=['#126b8e','#bf7c24','#617b44']
    for ax,metric,title in [(axs[0,0],'has_cycle','a   Protected loops form'),
                            (axs[0,1],'connected','b   Spanning protected cores become rarer')]:
        actual=np.array([r['means'][metric] for r in rows]);err=np.array([r['bounded_operator_error_upper'] for r in rows])
        ax.errorbar(ns,actual,yerr=err,color=colors[0],marker='o',lw=2,capsize=4,label='Evolved state')
        for label,text,color,style in [('parent_thinning','Random thinning of same parent',colors[1],'--'),
                                       ('gauss_parent_resources','Gauss + parent + resource null',colors[2],':')]:
            ax.plot(ns,[r['controls'][label]['means'][metric] for r in rows],style,
                    color=color,marker='s',ms=4,label=text)
        ax.set(xticks=ns,xlabel='Vertices N   (Q = N, T = 10)',ylabel='Probability',ylim=(0,.36))
        ax.set_title(title,loc='left',fontweight='bold',pad=12)
        ax.grid(axis='y',alpha=.2)
    axs[0,0].legend(fontsize=8,loc='upper left')
    ax=axs[1,0];width=.22
    for c,(text,color) in enumerate(zip(['Core distance 2','Core distance 3+','Different core components'],
                                        [colors[0],colors[1],'#b84d5b'])):
        y=100*np.array([r['shortcut']['contrasts_per_Q'][c] for r in rows])
        err=100*np.array([r['bounded_operator_error_upper'] for r in rows])
        ax.bar(ns+(c-1)*width,y,width=width,color=color,label=text)
        ax.errorbar(ns+(c-1)*width,y,yerr=err,fmt='none',ecolor='#29333a',capsize=3,lw=.8)
    ax.axhline(0,color='#333333',lw=.8)
    ax.set(xticks=ns,xlabel='Vertices N   (Q = N, T = 10)',
           ylabel='Neutral edges / Q: excess over\nuniform completion (percentage points)')
    ax.set_title('c   Cross-component connections are enhanced at N = 5, 6',loc='left',fontweight='bold',pad=12,fontsize=10)
    ax.legend(fontsize=8,loc='lower left');ax.grid(axis='y',alpha=.2)
    ax=axs[1,1];r=rows[-1];x=np.arange(6)
    d=json.loads((P/'RELATIONAL_N6_Q6_B128.json').read_text())
    def balls(a):return [6*a[f'ball_{j}_fraction'] for j in x]
    ax.errorbar(x,balls(r['means']),yerr=6*r['bounded_operator_error_upper'],
                color=colors[0],marker='o',lw=2,capsize=3,label='Protected core')
    ax.plot(x,balls(r['full_graph_means']),color='#b84d5b',marker='s',lw=2,label='Full evolved graph')
    ax.plot(x,balls(d['full_graph_uniform_fixed_Q_means']),'--',color='#77838b',marker='.',label='Uniform fixed-Q full graph')
    ax.set(xticks=x,xlabel='Intrinsic radius r',ylabel='Mean ball volume',ylim=(.8,6.3))
    ax.set_title('d   Tiny graphs provide no dimension-scaling window',loc='left',fontweight='bold',pad=12,fontsize=10)
    ax.legend(fontsize=8,loc='lower right');ax.grid(axis='y',alpha=.2)
    fig.suptitle('CP-12 relational record gate: correlations without established spatial locality',fontsize=14,fontweight='bold')
    fig.savefig(P/'relational_gate.png',dpi=180)
    fig.savefig(P/'relational_gate.svg')
    plt.close(fig)


if __name__=='__main__':main()
