"""Scientific overview; conditional diagnostics are explicitly labelled."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent


def load(name):return json.loads((P/name).read_text())


def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','font.family':'DejaVu Sans'})
    fig,axs=plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
    blue='#146d8a';gold='#bb8129';red='#b94c60';green='#587b41'
    ax=axs[0,0];d=load('EARLY_TIME.json');t=np.array([r['T'] for r in d['rows']])
    for key,text,color,coef,power in [('protected_triangle','Protected triangle',blue,8/5,6),
        ('occupied_flux_triangle','Entirely occupied flux triangle',red,8/1125,12)]:
        y=[r['values'][key] for r in d['rows']]
        ax.loglog(t,y,'o-',color=color,label=text,lw=2)
        ax.loglog(t,coef*t**power,':',color=color,alpha=.7,label=f'Exact leading $T^{{{power}}}$ term')
    ax.set(xlabel='Time T',ylabel='Probability',xticks=t)
    ax.set_xticklabels(['.02','.04','.08','.16','.32'])
    ax.set_title('a   Onset orders are realized (N = Q = 4)',loc='left',fontweight='bold',pad=12)
    ax.legend(fontsize=8);ax.grid(alpha=.2,which='major')
    ax=axs[0,1];ns=np.array([4,5,6]);rows=[]
    for n,b in [(4,256),(5,128),(6,128)]:rows.append(load(f'CHANNELS_N{n}_Q{n}_B{b}.json')['rows'][-1])
    for key,text,color in [('protected_cycle_any','Any protected loop',blue),
        ('occupied_flux_cycle_any','Entirely occupied flux support',gold),
        ('electric_circulation_cycle_any','Electric circulation, including holes',red)]:
        y=np.array([r['values'][key] for r in rows])
        iv=np.array([r['positive_score_intervals'][key] for r in rows])
        ax.errorbar(ns,y,yerr=np.vstack([y-iv[:,0],iv[:,1]-y]),fmt='o-',lw=2,capsize=3,color=color,label=text)
    ax.set(xticks=ns,xlabel='Vertices N   (Q = N, T = 10)',ylabel='Unconditioned probability',ylim=(0,.35))
    ax.set_title('b   A hole does not imply absence of circulation',loc='left',fontweight='bold',pad=12,fontsize=10)
    ax.legend(fontsize=8,loc='upper left');ax.grid(axis='y',alpha=.2)
    ax=axs[1,0];d=load('CIRCULATION_SURVIVAL.json')
    c=next(c for c in d['cases'] if c['preparation']=='neutral_sea_wait_10')
    ch=next(c for c in c['channels'] if c['initial_channel']=='occupied_flux')
    lag=np.array([r['lag'] for r in ch['rows']])
    for key,text,color in [('protected','Same triangle remains protected',blue),
        ('electric_circulation','Same triangle carries circulation',red),
        ('occupied_flux','Same triangle remains entirely occupied',gold)]:
        y=np.array([r['destinations'][key]['conditional_probability'] for r in ch['rows']])
        ax.plot(lag,y,'o-',color=color,label=text,lw=2)
    ax.set(xlabel='Lag',ylabel='Conditional probability',ylim=(0,1.04))
    ax.set_title('c   Diagnostic flux filter at T = 10 (N = Q = 4)',loc='left',fontweight='bold',pad=12,fontsize=10)
    ax.legend(fontsize=8,loc='lower left');ax.grid(axis='y',alpha=.2)
    ax=axs[1,1];d=load('CONDITIONAL_DENSITY.json');x=np.arange(3)
    y=np.array([r['density_susceptibility'] for r in d['rows']])
    iv=np.array([r['density_susceptibility_interval'] for r in d['rows']])
    ax.bar(x,y,color=[blue,green,gold],width=.58)
    ax.errorbar(x,y,yerr=np.vstack([y-iv[:,0],iv[:,1]-y]),fmt='none',ecolor='#26343b',capsize=5)
    ax.axhline(0,color='#333333',lw=.8)
    ax.set(xticks=x,xticklabels=['Main state\nwait 10','Spectral\nreference','Maximum-entropy\nreference'],
           ylabel='Density susceptibility',ylim=(-.020,.010))
    ax.set_title('d   Density response in that same filtered state',loc='left',fontweight='bold',pad=12,fontsize=10)
    ax.text(.02,.05,'Lag 0.5; same-triangle protection >99.3%\nPhase source at 0; mean density at 1, 2, 3',
            transform=ax.transAxes,fontsize=8,color='#39464d')
    ax.grid(axis='y',alpha=.2)
    fig.suptitle('CP-12 flux loops: formation, persistence and conditional density response',fontsize=14,fontweight='bold')
    fig.savefig(P/'flux_loop_gate.png',dpi=180);fig.savefig(P/'flux_loop_gate.svg')
    plt.close(fig)


if __name__=='__main__':main()
