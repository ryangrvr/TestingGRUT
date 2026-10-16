"""Scientific plot of actual supplied control data; no conceptual GRUT data."""
from pathlib import Path
import argparse
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

P=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--png'); args=ap.parse_args()
    r=json.loads((P/'RESULTS.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
        'axes.spines.top':False,'axes.spines.right':False,
        'svg.hashsalt':'GRUT-protected-transport-20261016'})
    fig,axes=plt.subplots(1,2,figsize=(11.2,4.9),layout='constrained')
    blue='#2367a0'; purple='#8547a4'; dark='#273442'; gray='#77818b'
    fixed=[x for x in r['modules'] if x['kind']=='cycle' and x['N']<=64]
    latch=[x for x in r['latch'] if x['N']<=64]
    ax=axes[0]
    ax.semilogy([x['diameter'] for x in fixed],
        [x['far_E0_Feshbach_coupling'] for x in fixed],
        'o-',color=blue,label='Static cycle: complete virtual sum')
    ax.semilogy([x['distance'] for x in latch],
        [x['Schur'][0]['far'] for x in latch],
        's-',color=purple,label='GC1: virtual chord retained')
    floor=r['latch'][-1]['constants']['shortcut_floor_all_low_energies']
    ax.axhline(floor,color=purple,linestyle='--',linewidth=1,
        label='GC1 analytic floor, all low energies')
    ax.set(xlabel='Initial cycle distance',ylabel=r'$|H_F(0)_{j0}|$',
        title='Same reference metric; different virtual geometry',
        ylim=(1e-34,1e-3),xticks=[4,8,16,32])
    ax.grid(axis='y',alpha=.15); ax.legend(loc='lower left',frameon=False,fontsize=8.5)
    ax=axes[1]
    x=[z['distance'] for z in r['latch']]
    ax.loglog(x,[z['Schur'][0]['near'] for z in r['latch']],
        'o-',color=dark,label='GC1 nearby coupling')
    ax.loglog(x,[z['Schur'][0]['far'] for z in r['latch']],
        's-',color=purple,label='GC1 antipodal coupling')
    ax.axhline(floor,color=purple,linestyle='--',linewidth=1)
    ax.set(xlabel='Initial cycle distance (N/2)',ylabel=r'$|H_F(0)_{j0}|$',
        title='A uniform gap does not remove the distant floor',
        ylim=(3e-7,3e-3),xticks=x)
    ax.set_xticklabels([str(z) for z in x]); ax.minorticks_off()
    ax.text(.04,.53,'Gap ≥ 4.5 at every N\nOriginal geometry ≥ 97.59% at every time',
        transform=ax.transAxes,fontsize=9,color=dark,
        bbox=dict(boxstyle='round,pad=.5',facecolor='#f0f3f6',edgecolor='none'))
    ax.grid(axis='y',alpha=.15); ax.legend(loc='lower right',frameon=False,fontsize=8.5)
    fig.suptitle('Protected transport must include virtual geometry\n'
        'Supplied control models — no emergent spatial phase',
        fontsize=13,color=dark)
    fig.savefig(P/'PROTECTED_TRANSPORT.svg',metadata={'Date':None,
        'Creator':'GRUT Builder: exact supplied-control data'})
    target=P/'PROTECTED_TRANSPORT.svg'
    target.write_text('\n'.join(line.rstrip() for line in target.read_text().splitlines())+'\n')
    if args.png: fig.savefig(args.png,dpi=160)
    plt.close(fig)


if __name__=='__main__': main()
