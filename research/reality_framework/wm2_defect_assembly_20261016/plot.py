"""Standalone scientific figure from exact construction/control data."""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--png');args=ap.parse_args()
    r=json.loads((P/'RESULTS.json').read_text());e=json.loads((P/'EXACT_CONFIGURATION.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
        'axes.spines.top':False,'axes.spines.right':False,
        'svg.hashsalt':'WM2-S1-20261016'})
    fig,axes=plt.subplots(1,3,figsize=(13.3,4.8),layout='constrained')
    dark='#263543';blue='#226ca3';purple='#8843a3';orange='#c17524';gray='#7d8790'
    barrier=e['activation_barriers'][-1]['minimum_maximum_defect_count_by_initial_distance']
    d=np.array([int(x) for x in barrier]);height=np.array(list(barrier.values()))
    ax=axes[0];ax.plot(d,d,'--',color=gray,label='Proved lower bound')
    ax.plot(d,height,'o',color=blue,label='Complete minimax search, N=7')
    ax.set(xlabel='Generated protected distance',ylabel='Minimum peak free-event defects',
        title='A distant bond needs a growing barrier',xticks=d,yticks=d)
    ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=8.2,loc='upper left')
    blocks=r['Schur'][-1]['all_geometry_receiver_blocks']
    ds=np.array([x['min_distance'] for x in blocks]);actual=np.array([x['operator_norm'] for x in blocks])
    bound=1/(13*6)*(5/31)**ds
    ax=axes[1];ax.semilogy(ds,bound,'--',color=purple,label='Uniform bound, all low energies')
    ax.semilogy(ds,actual,'s-',color=blue,label='Complete Schur block at E=0, N=6')
    ax.set(xlabel='Minimum source–receiver distance',ylabel='Receiver/source operator norm',
        title='All virtual words and receiver records',xticks=ds)
    ax.grid(axis='y',alpha=.15);ax.legend(frameon=False,fontsize=8,loc='lower left')
    rows=r['formation'];ns=np.array([x['N'] for x in rows]);full=np.array([x['full_assembly_probability_upper_all_times'] for x in rows])
    half=np.array([x['extensive_assembly_probability_upper_all_times'] for x in rows])
    ax=axes[2];ax.loglog(ns,half,'o-',color=orange,label='At least half assembled')
    ax.loglog(ns,full,'s-',color=purple,label='Fully assembled')
    ax.set(xlabel='Number of events N',ylabel='Formation probability: all-time upper bound',
        title='Closed seeded formation fails',xticks=[2,4,8,16,32,64,128])
    ax.set_xticklabels(['2','4','8','16','32','64','128']);ax.minorticks_off()
    ax.grid(axis='y',alpha=.15);ax.legend(frameon=False,fontsize=8.5,loc='lower left')
    fig.suptitle('WM2-S1: generated adjacency and protected transport; no appreciable formed phase\n'
        'Single-end grammar and dimension are supplied inputs',fontsize=13,color=dark)
    fig.savefig(P/'WM2_CONSTRUCTION.svg',metadata={'Date':None,'Creator':'GRUT Builder WM2-S1 exact and finite controls'})
    target=P/'WM2_CONSTRUCTION.svg';target.write_text('\n'.join(x.rstrip() for x in target.read_text().splitlines())+'\n')
    if args.png:fig.savefig(args.png,dpi=160)
    plt.close(fig)


if __name__=='__main__':main()
