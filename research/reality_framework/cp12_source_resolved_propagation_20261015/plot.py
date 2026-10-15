"""Static scientific figure from recorded full-H receiver intervals."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent


def read(name):return json.loads((P/name).read_text())


def pair_curve(ax,data,receiver,color,label):
    rows=data['rows'];t=np.array([r['lag'] for r in rows])
    values=[r['fixed_receiver_pairs'][str(receiver)]['density_response'] for r in rows]
    mid=np.array([abs(v['value']) for v in values]);I=np.array([v['absolute_interval'] for v in values])
    resolved=I[:,0]>0
    ax.plot(t[resolved],mid[resolved],'-o',color=color,label=label,ms=4,lw=1.7)
    if resolved.any():ax.errorbar(t[resolved],mid[resolved],
        yerr=np.stack([mid[resolved]-I[resolved,0],I[resolved,1]-mid[resolved]]),
        fmt='none',color=color,capsize=3,lw=1)
    if (~resolved).any():
        upper=I[~resolved,1]
        ax.errorbar(t[~resolved],upper,yerr=upper*.35,uplims=True,fmt='v',color=color,ms=4)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--png',type=Path);args=ap.parse_args()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
        'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'grut-source-20261015'})
    colors=['#176d95','#c34a2c','#6a5898']
    fig,axes=plt.subplots(1,3,figsize=(13.1,4.3),layout='constrained')
    n=np.array([4,5,6]);neutral=[read(f'NEUTRAL_N{i}_Q{i+1}_B{192 if i==6 else 256}.json')['rows'][0] for i in n]
    refined=[read(f'REFINED_RECORD_N{i}_Q{i+1}_B{192 if i==6 else 256}_W5.json') for i in n]
    ring=np.array([d['formation_probability'] for d in neutral])
    chord=np.array([next(r['class_probability'] for r in d['classes'] if r['chord_core_distance']==i//2) for i,d in zip(n,refined)])
    endpoint=chord*2/n
    ax=axes[0]
    for y,c,label in [(ring,colors[0],'Any spanning ring'),(chord,colors[2],'Longest-chord ring class'),(endpoint,colors[1],'Fixed source on that chord')]:
        ax.plot(n,100*y,'-o',color=c,label=label,ms=4)
    ax.axhline(1,color='#888888',lw=.8,ls=':');ax.text(4.05,1.12,'1% finite-test benchmark',fontsize=8,color='#666666')
    ax.set(yscale='log',xticks=n,xlabel='Number of events N',ylabel='Original-neutral probability (%)',
        title='Formation at waiting time 5')
    ax.legend(frameon=False,fontsize=8,loc='upper right')
    d5=read('RESPONSE_N5_Q6_B256_W5_G200000_F2_L10.json')
    d6=read('RESPONSE_N6_Q7_B192_W5_G200000_F3_L8.json')
    pair_curve(axes[1],d5,2,colors[1],'Receiver 2: ring d=2, host d=1')
    pair_curve(axes[1],d5,3,colors[0],'Receiver 3: ring d=2, host d=2')
    pair_curve(axes[2],d6,3,colors[1],'Receiver 3: ring d=3, host d=1')
    pair_curve(axes[2],d6,2,colors[0],'Receiver 2: ring d=2, host d=2')
    pair_curve(axes[2],d6,1,colors[2],'Receiver 1: ring d=1, host d=1')
    for ax,title in zip(axes[1:],['Five sites: equal ring distance','Six sites: farther receiver bypass']):
        ax.set(yscale='log',xlabel='Lag after source pulse',ylabel=r'$|\chi_{0j}|$',title=title)
        ax.set_xticks([.0625,.125,.25,.5],['1/16','1/8','1/4','1/2'])
        ax.legend(frameon=False,fontsize=8,loc='lower right')
    for ax in axes:ax.grid(axis='y',alpha=.17)
    fig.suptitle('Full frozen CP-12 dynamics • auxiliary Q=N+1 sequence • recorded core and initial host',fontsize=12)
    fig.savefig(P/'SOURCE_RESPONSE.svg',metadata={'Date':None})
    svg=P/'SOURCE_RESPONSE.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    if args.png:fig.savefig(args.png,dpi=170)


if __name__=='__main__':main()
