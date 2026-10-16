"""Scientific figure from actual full-host data; no conceptual front data."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent


def read(name):return json.loads((P/name).read_text())


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--png',type=Path);a=ap.parse_args()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
        'axes.spines.right':False,'svg.hashsalt':'grut-active-host-20261016'})
    fig,axs=plt.subplots(1,3,figsize=(13,4.2),layout='constrained');colors=['#176d95','#b14d2d']
    for n,b,c in [(5,256,colors[0]),(6,192,colors[1])]:
        d=read(f'HOST_N{n}_Q{n+1}_B{b}.json');rows=d['rows']
        axs[0].plot([r['time'] for r in rows],[100*r['total_variation_from_uniform'] for r in rows],'-o',color=c,label=f'N={n}',ms=4)
    axs[0].set(xlabel='Time from original neutral ray',ylabel='Total variation from uniform (%)',title='Genuine but limited host bias')
    axs[0].legend(frameon=False)
    d=read('HOST_N6_Q7_B192.json');x=np.arange(4)
    for shift,row,c,label in [(-.18,d['rows'][0],'#91b7c9','Original neutral'),(.18,d['rows'][-1],colors[1],'Full H, time 5')]:
        geo=row['geometry'];p=[geo['disconnected_probability']]+[geo['diameter_distribution'].get(str(k),0.) for k in [2,3,4]]
        axs[1].bar(x+shift,100*np.array(p),width=.34,color=c,label=label)
    axs[1].set(xticks=x,xticklabels=['Disconnected','2','3','4'],xlabel='Connected host diameter',ylabel='Original-preparation probability (%)',title='Six sites: diameter distribution')
    axs[1].tick_params(axis='x',labelsize=8);axs[1].legend(frameon=False,fontsize=8)
    for n,b,c in [(5,256,colors[0]),(6,192,colors[1])]:
        d=max((json.loads(p.read_text()) for p in P.glob(f'HOST_RESPONSE_N{n}_Q{n+1}_B{b}_L*.json')),key=lambda d:d['level'])
        r=next(r for r in d['rows'] if r['lag']==.25)
        ss=sorted(r['shells'].items(),key=lambda v:int(v[0]));xx=np.array([int(k) for k,v in ss])
        mid=np.array([abs(v['response']['value']) for k,v in ss]);I=np.array([v['response']['absolute_interval'] for k,v in ss])
        ok=I[:,0]>0
        axs[2].plot(xx[ok],mid[ok],'-o',color=c,label=f'N={n}',ms=4)
        if ok.any():axs[2].errorbar(xx[ok],mid[ok],yerr=np.stack([mid[ok]-I[ok,0],I[ok,1]-mid[ok]]),fmt='none',color=c,capsize=3)
        if (~ok).any():axs[2].errorbar(xx[~ok],I[~ok,1],yerr=I[~ok,1]*.3,uplims=True,fmt='v',color=c,ms=4)
    axs[2].set(yscale='log',xticks=[1,2,3,4],xlabel='Initial active-host distance',ylabel=r'$|\chi_{0j}|$ at lag 1/4',title='Host-only source response')
    axs[2].legend(frameon=False,fontsize=8)
    for ax in axs:ax.grid(axis='y',alpha=.17)
    fig.suptitle('Frozen CP-12 • auxiliary Q=N+1 • complete host record costs reported separately',fontsize=12)
    output=P/'ACTIVE_HOST.svg';fig.savefig(output,metadata={'Date':None})
    output.write_text('\n'.join(line.rstrip() for line in output.read_text().splitlines())+'\n')
    if a.png:fig.savefig(a.png,dpi=160)


if __name__=='__main__':main()
