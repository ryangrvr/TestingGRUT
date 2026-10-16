"""Deterministic scientific figure from saved S2 evidence, no generated imagery."""
from pathlib import Path
import json
import io
import tempfile
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent


def main():
    r=json.loads((P/'RESULTS.json').read_text())
    plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans',
        'svg.hashsalt':'wm2-s2-storage-20261016','axes.spines.top':False,
        'axes.spines.right':False,'axes.titlesize':13,'axes.titleweight':'bold'})
    fig,axes=plt.subplots(1,3,figsize=(16,5.4))
    blue='#176b91';gold='#b36a12';red='#aa3745'
    f=r['formation'];a=axes[0]
    a.plot([x['N'] for x in f],[x['finite_window_macroscopic_probability'] for x in f],
        'o-',color=blue,label='Radial finite-window calculation')
    a.axhline(.375,color=gold,linestyle='--',label='All-size lower bound: 3/8')
    a.set_xscale('log',base=2);a.set_ylim(0,.85);a.set_xlabel('Number of events N')
    a.set_ylabel('Probability of attaching at least half')
    a.set_title('A  Formation at exact resonance',loc='left')
    a.legend(fontsize=9,loc='lower left',frameon=False)
    a.text(.02,.94,'One uniform observation in [0, T_N]\nT_N grows as N log N',
        transform=a.transAxes,va='top',fontsize=10,color='#444444')
    a=axes[1];ns=np.unique(np.geomspace(8,32768,120).astype(int))
    from model import detuning_bound
    for d,col in [(.01,blue),(.1,gold),(1.,red)]:
        a.loglog(ns,[detuning_bound(int(n),d)['all_time_macroscopic_upper'] for n in ns],
            color=col,label=f'Fixed detuning {d:g}')
    a.set_xlabel('Number of events N');a.set_ylabel('All-time formation upper bound')
    a.set_title('B  Detuning restores the obstruction',loc='left')
    a.set_ylim(1e-9,1.3);a.legend(fontsize=9,frameon=False,loc='lower left')
    a.text(.98,.94,'Analytic bounds, not sampled maxima',ha='right',va='top',
        transform=a.transAxes,fontsize=9,color='#444444')
    a=axes[2];ds=np.arange(1,13);c=r['constants']
    a.semilogy(ds,[c['Sigma_C0']*c['Sigma_ratio']**int(d) for d in ds],
        '--',color=gold,label='Analytic bound, all N')
    a.set_xlabel('Generated operational distance d')
    a.set_ylabel('All-record virtual operator norm')
    a.set_title('C  Protection includes every fuel state',loc='left')
    a.set_xticks([1,3,5,7,9,11]);a.legend(fontsize=9,frameon=False,loc='upper right')
    a.text(.98,.76,'Bound: (1/680)(6/91)^d\nE <= 2.5; gap >= 42.5',ha='right',
        transform=a.transAxes,fontsize=10,color='#444444')
    checked=sum(len(x['rows']) for x in r['Schur'])
    a.text(.03,.08,f'{checked} all-record matrix blocks checked\nNo fuel-state truncation',
        transform=a.transAxes,fontsize=10,color='#444444')
    for a in axes:a.grid(alpha=.16,which='major')
    fig.suptitle('WM2-S2: a resonant construction, not a robust spatial phase',
        x=.06,ha='left',fontsize=17,fontweight='bold')
    fig.text(.06,.025,'Separate branching control: five protected moves cross growing raw tree distance and restore the whole tree.  '
        'No selected dimension or empirical prediction.',fontsize=10,color='#444444')
    fig.tight_layout(rect=(.035,.075,1,.9),w_pad=2.5)
    target=P/'S2_RESULTS.svg';buffer=io.StringIO()
    fig.savefig(buffer,format='svg',metadata={'Date':None},bbox_inches='tight')
    text='\n'.join(line.rstrip() for line in buffer.getvalue().splitlines())+'\n'
    ET.fromstring(text)
    temporary=target.with_suffix('.svg.tmp');temporary.write_text(text);temporary.replace(target)
    fig.savefig(Path(tempfile.gettempdir())/'wm2_s2_results_preview.png',dpi=140,bbox_inches='tight')


if __name__=='__main__':main()
