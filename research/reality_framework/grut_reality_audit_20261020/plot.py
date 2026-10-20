"""Scientific figure for the separate full-field numerical audit."""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
summary=json.loads((ROOT/'SUMMARY.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'svg.hashsalt':'grut-reality-audit-20261020'})
fig,axes=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
fig.suptitle('Full-field tests support the bulk predictions; empirical validation remains open',fontsize=13)
def trace(name):return np.array(json.loads((ROOT/'results'/(name+'.json')).read_text())['trace'])

row=trace('dense_k02'); selected=row[:,0]<=100
t=row[selected,0];signal=row[selected,3]/1e-5
expected=summary['cases']['dense_k02']['expected']
design=np.column_stack([np.ones(len(t)),*[func(w*t) for w in expected for func in [np.cos,np.sin]]])
coefficients=np.linalg.lstsq(design,signal,rcond=None)[0]
ax=axes[0]
ax.plot(t,signal,label='Full-field evolution',lw=1.4,color='#2563eb')
ax.plot(t[::12],(design@coefficients)[::12],'o',ms=2,mfc='none',color='#d97706',label='Predicted frequencies; fitted amplitudes')
ax.set(title='Dense bulk, k = 0.2',xlabel=r'Dimensionless time $\tau$',ylabel='Charged-amplitude mode / initial seed')
ax.legend(frameon=False,fontsize=8,loc='lower left')
ax=axes[1]
for name,color,label in [('mixed_k02','#16a34a','Mixed bulk, k = 0.2'),('mixed_k04','#dc2626','Mixed bulk, k = 0.4'),('mixed_stable_k08','#6b7280','Above-band control, k = 0.8')]:
    rows=trace(name); initial=rows[0,9]
    ax.plot(rows[:,0],rows[:,9]/initial,color=color,lw=1.5,label=label)
    case=summary['cases'][name]
    if 'expected_growth' in case:
        start,end=case['fit_window']; anchor=rows[np.argmin(abs(rows[:,0]-start))]
        times=np.linspace(start,end,30)
        ax.plot(times,anchor[9]/initial*np.exp(case['expected_growth']*(times-start)),'--',lw=.9,color='#111827')
ax.set(title='Coherent mixed bulk: predicted growth band',xlabel=r'Dimensionless time $\tau$',ylabel='Mode norm / initial mode norm',yscale='log',xlim=(0,160))
ax.legend(frameon=False,fontsize=8,loc='upper left')
ax.text(.03,.04,'Dashed: prior growth prediction, anchored for display',transform=ax.transAxes,fontsize=8)
for ax in axes:ax.grid(alpha=.18)
svg=ROOT/'AUDIT_FIGURE.svg'
fig.savefig(svg,metadata={'Date':None})
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--qa',type=Path,default=ROOT/'AUDIT_FIGURE.png')
qa=parser.parse_args().qa
qa.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(qa,dpi=145)
print(qa)
