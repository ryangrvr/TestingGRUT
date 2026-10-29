"""Compact visualization of the supplied-pulse controls; no extrapolations."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
j=json.loads((ROOT/'results.json').read_text());d=np.load(ROOT/'field_diagnostics.npz')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','figure.facecolor':'white'})
fig,ax=plt.subplots(1,3,figsize=(13.5,4.4),gridspec_kw={'width_ratios':[1,1.05,1.25]})
colors=['#a2a9bc','#617da7','#1a5773'];hs=[.1,.05,.025]
y=[j['early'][str(h)]['20']['pressure_integral_mean_mX'] for h in hs]
ax[0].axhline(0,color='#44546a',lw=.9);ax[0].bar(range(3),y,color=colors,width=.62)
ax[0].set_ylim(-1.2,1.05)
ax[0].set_xticks(range(3),[str(h) for h in hs]);ax[0].set_xlabel('Radial spacing  h');ax[0].set_ylabel(r'Mean $\int_{r<20}p\,d^3x$  [$m_X$]');ax[0].set_title('Core pressure: t = 100–1000')
for i,v in enumerate(y):ax[0].annotate(f'{v:+.2f}',(i,v),xytext=(0,6 if v>=0 else -15),textcoords='offset points',ha='center',fontsize=9)
for i,name in enumerate(('axis_spec005','axis_spec005_t1000')):
    modes=json.loads((ROOT/(name+'.json')).read_text())['modes'];bottom=0
    for l,c in [(0,'#a9bdcf'),(2,'#236b81'),(4,'#e8ad63')]:
        v=next(q['raw_neutral_energy'] for q in modes if q['ell']==l)
        ax[1].bar(i,v,bottom=bottom,color=c,width=.55,label=f'ℓ = {l}' if i==0 else None);bottom+=v
    ax[1].text(i,bottom+.14,f"{100*next(q['raw_neutral_energy'] for q in modes if q['ell']==2)/bottom:.1f}% ℓ=2",ha='center',fontsize=9)
ax[1].set_xticks([0,1],['0–300','0–1000']);ax[1].set_xlabel(r'Observation window  $m_X t$');ax[1].set_ylabel(r'Neutral energy crossing $r=60$  [$m_X$]');ax[1].set_title('Axisymmetric 5% pulse');ax[1].legend(frameon=False,fontsize=8,loc='upper left')
base=d['open_2000_history'];
for name,r,c in [('reflect_80',80,'#d17e3c'),('reflect_120',120,'#9365a3'),('reflect160_retry',160,'#28788b')]:
    a=d[name+'_history'];ax[2].plot(a['tau'],100*(a['E20']-base['E20']),label=f'R = {r}',lw=1.1,color=c)
ax[2].axhline(0,color='#44546a',lw=.7);ax[2].set_xlabel(r'Time  $m_X t$');ax[2].set_ylabel(r'Change in core $E_{r<20}$  [$m_X$]');ax[2].set_title('Coherent reflecting-cell return');ax[2].legend(frameon=False,fontsize=8)
fig.suptitle('M1 field diagnostics: controlled tests, not a cosmological prediction',fontsize=13,fontweight='bold',y=.99)
fig.text(.5,.015,'Pressure is a finite-region, finite-time mean.  Reflecting cells are 48 million–384 million times denser than the assumed core density.',ha='center',fontsize=9,color='#46576a')
fig.tight_layout(rect=[0,.055,1,.94]);fig.savefig(ROOT/'diagnostics.png',dpi=160)
