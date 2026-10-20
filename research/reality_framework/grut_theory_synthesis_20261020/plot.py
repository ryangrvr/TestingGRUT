"""Plot the derived classical phase structure; all axes are dimensionless."""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
result=json.loads((ROOT/'DERIVATION_RESULTS.json').read_text())
wall=json.loads((ROOT/'wall_profile.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'svg.hashsalt':'grut-theory-20261020'})
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
fig.suptitle('Candidate 1: a derived binding regime and tests of its phase structure',fontsize=15)
a,b=result['a'],result['b']
mu=result['benchmark']['omega0']
f0=result['benchmark']['F0']
f=np.linspace(0,2.5,600)
v=np.where(f*f<=a,.5*f*f+.25*(b-1/a)*f**4,.25*a+.25*b*f**4)
ax=axes[0,0]
ax.plot(f,v,label='Local potential minimum',color='#2563eb')
ax.plot(f,.5*mu*mu*f*f,'--',label=r'Sharp bound $\mu_0^2 F^2/2$',color='#d97706')
ax.scatter([f0],[.5*mu*mu*f0*f0],s=30,color='#111827',zorder=4)
ax.set(title='A. Charge binding is allowed when ab < 1',xlabel='Charged amplitude F',ylabel='Potential density')
ax.legend(frameon=False,fontsize=8)
ax=axes[0,1]
ax.plot(wall['x'],wall['F'],label='F',color='#2563eb')
ax.plot(wall['x'],wall['S'],label='S',color='#d97706')
ax.set(title='B. Interface between charged bulk and vacuum',xlabel='Distance across interface',ylabel='Field amplitude',xlim=(-8,8))
ax.text(.04,.08,r'$\sigma\simeq0.68001349$',transform=ax.transAxes)
ax.legend(frameon=False)
ax=axes[1,0]
disp=result['dispersion']
k=np.array([r['k'] for r in disp])
for key,label,color in [('nu_minus_squared','Phase / sound','#2563eb'),('nu_plus_squared','Amplitude','#d97706'),('nu_S_squared','Mediator','#16a34a')]:
    ax.plot(k,np.sqrt(np.maximum(0,[r[key] for r in disp])),label=label,color=color)
ax.plot(k,k/np.sqrt(3),'--',color='#6b7280',lw=1,label=r'$k/\sqrt{3}$ limit')
ax.set(title='C. Full linear spectrum of uniform dense bulk',xlabel='Wave number k',ylabel='Frequency',xlim=(0,2),ylim=(0,3))
ax.legend(frameon=False,fontsize=8)
ax=axes[1,1]
mixed=result['mixed_branch_prediction']
rows=mixed['growth_scan']
k=np.array([r['k'] for r in rows])
ax.plot(k,[r['growth'] for r in rows],color='#dc2626',label='Full coupled linear equations')
ax.plot(k,k*np.sqrt(-mixed['long_wave_cs_squared']),'--',color='#6b7280',label='Long-wave growth limit')
ax.axvline(mixed['instability_cutoff'],color='#111827',lw=.8,ls=':')
ax.set(title='D. A coherent dilute preparation is unstable',xlabel='Wave number k',ylabel='Growth rate',xlim=(0,1),ylim=(-.002,.14))
ax.text(.03,.08,'F = 0.5; stationary mixed branch\nInfinite uniform medium; no cosmology inferred',transform=ax.transAxes,fontsize=8)
ax.legend(frameon=False,fontsize=8)
for ax in axes.flat:ax.grid(alpha=.18)
svg=ROOT/'THEORY_FIGURE.svg'
fig.savefig(svg,metadata={'Date':None})
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--qa',type=Path,default=ROOT/'THEORY_FIGURE.png')
qa=parser.parse_args().qa
qa.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(qa,dpi=140)
print(qa)
