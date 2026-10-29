"""Visual summary of the late-source and spectral-sensitivity results."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
def read(name):return json.loads((ROOT/(name+'.json')).read_text())

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.5,'axes.titlesize':11})
fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
fig.patch.set_facecolor('#f7faf9');blue='#126a78';orange='#d57438';gray='#8b9a9b';dark='#28434a'
for ax in axes.flat:
    ax.set_facecolor('white');ax.grid(alpha=.18);ax.spines[['top','right']].set_visible(False)

d=read('late_h0025_analysis');rows=d['earlier_core_checkpoints']+d['core_checkpoints'];t=np.array([x['mX_t_absolute'] for x in rows])
y=np.array([x['excess_inside20_over_same_grid_stationary_mX'] for x in rows]);u=np.linspace(t[0],t[-1],300)-t[0]
f=d['phenomenological_fits'];a=f['zero_floor'];b=f['free_nonnegative_floor']
axes[0,0].plot(t,y,'o',color=blue,label='Same-grid core excess')
axes[0,0].plot(u+t[0],a['amplitude_at_4000_mX']*np.exp(-u/a['decay_time_mX_inverse']),color=blue,lw=1.8,ls='--',label='Zero-floor fit')
axes[0,0].plot(u+t[0],b['floor_mX']+b['amplitude_at_4000_mX']*np.exp(-u/b['decay_time_mX_inverse']),color=orange,lw=1.8,ls=':',label='Free-floor fit')
axes[0,0].set(xlabel=r'$m_Xt$',ylabel=r'Inside-$20$ excess [$m_X$]',title='A  ·  Excitation persists through 20,000')
axes[0,0].legend(frameon=False,fontsize=8,loc='upper right')

src=read('late_h0025_t400')['species'];keys=['plus','minus','neutral'];names=[r'direct $X$',r'direct $\bar X$',r'neutral $h$']
vals=[src[k]['positive_projection_number'] for k in keys]
bars=axes[0,1].bar(names,vals,color=[blue,orange,gray],width=.65)
axes[0,1].set_yscale('log');axes[0,1].set_ylim(.5,120)
for x,v in zip(bars,vals):axes[0,1].text(x.get_x()+x.get_width()/2,v*1.12,f'{v:.2f}',ha='center',fontsize=9)
axes[0,1].set(ylabel='Positive-projected source per core',title=r'B  ·  New source: $10^4 < m_Xt < 2\times10^4$')

dis=read('mixing_discrete_sweep')['rows'];cont=read('mixing_cont_nv8001')['rows'];ts=[x['mX_t'] for x in dis]
axes[1,0].loglog(ts,[x['plus']['density_cv'] for x in dis],color=blue,ls='--',marker='o',lw=1.5,label=r'$X$, bin centers')
axes[1,0].loglog(ts,[x['minus']['density_cv'] for x in dis],color=orange,ls='--',marker='o',lw=1.5,label=r'$\bar X$, bin centers')
axes[1,0].loglog([x['mX_t'] for x in cont],[x['density_cv'][0] for x in cont],color=blue,lw=2,label=r'$X$, continuous bins')
axes[1,0].loglog([x['mX_t'] for x in cont],[x['density_cv'][1] for x in cont],color=orange,lw=2,label=r'$\bar X$, continuous bins')
axes[1,0].axhline(.001,color=dark,ls=':',alpha=.7)
axes[1,0].set(xlabel=r'$m_Xt$',ylabel='Density std. dev. / mean',title='C  ·  Poisson model depends on bin filling')
axes[1,0].legend(frameon=False,fontsize=7.5,loc='lower left')

v=read('SCENARIO_AUDIT')['energy_decompositions_at_scale_factor_10'];labels=['through 10,000','through 20,000'];data=[v['parent'],v['through_20000']]
bottom=np.zeros(2)
for key,color,label in [('cold_excess_percent',blue,'Cold excess'),('gas_kinetic_percent',orange,'Gas kinetic'),('held_remainder_percent',gray,'Held remainder')]:
    value=np.array([x[key] for x in data]);axes[1,1].bar(labels,value,bottom=bottom,color=color,label=label,width=.58);bottom+=value
for x,y in enumerate(bottom):axes[1,1].text(x,y+.015,f'{y:.4f}%',ha='center',fontsize=9)
axes[1,1].axhline(1,color=dark,ls=':',alpha=.7);axes[1,1].set_ylim(0,1.08)
axes[1,1].set(ylabel='Excess energy / reference mass (%)',title='D  ·  Conditional expansion ledger')
axes[1,1].legend(frameon=False,fontsize=8,loc='lower center',bbox_to_anchor=(.5,-.28),ncol=3)
fig.suptitle('M1 late-source continuation and matching controls',fontsize=15,fontweight='bold',color='#173942')
fig.savefig(ROOT/'LATE_CORE.png',dpi=180,facecolor=fig.get_facecolor())
svg=ROOT/'LATE_CORE.svg';fig.savefig(svg,facecolor=fig.get_facecolor());svg.write_text('\n'.join(x.rstrip() for x in svg.read_text().splitlines())+'\n')
print('Wrote LATE_CORE.png and LATE_CORE.svg')
