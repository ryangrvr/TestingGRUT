"""Compact, explicitly conditional visualization of the saved plume results."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from local_rates import pulse_series

ROOT=Path(__file__).resolve().parents[1]
def read(name):return json.loads((ROOT/(name+'.json')).read_text())

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':10})
fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
fig.patch.set_facecolor('#f8faf9')
blue='#126a78';orange='#d46f38';dark='#334b52'
for ax in axs.flat:
    ax.set_facecolor('white');ax.grid(alpha=.18);ax.spines[['top','right']].set_visible(False)

z=np.load(ROOT/'source_t400.npz');t,e,p,m=pulse_series(z,5.)
for f,color,label in ((p.sum(axis=1)*(e[1]-e[0]),blue,r'$X$'),(m.sum(axis=1),orange,r'$\bar X$')):
    n=cumulative_trapezoid(f,x=t,initial=0);axs[0,0].plot(t,n/n[-1],color=color,lw=2.4,label=label)
axs[0,0].axvline(1000,color=dark,alpha=.6,ls=':')
axs[0,0].set(xlim=(0,10000),ylim=(0,1.04),xlabel=r'$m_Xt$',ylabel='Fraction of each source emitted',title='A  ·  Source times (400-wide window)')
axs[0,0].legend(frameon=False,loc='lower right')

rr=read('SOURCE_DIAGNOSTICS')['matching_radius_controls']
axs[0,1].plot([x['actual_radius'] for x in rr],[x['collinear_primary_NWA_events_outside_surface'] for x in rr],
              color=blue,marker='o',lw=2.4)
axs[0,1].set(xlabel=r'Matching radius [$m_X^{-1}$]',ylabel='Events per core beyond radius',title='B  ·  Collinear test (radial step 0.05)')
axs[0,1].text(.04,.92,'Smaller-radius interior not computed',transform=axs[0,1].transAxes,fontsize=9,va='top',color=dark)

mix=read('mixing_seed37');times=np.array(mix['snapshot_times_mX_inverse']);primary=mix['regions'][0]['primary']
for species,color,label in (('plus',blue,r'$X$'),('minus',orange,r'$\bar X$')):
    p=primary[species];y=[x['arrived_number']/p['initial_number'] for x in p['snapshots']]
    axs[1,0].semilogx(times,y,color=color,marker='o',lw=2.2,label=label)
axs[1,0].set(xlabel=r'$m_Xt$',ylabel='Fraction crossing half-spacing',ylim=(0,1.05),title='C  ·  Ballistic first passage')
axs[1,0].legend(frameon=False,loc='lower right')

shot=read('shot_noise_dr2p5')['rows'];tt=[x['mX_t'] for x in shot]
for key,color,label in (('poisson_cv_plus',blue,r'$X$'),('poisson_cv_minus',orange,r'$\bar X$')):
    axs[1,1].loglog(tt,[x[key] for x in shot],color=color,marker='o',lw=2.2,label=label)
axs[1,1].axhline(1,color=dark,ls=':',alpha=.6)
axs[1,1].set(xlabel=r'$m_Xt$',ylabel='Pointwise density std. dev. / mean',title='D  ·  Poisson-center model')
axs[1,1].legend(frameon=False,loc='upper right')

fig.suptitle('M1 local plume: bounded, collisionless model checks',fontsize=15,fontweight='bold',color='#173942')
fig.savefig(ROOT/'LOCAL_CLOUD.png',dpi=180,facecolor=fig.get_facecolor())
svg=ROOT/'LOCAL_CLOUD.svg'
fig.savefig(svg,facecolor=fig.get_facecolor())
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')

if __name__=='__main__':print('Wrote LOCAL_CLOUD.png and LOCAL_CLOUD.svg')
