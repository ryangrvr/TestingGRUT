"""Scientific plots from measured diagnostics, never generated imagery."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
P=Path(__file__).resolve().parent;R=P/'raw'
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','axes.titleweight':'bold','svg.fonttype':'none'})
def read(name):return np.genfromtxt(R/(name+'.csv'),delimiter=',',names=True,ndmin=1)
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
colors=['#8cb4d1','#397fb0','#143f66']
ax=axes[0,0]
for tag,col in zip(['coarse','medium','fine'],colors):
 a=read('single_generic_'+tag);ax.plot(a['t'],a['c0_Q12']/a[0]['Q'],label=tag,color=col)
ax.set(xlabel=r'$\tau$',ylabel='Core charge / initial total',ylim=(.993,1.0003),title='Perturbed isolated object · unrestricted 3D')
ax.text(.03,.08,'R = 12 sphere; 90% criterion lies below axis',transform=ax.transAxes,fontsize=9)
ax.legend(frameon=False,loc='lower right')
ax=axes[0,1]
for kind,ls in [('reference','-'),('generic','--')]:
 a=read('pair_'+kind+'_fine');m=a['t']>=120
 for i,c in [(0,'#17628f'),(1,'#c46c24')]:ax.plot(a['t'][m],a[f'c{i}_Q12'][m]/a[0]['Q'],color=c,ls=ls,label=f'{kind}, '+('smaller core' if i==0 else 'larger core'))
ax.set(xlabel=r'$\tau$',ylabel='Each core charge / initial total',title='Outgoing cores · finest Cartesian grid',ylim=(.27,.69))
ax.legend(frameon=False,fontsize=8)
ax=axes[1,0]
for tag,col in zip(['coarse','medium','fine'],colors):
 a=read('impact_'+tag);a=a[a['t']>=120];ax.plot(a['t'],a['c0_x']/.02,color=col,label=tag);ax.plot(a['t'],a['c1_x']/.02,color=col,ls='--')
ax.set(xlabel=r'$\tau$',ylabel='Core x / per-source offset δ',title='Matched transverse seed · full nonlinear 3D')
ax.text(.03,.08,'Solid: smaller core. Dashed: larger core.',transform=ax.transAxes,fontsize=8)
ax.legend(frameon=False,fontsize=8,loc='center left')
ax=axes[1,1]
g=json.loads((P/'GEOMETRIC_COMPARISON.json').read_text());mean=g['target']['mean'];sd=g['target']['sigma'];lower=g['H0_conservative_lower_bound']
ax.errorbar(mean,1,xerr=sd,fmt='o',capsize=5,color='#17628f',label='Megamaser measurement (68%)')
ax.plot([mean-3*sd,mean+3*sd],[1,1],lw=3,alpha=.3,color='#17628f')
ax.annotate('',xy=(133,0),xytext=(lower,0),arrowprops={'arrowstyle':'->','color':'#b34738','lw':3})
ax.plot(lower,0,'|',ms=18,mew=2,color='#b34738')
ax.text(95,.53,'Incompatible bounds',ha='center',color='#b34738',weight='bold')
ax.text(lower,-.22,f'≥ {lower:.2f}',ha='center',color='#b34738')
ax.set(xlim=(58,138),ylim=(-.45,1.45),yticks=[0,1],yticklabels=['Specified GRUT preparation','Geometric H₀ data'],xlabel=r'$H_0$ [km s$^{-1}$ Mpc$^{-1}$]',title='Conserved-charge empirical exclusion')
ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
fig.suptitle('Finite-object evidence and an empirical test have different outcomes',fontsize=15,weight='bold')
fig.savefig(P/'RESULTS_FIGURE.png',dpi=170);fig.savefig(P/'RESULTS_FIGURE.svg')
plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
for ax,typ in zip(axes,['cartesian','original']):
 for tag,h,col in zip(['coarse','medium','fine'],[.5,1/3,.25],colors):
  if typ=='original':
   d=json.loads((R/('large_modes_h'+format(h,'.6f')+'.json')).read_text());t=np.array([r['time'] for r in d['trace']]);v=np.array([r['relative_symmetry_residual'] for r in d['trace']])
  else:
   a=read('pair_reference_'+tag);t=a['t'];v=np.array([a['res0_norm']/a['gen0_norm'],a['res1_norm']/a['gen1_norm']]).T
  m=t>0;ax.semilogy(t[m],v[m,0],color=col,label=tag);ax.semilogy(t[m],v[m,1],color=col,ls='--')
 ax.set(xlabel=r'$\tau$',ylabel='Relative generator residual',ylim=(1e-4,3),title='Cartesian · R36 diagnostic' if typ=='cartesian' else 'Original m=1 solver · large domain')
 ax.legend(frameon=False,fontsize=9)
fig.suptitle('Spatial symmetry errors remain material · solid translation, dashed rotation',fontsize=13,weight='bold')
fig.savefig(P/'SYMMETRY_FIGURE.png',dpi=170);fig.savefig(P/'SYMMETRY_FIGURE.svg');plt.close(fig)

fig,axes=plt.subplots(1,3,figsize=(12,4.7),layout='constrained')
volumes={}
for name in ['single_generic_fine','pair_generic_fine','production_fine']:
 meta=json.loads((R/(name+'.csv.spatial.json')).read_text());a=np.fromfile(R/(name+'.csv.spatial.f32'),dtype=np.float32).reshape(meta['n'],meta['n'],meta.get('nz',meta['n']),3);volumes[name]=(meta,a)
norm=LogNorm(vmin=1e-4,vmax=max(4,max(float(a[...,0].max()) for meta,a in volumes.values())))
for ax,name,title in zip(axes,['single_generic_fine','pair_generic_fine','production_fine'],['Perturbed isolated core · τ160','Perturbed outgoing cores · τ160','Condensate fragmentation · τ240']):
 meta,a=volumes[name];n=meta['n'];nz=meta.get('nz',n);Lz=meta.get('Lz',meta['L']);h=meta['h'];st=meta['stride']
 # The stored volume is subsampled from the actual grid, not an evolved 2D slice.
 y=np.argmin(abs(-meta['L']+(np.arange(n)*st+.5)*h));F=a[:,y,:,0]
 xc=-meta['L']+(np.arange(n)*st+.5)*h;zc=-Lz+(np.arange(nz)*st+.5)*h
 im=ax.pcolormesh(xc,zc,F.T,shading='nearest',cmap='magma',norm=norm);ax.set_aspect('equal');ax.set_xlim(-meta['L'],meta['L']);ax.set_ylim(-Lz,Lz)
 ax.set(xlabel='x',ylabel='z',title=title)
 ax.text(.02,.02,f"y≈{xc[y]:.2f}; displayed stride {st}",transform=ax.transAxes,color='white',fontsize=8)
fig.colorbar(im,ax=axes,label=r'$u^2+w^2$ (central y slice)',shrink=.8)
fig.suptitle('Spatial diagnostics of full 3D evolutions',fontsize=14,weight='bold')
fig.savefig(P/'SPATIAL_FIGURE.png',dpi=170);fig.savefig(P/'SPATIAL_FIGURE.svg')
for name in ['RESULTS_FIGURE.svg','SPATIAL_FIGURE.svg','SYMMETRY_FIGURE.svg']:
 path=P/name;path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
