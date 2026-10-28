"""Scientific summary of the actual changed potential, profile and calibration."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent
d=json.loads((ROOT/'DEFORMATION_RESULTS.json').read_text());m=json.loads((ROOT/'DEFORMED_MATCHED_RESULTS.json').read_text());old=json.loads((ROOT/'MATCHED_RESULTS.json').read_text())
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(2,2,figsize=(10.6,7.3),layout='constrained')
u=np.linspace(-1.5,1.5,601);sym=.03/4*(u*u-1)**2;changed=sym-.002*(u-1)**3
ax[0,0].plot(u,sym,'--',color='#777777',label='Symmetric potential')
ax[0,0].plot(u,changed,color='#235789',lw=2,label='M1: one stationary vacuum')
ax[0,0].set(xlabel='Binding field φ/v',ylabel='Vφ / v⁴',title='One supplied deformation removes the second vacuum',ylim=(-.001,.055))
ax[0,0].legend(frameon=False,fontsize=8)
rows=[x for x in d['rows'] if x['tag'].startswith('grid_')]
ax[0,1].plot([x['epsilon'] for x in rows],[x['E_over_Qm'] for x in rows],'o-',color='#235789')
ax[0,1].axhline(1,color='#c45a22',ls='--');ax[0,1].axvline(d['single_stationary_vacuum_threshold'],color='#32865d',ls=':')
ax[0,1].axvspan(d['single_stationary_vacuum_threshold'],.003,color='#32865d',alpha=.12)
ax[0,1].set(xlabel='New parameter ε',ylabel='Object energy / free constituent mass',title='The tested objects stay bound',ylim=(.91,1.01))
new=np.loadtxt(ROOT/'PROFILE_M1.csv',delimiter=',',skiprows=1)
prior=np.loadtxt(RF/'grut_portal_debris_window_20261023/data/profile_Q200000_lambda0.03_refined.csv',delimiter=',',skiprows=1)
for col,color,label in ((1,'#235789','Charged amplitude F'),(3,'#32865d','Binding field S')):
    ax[1,0].plot(prior[:,0],prior[:,col],ls='--',color=color,alpha=.6)
    ax[1,0].plot(new[:,0],new[:,col],color=color,label=label)
ax[1,0].set(xlabel='Radius × mX',ylabel='Dimensionless field',xlim=(0,30),title='New profiles; dashed curves are the old potential')
ax[1,0].legend(frameon=False,fontsize=8,loc='upper right')
base=next(x for x in old['rows'] if x['tag']=='b3_g0.5_baseline')
central=next(x for x in m['rows'] if x['tag']=='M1_b3_g0.5_matched')
D=m['target_DM_per_visible_entropy_GeV'];ep_old=old['object_energy_per_charge_GeV'];ep_new=m['object_energy_per_charge_GeV']
ratios=[base['absolute_yield_magnitude']*ep_old/D,base['absolute_yield_magnitude']*ep_new/D,central['target_mass_ratio']]
ax[1,1].plot(range(3),ratios,'o-',color='#235789',lw=1.5)
ax[1,1].axhline(1,color='#32865d',ls='--')
for i,y in enumerate(ratios):ax[1,1].annotate(f'{y:.5f}',(i,y),xytext=(0,9),textcoords='offset points',ha='center',fontsize=9)
ax[1,1].set(xticks=range(3),xticklabels=['Old\npotential','One-vacuum\nunscaled','M1\nmatched'],ylabel='Conditional mass / target D',title='Same source, changed object mass, recomputed scale',ylim=(.965,1.012),xlim=(-.25,2.25))
fig.suptitle('GRUT-M1: a tested conditional candidate; final formation is still to compute',fontsize=12)
fig.savefig(ROOT/'M1.png',dpi=170);plt.close(fig)
print('Wrote M1.png.')
