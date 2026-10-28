"""Scientific figure: computed conditional shapes and fixed-yield mass ledger."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
o=json.loads((ROOT/'OBSERVABLE_PREDICTION.json').read_text())
a=json.loads((ROOT/'ABUNDANCE_LEDGER.json').read_text())
s=np.load(ROOT/'RECOIL_SHAPES.npz')
edges=s['edges_keV'];centers=(edges[:-1]+edges[1:])/2;width=np.diff(edges)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(10.8,7.2),layout='constrained')
for ax,name in zip(axes.flat,['Xe131','Ar40','Ge74']):
    ax.semilogy(centers,s[name+'_object_probability']/width,color='#235789',lw=2,label='Supplied object profile')
    ax.semilogy(centers,s[name+'_point_probability_at_joint_mass']/width,color='#c45a22',ls='--',label='Common point-WIMP fit')
    ax.set(title=name,xlabel='Recoil energy (keV)',ylabel='Normalized density (1/keV)',xlim=(5,40),ylim=(1e-7,1))
    zero=next(z for z in o['first_zero_target_energies'] if z['target']==name)
    if zero['inside_primary_window']:
        ax.axvline(zero['ER_keV'],color='#32865d',ls=':',lw=1.4)
        ax.text(.48,.85,'Profile zero: 26.7 keV',transform=ax.transAxes,color='#286b4b',fontsize=9)
axes[0,0].legend(loc='lower left',fontsize=8,frameon=False)
ax=axes[1,1]
eta=np.linspace(0,1,401);m=a['mX_GeV'];ep=a['object_energy_per_charge_GeV'];Y=a['YQ_input'];D=a['target_DM_per_visible_entropy_GeV']
needed=(D-Y*(eta*ep+(1-eta)*m))/D*100
lo,hi=a['core_charge_fraction_interval_at_1pct_total_DM_pairs']
ax.plot(eta,needed,color='#235789',lw=2)
ax.axhspan(0,1,color='#e5eee9');ax.axvspan(lo,hi,color='#32865d',alpha=.16)
ax.axhline(1,color='#c45a22',ls='--',lw=1.2)
ax.axhline(0,color='#777777',lw=.8)
ax.set(title='Cold mass consistency at the fixed yield',xlabel='Supplied fraction of charge in reference objects η',
       ylabel='Pairs needed to fill target mass (%)',xlim=(0,1),ylim=(-2,8))
ax.text(.39,.76,'All charge in objects:\n6.35% pairs needed',transform=ax.transAxes,fontsize=9)
ax.text(.38,.12,'Allowance: at most 1% of fixed DM target',transform=ax.transAxes,fontsize=8)
fig.suptitle('Conditional calculations: no detector response or physical final population',fontsize=12)
fig.savefig(ROOT/'FIVE_OUTPUTS.png',dpi=170)
plt.close(fig)
print('Wrote FIVE_OUTPUTS.png from the calculated results.')
