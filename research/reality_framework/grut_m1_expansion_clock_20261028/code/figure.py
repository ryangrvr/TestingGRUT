"""Plot two supplied source windows under the same charge-dependent model."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicHermiteSpline

ROOT=Path(__file__).resolve().parents[1]
RF=ROOT.parent
ref=187479.9981052459
cases=[('profiled96_to10','Original source through $1000/m_X$','#187f83','parent'),
       ('profiled_h0025_96_to10','Finer source through $10000/m_X$','#d56c36','h0025')]
fig,axs=plt.subplots(2,2,figsize=(11.6,7.6),layout='constrained')
fig.patch.set_facecolor('#f7f9f8')
for ax in axs.flat:
    ax.set_facecolor('#ffffff');ax.grid(color='#d7dfdc',alpha=.7,lw=.6);ax.spines[['top','right']].set_visible(False)
for tag,label,color,source in cases:
    d=json.loads((ROOT/(tag+'.json')).read_text());r=d['snapshots']
    tau=np.array([x['tau'] for x in r]);valid=tau>0;tau=tau[valid];r=[x for x in r if x['tau']>0]
    N=np.array([x['Nminus'] for x in r]);visible=N>1e-7
    axs[0,0].plot(tau[visible],N[visible],marker='o',ms=2.7,color=color,label=label)
    axs[0,1].plot(tau,[x['Qcore']/2000 for x in r],marker='o',ms=2.7,color=color)
    axs[1,0].plot(tau,[100*x['gas_kinetic_mX']/ref for x in r],marker='o',ms=2.7,color=color)
    cold=np.array([100*x['stationary_cold_excess_fraction'] for x in r]);
    axs[1,1].plot(tau,cold,marker='o',ms=2.7,color=color,label=label)
    pulse=(json.loads((RF/'grut_m1_kinetics_20261028'/'PHYSICS_CHECKS.json').read_text())['source_energy_ledger']['pulse_initial_lattice_energy_mX'] if source=='parent' else json.loads((ROOT/'sources'/f'ringdown_{source}.json').read_text())['initial_pulse_energy_mX'])
    budget=np.array([100*((pulse-x['work_redshift_mX']+x['interpolated_generator_energy_defect_mX'])/ref-1) for x in r])
    axs[1,1].plot(tau,budget,ls='--',lw=1.5,color=color,alpha=.8)
for ax in axs.flat:
    ax.set_xscale('log');ax.set_xlabel('Optical time $\\tau$')
axs[0,0].set_yscale('log');axs[0,0].set(ylabel='Antiparticles / reference core',title='A  Residual antiparticles')
axs[0,1].set(ylabel='Core charge / reference charge (%)',title='B  Core charge retained')
axs[1,0].set(ylabel='Gas kinetic energy / reference mass (%)',title='C  Expansion removes gas heat')
axs[1,1].axhline(1,color='#66737b',lw=1,ls=':')
axs[1,1].set(ylabel='Excess / reference mass (%)',title='D  Cold mass and full energy ledger')
axs[1,1].text(.025,.04,'Solid: cold composition only\nDashed: includes unassigned energy',transform=axs[1,1].transAxes,fontsize=8,color='#3c4746')
axs[0,0].legend(loc='lower left',fontsize=8,frameon=False)
fig.suptitle('M1 supplied-pulse continuation: expanding, quasistatic cores',fontsize=15,fontweight='bold',color='#213a37')
fig.text(.5,-.008,'Illustrative 20 MeV start; last point corresponds to a factor 10 expansion. Not a derived formation or final cosmological abundance.',ha='center',fontsize=8.3,color='#53635f')
fig.savefig(ROOT/'EXPANSION.png',dpi=180,bbox_inches='tight',facecolor=fig.get_facecolor())
fig.savefig(ROOT/'EXPANSION.svg',bbox_inches='tight',facecolor=fig.get_facecolor())

if __name__=='__main__':pass
