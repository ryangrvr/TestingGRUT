"""Build derived comparison tables and the scientific summary figure."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent
def read(name):return json.loads((ROOT/(name+'.json')).read_text())
def main():
 gas=read('elastic_s96_n512_a16');h=read('holder_model');p=read('perturbations')
 previous=ROOT/'CLASS_CONTROL_RESULTS.json';prior={x['tag']:x for x in json.loads(previous.read_text())['attempts']} if previous.exists() else {}
 cls={'scope':'CLASS collisionless, adiabatic back-extrapolation controls; no formation matching or self-scattering; not the M1 spectrum.','rows':[],'attempts':[]}
 for tag,cold in [('warmrk20','coldrk20'),('warmrk120','coldrk_refined')]:
  a=np.load(ROOT/(tag+'.npz'))['spectrum'];b=np.load(ROOT/(cold+'.npz'))['spectrum']
  for k in [1.,10.]:
   ratio=float(np.exp(np.interp(np.log(k),np.log(a[:,0]),np.log(a[:,1]))-np.interp(np.log(k),np.log(b[:,0]),np.log(b[:,1]))))
   cls['rows'].append({'tag':tag,'k_h_Mpc_inverse':k,'power_ratio':ratio})
 for tag in ['coldmicro','warmmicro','coldrk20','warmrk20','coldrk_refined','warmrk120']:
  j=read(tag);log=ROOT/'.work/class'/(tag+'.log');cls['attempts'].append({'tag':tag,'status':j['status'],'seconds':j['seconds'],'log_tail':log.read_text()[-2200:] if log.exists() else prior.get(tag,{}).get('log_tail','not retained')})
 (ROOT/'CLASS_CONTROL_RESULTS.json').write_text(json.dumps(cls,indent=2)+'\n')
 tags=['linear_refined','linear_more','linear_fit','linear_fit_refined','linear_collisionless','linear_isothermal','linear_rate_half','linear_rate_double']
 table={t:{str(int(x['k_Mpc_inverse'])):x['T2_z0'] for x in read(t)['rows']} for t in tags}
 table['scope']='Conditional BGK test cosmology with supplied formation perturbations; linearly extrapolated z=0 matter power ratio, not a full M1 or nonlinear prediction.'
 (ROOT/'TRANSFER_COMPARISON.json').write_text(json.dumps(table,indent=2)+'\n')
 cand={'name':'M1 cubic-holder continuation, conditional construction','parent_commit':'94f660d0e04e7b40b668ffb316793a575f3e3813','parent_point':'../grut_matched_candidate_20261027/CANDIDATE_M1.json','parent_H_GeV':h['parent_H_GeV'],
  'mu_eff_initial_over_H':h['mu_eff_initial_over_H'],'xi_muR':h['xi_muR'],'alpha_R':h['alpha_R'],'physical_scalar_mass2_over_H2':h['canonical_mass2_over_H2'][1:],
  'charge_yield_relative_to_parent':p['holder_charge_relative_to_original'],'fixed_dimensionless_family_refits':p['fixed_dimensionless_family_calibration'],
  'gas_source':'Inherited conditional 96-bin checkpoint; no re-formation calculation with the new holder','gas_fraction_of_total_matter':read('linear_fit')['background_parameters']['gas_fraction_of_total_matter'],
  'status':'Classical fixed-R scalar potential bounded. Spectator production and conditional gas transport calculated. No UV completion, actual CDI spectrum, full formation ensemble, or full interacting cosmological prediction.'}
 (ROOT/'STABILIZED_CANDIDATE.json').write_text(json.dumps(cand,indent=2)+'\n')
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
 fig,axes=plt.subplots(1,3,figsize=(14,4.6),layout='constrained')
 for tag,label,color in [('reference','Original preparation','#7c858d'),('holder','Cubic curvature holder','#136b80')]:
  z=np.load(ROOT/(tag+'.npz'));v=z['obs'];keep=v[:,0]>=.1;axes[0].plot(v[keep,0],-v[keep,4],label=label,color=color,lw=1.5)
 axes[0].set(xscale='log',xlabel=r'$H_i t$',ylabel=r'$-N_X$',title='Spectator charge production');axes[0].legend(frameon=False,fontsize=8);axes[0].text(.97,.15,'Late charge: +1.81%',ha='right',transform=axes[0].transAxes)
 z=np.load(ROOT/'elastic_s96_n512_a16.npz');axes[1].loglog(z['T_MeV'],z['viscosity_proxy_over_H'],color='#136b80',label='Viscosity-rate proxy');axes[1].axhline(1,color='#8b8b8b',ls='--',lw=1);axes[1].invert_xaxis();axes[1].set(xlabel='Photon temperature [MeV]',ylabel=r'$\Gamma_{\rm visc}/H$',title='Homogeneous elastic transport');axes[1].text(.06,.09,f"{gas['collisions_per_particle_integrated']:.2f} collisions per particle\nCrossing: {gas['viscosity_proxy_crossing_temperature_MeV']:.3f} MeV",transform=axes[1].transAxes)
 for tag,label,color,ls in [('linear_fit','Interacting fit','#136b80','-'),('linear_collisionless','Collisionless control','#c66d3e','--'),('linear_isothermal','Isothermal formation seed','#767b4d',':')]:
  rows=read(tag)['rows'];axes[2].semilogx([x['k_Mpc_inverse'] for x in rows],[x['T2_z0'] for x in rows],label=label,color=color,ls=ls,lw=1.6)
 axes[2].set(xlabel=r'$k$ [Mpc$^{-1}$]',ylabel=r'$P_{\rm mix}/P_{\rm cold}$',title='Conditional BGK test cosmology',ylim=(.74,1.01));axes[2].legend(frameon=False,fontsize=8,loc='lower left')
 fig.suptitle('GRUT–M1: computed extensions and their stated domains',fontsize=14)
 fig.savefig(ROOT/'summary.png',dpi=180);plt.close(fig)
 print(json.dumps({'CLASS':cls['rows'],'figure':'summary.png'},indent=2))
if __name__=='__main__':main()
