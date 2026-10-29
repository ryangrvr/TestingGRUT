"""M1 checkpoint gas collision and viscosity-weighted rate diagnostics.

Uses the parent tree-level X-X amplitude and identical-particle counting.
The viscosity-weighted integral is a proxy, not a solved transport equation.
"""
from pathlib import Path
import sys,json,time
import numpy as np
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parents[1];RF=ROOT.parent
sys.path.insert(0,str(RF/'grut_m1_expansion_clock_20261028'/'code'))
import expansion,transport

def rates(e,pop,nangle):
 mu,wm=leggauss(nangle);ct,wt=leggauss(nangle);mom=np.sqrt(e*e-1);total=viscosity=0.
 for i in range(len(e)):
  j=np.arange(i,len(e));s=2+2*(e[i]*e[j][None,:]-mom[i]*mom[j][None,:]*mu[:,None])
  p2=np.maximum(s/4-1,0);beta=np.sqrt(np.maximum(1-4/s,0));vm=np.sqrt(np.maximum(s*(s-4),0))/(2*e[i]*e[j][None,:])
  t=-2*p2[:,:,None]*(1-ct[None,None,:]);u=-2*p2[:,:,None]*(1+ct[None,None,:])
  amp=4*.001+transport.GH**2*(1/(t-6)+1/(u-6))
  integrand=(wm[:,None,None]/2)*wt[None,None,:]*vm[:,:,None]*amp**2/(64*np.pi*s[:,:,None])
  K=integrand.sum(axis=(0,2));Kv=(integrand*(1-ct[None,None,:]**2)).sum(axis=(0,2))
  w=pop[i]*pop[j];w[0]*=.5;total+=float(w@K);viscosity+=float(w@Kv)
 return total,viscosity

def main():
 rows=[];t0=time.monotonic();epsilon=expansion.clock()['H0_mX2_over_ncore0'];a=10.
 for bins in (64,96):
  z=np.load(RF/'grut_m1_core_matching_20261029'/f'late_kin_h0025_{bins}.npz');e=z['E'];n=len(e);v=z['Y'][:,-1];plus=np.maximum(v[:n],0);N=float(plus.sum())
  for angles in (8,16):
   events,visc=rates(e,plus,angles)
   rows.append({'shift_bins':bins,'kinetic_nodes':n,'incident_angle_nodes':angles,'cm_angle_nodes':angles,'plus_per_core':N,
     'events_per_core_per_density_tau':events,'viscosity_weighted_events_per_density_tau':visc,
     'total_collisions_per_particle_per_Hubble':2*events/(N*a*epsilon),
     'viscosity_weighted_collisions_per_particle_per_Hubble':2*visc/(N*a*epsilon),
     'viscosity_to_total_ratio':visc/events})
 original=1.62913547338985;ref=rows[0]['events_per_core_per_density_tau']
 out={'a_over_a0':10,'epsilon':epsilon,'rows':rows,'elapsed_seconds':time.monotonic()-t0,
      'reference_full_parent_operator_events64':original,'reference_relative_error':ref/original-1,
      'scope':'At one instant only: homogeneous populations, the parent tree-level elastic X-X kernel and its symmetric identical-event counting. Viscosity weighting sin²(theta_CM) is an angular-transfer proxy, not a computed perturbation relaxation eigenvalue. Full formation, finite-density screening, gas-core interactions and later evolution remain conditional.'}
 (ROOT/'collision_rate.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
