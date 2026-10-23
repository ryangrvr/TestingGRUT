"""Portal production/decay and specified hot/cold cosmological histories."""
import json
import numpy as np
from scipy.integrate import cumulative_simpson,simpson,solve_ivp
from scipy.special import kv
from model import ROOT,Object,G,LAM,HD,MH,VH,MPLANCK,GEV2CM2

ME,MMU,MPI,ALPHA=.00051099895,.1056583755,.1349768,1/137.035999084
HBAR=6.582119569e-25
DM_ENTROPY=4.4e-10


def theta(v,k):return k*v*VH/(MH*MH-2*LAM*v*v)


def gamma_hss(v,k):
 ms2=2*LAM*v*v
 c=k*VH*(MH*MH+2*ms2)/(MH*MH-ms2)
 return c*c/(32*np.pi*MH)*np.sqrt(1-4*ms2/(MH*MH))


def gamma_sm(m,cgamma=10.,hadronic_factor=1.,upper_tail=False):
 m=np.asarray(m)
 out=np.zeros_like(m)
 for ml in (ME,MMU):
  out+=ml*ml*m/(8*np.pi*VH*VH)*np.maximum(1-4*ml*ml/m**2,0)**1.5
 out+=cgamma**2*ALPHA**2*m**3/(256*np.pi**3*VH**2)
 # Leading chiral scalar-current response; no hadronic decay opens for our s masses.
 out+=hadronic_factor*3/(32*np.pi*VH**2*m)*(2*m*m/9+11*MPI*MPI/9)**2*np.sqrt(np.maximum(1-4*MPI*MPI/m**2,0))
 return np.where(m>=1.,m,out) if upper_tail else out


def spectral_production(m,v,k):
 s=m*m;ms2=2*LAM*v*v;mx2=G*G*v*v;md2=(G*G+HD*HD)*v*v
 gx=G**4*v*v/(4*np.pi*m)*np.sqrt(np.maximum(1-4*mx2/s,0))
 gd=HD**4*v*v/(8*np.pi*m)*np.sqrt(np.maximum(1-4*md2/s,0))
 gh=k*k*VH*VH/(32*np.pi*m)*((s+2*ms2)/(s-ms2))**2*np.sqrt(np.maximum(1-4*ms2/s,0))
 mixing=(k*v*VH)**2
 return gh/(s-MH*MH)**2+mixing*(gx+gd)/((s-MH*MH)**2*(s-ms2)**2)


def cold(obj,v,k,cgamma=10.,hadronic_factor=1.,upper_tail=False,points=401,mass_points=5001):
 ti=.020;ge=10.75;ms=np.sqrt(2*LAM)*v
 n=np.linspace(0,np.log(ti/1e-5),points);t=ti*np.exp(-n)
 s=2*np.pi**2*ge*t**3/45;rho=np.pi**2*ge*t**4/30
 h=np.sqrt(rho/(3*MPLANCK**2))
 mix=theta(v,k)
 width=float(mix*mix*gamma_sm(ms,cgamma,hadronic_factor,False))
 qinv=width*ms**3*t/(2*np.pi**2)*kv(2,ms/t)
 ninv=width*ms*ms*t/(2*np.pi**2)*kv(1,ms/t)
 start=2*G*v
 thresholds=[2*G*v,2*np.sqrt(G*G+HD*HD)*v,2*ms,2*MMU,2*MPI,1.]
 mg=np.unique(np.r_[np.linspace(start,2.,mass_points),[x for x in thresholds if start<=x<=2.]])
 spec=gamma_sm(mg,cgamma,hadronic_factor,upper_tail)*spectral_production(mg,v,k)
 # Inclusive equilibrium SM reaction density from its scalar spectral function.
 qann=t/np.pi**3*simpson(mg[None,:]**5*spec[None,:]*kv(2,mg[None,:]/t[:,None]),x=mg,axis=1)
 nann=t/np.pi**3*simpson(mg[None,:]**4*spec[None,:]*kv(1,mg[None,:]/t[:,None]),x=mg,axis=1)
 q=qinv+qann
 budget=cumulative_simpson(q/(s*h),x=n,initial=0)
 rhomax=budget*s # all injected energy immediately treated as pressureless and stable
 mass=obj.mass(v);gq=-k*VH*mass/(2*LAM*v*v)
 # Count every collision in the included coherent elastic channels as
 # potentially damaging, with form factors and speeds at their upper bounds.
 # This does not bound uncomputed inelastic transitions of the soliton.
 rate=np.zeros_like(t)
 for mf in (ME,MMU):
  dens=4*mf*mf*t/(2*np.pi**2)*kv(2,mf/t)
  mom2=12*mf**3*t*t/(2*np.pi**2)*kv(3,mf/t)
  rate+=gq*gq*(mf/VH)**2/(np.pi*MH**4)*(mf*mf*dens+.5*mom2)
 dens=3*MPI*MPI*t/(2*np.pi**2)*kv(2,MPI/t)
 mom2=9*MPI**3*t*t/(2*np.pi**2)*kv(3,MPI/t)
 mom4=45*MPI**4*t**3/(2*np.pi**2)*kv(4,MPI/t)
 a,b=11*MPI*MPI/9,8/9
 rate+=hadronic_factor*gq*gq/(4*np.pi*MH**4*VH**2)*(a*a*dens+2*a*b*mom2+b*b*mom4)
 ah=ALPHA*cgamma/(8*np.pi*VH)
 photon_mom4=720*t**7/np.pi**2
 rate+=gq*gq*(4*ah)**2/(4*np.pi*MH**4)*photon_mom4
 rate+=64*8.6e-11*s*obj.sigma_n(v,k)/GEV2CM2
 prob=simpson(rate/h,x=n)
 ysingle=simpson(ninv/(s*h),x=n)
 ypair=simpson(2*nann/(s*h),x=n)
 worst=budget[-1]
 # At the nominal decay time, use the low-temperature radiation relation.
 td=np.sqrt(width*MPLANCK/(2*np.sqrt(np.pi**2*3.36/90)))
 filling=DM_ENTROPY/mass*s[0]*4*np.pi/3*(obj.row['R90_over_inverse_gv']/(G*v))**3
 row={'Q':obj.charge,'v_phi_GeV':v,'kappa':k,'maximum_bath_temperature_GeV':ti,
      's_mass_GeV':ms,'theta':mix,'gamma_h_to_ss_GeV':gamma_hss(v,k),
      's_total_width_in_stated_SM_response_GeV':width,'lifetime_seconds':HBAR/width,
      'gamma_coefficient_envelope':cgamma,'chiral_width_multiplier':hadronic_factor,'SM_width_bound_above_1GeV':upper_tail,
      'produced_single_scalar_yield':float(ysingle),'produced_pair_particle_yield':float(ypair),
      'all_energy_survives_matter_GeV_per_entropy_BOUND':float(worst),
      'injected_energy_over_reference_DM_energy_BOUND':float(worst/DM_ENTROPY),
      'max_hidden_tadpole_over_vacuum_curvature_BOUND':float(max(1.5*rhomax/(LAM*v**4))),
      'max_visible_tadpole_over_vacuum_curvature_BOUND':float(max(2*k*rho/(MH*MH*LAM*v*v))),
      'visible_any_collision_optical_depth_BOUND':float(prob),
      'collision_bound_scope':'Included coherent elastic lepton, pion, photon and baryon channels; inelastic object transitions are not calculated.',
      'decay_temperature_radiation_estimate_GeV':float(td),
      'injected_energy_over_radiation_at_decay_BOUND_ESTIMATE':float(4*worst/(3*td)),
      'initial_object_energy_R90_volume_fraction':float(filling),
      'max_mean_injected_energy_per_object_GeV_BOUND':float(mass*worst/DM_ENTROPY)}
 return row,np.column_stack([n,t,qinv,qann,budget,rate/h])


def hot(v,k,method='DOP853',rtol=1e-9):
 ti=80.;gvis=106.75;ghid=4.;width=gamma_hss(v,k)
 rhoi=np.pi**2*gvis*ti**4/30
 # Comoving energy variables avoid tiny-difference subtraction in the conservation check.
 def rhs(n,y):
  rv=rhoi*np.exp(-4*n)*y[0];rd=rhoi*np.exp(-4*n)*y[1]
  temp=(30*rv/(np.pi**2*gvis))**.25;h=np.sqrt((rv+rd)/(3*MPLANCK**2))
  nh=MH*MH*temp/(2*np.pi**2)*kv(2,MH/temp)
  q=width*MH*nh
  qvar=2*width*MH*temp/np.pi**2*kv(0,MH/temp)
  return [-q*np.exp(4*n)/(h*rhoi),q*np.exp(4*n)/(h*rhoi),qvar*np.exp(2*n)/(h*v*v)]
 n=np.linspace(0,np.log(80/.5),1001)
 sol=solve_ivp(rhs,(n[0],n[-1]),[1.,0.,0.],t_eval=n,method=method,rtol=rtol,atol=1e-15)
 if not sol.success:raise RuntimeError(sol.message)
 rv=rhoi*np.exp(-4*n)*sol.y[0];rd=rhoi*np.exp(-4*n)*sol.y[1]
 tv=(30*rv/(np.pi**2*gvis))**.25;th=(30*rd/(np.pi**2*ghid))**.25
 variance=sol.y[2]*v*v*np.exp(-2*n)
 cth=LAM/4+G*G/6+HD*HD/12
 thermal_ratio=cth*th*th/(LAM*v*v)
 collisionless_ratio=3*variance/(v*v)
 row={'v_phi_GeV':v,'kappa':k,'initial_SM_temperature_GeV':ti,
      'gamma_h_to_ss_GeV':width,'max_hidden_temperature_rapid_equilibrium_GeV':float(max(th)),
      'critical_temperature_leading_highT_GeV':float(v*np.sqrt(LAM/cth)),
      'max_thermal_curvature_ratio':float(max(thermal_ratio)),
      'max_collisionless_self_curvature_ratio':float(max(collisionless_ratio)),
      'maximum_hidden_to_visible_energy_ratio':float(max(rd/rv)),
      'comoving_energy_conservation_max_residual':float(max(abs(sol.y[0]+sol.y[1]-1))),
      'SM_temperature_at_peak_hidden_temperature_GeV':float(tv[np.argmax(th)])}
 return row,np.column_stack([n,tv,rd,th,variance,thermal_ratio,collisionless_ratio])


def run():
 det=json.loads((ROOT/'DETECTOR_RESULTS.json').read_text())
 result={'classification':'Conditional histories at supplied detector-normalized couplings; formation and earlier reheating are initial data. Total vacuum matching is supplied.',
         'reference_DM_mass_per_entropy_GeV':DM_ENTROPY,'cases':[]}
 for case in det['cases']:
  q,v,k=case['Q'],case['v_phi_GeV'],case['kappa_anchor'];obj=Object(q)
  c,raw=cold(obj,v,k)
  stress,_=cold(obj,v,k,cgamma=100,hadronic_factor=1000,upper_tail=True)
  h,hr=hot(v,k)
  row={'Q':q,'v_phi_GeV':v,'kappa_anchor':k,'cold':c,'cold_production_stress_envelope':stress,'hot':h}
  if q==39944 and v==1.:
   refine,_=cold(obj,v,k,points=801,mass_points=10001)
   alt,_=hot(v,k,method='Radau',rtol=2e-10)
   row['cold_energy_budget_refinement_relative']=refine['all_energy_survives_matter_GeV_per_entropy_BOUND']/c['all_energy_survives_matter_GeV_per_entropy_BOUND']-1
   row['hot_independent_method_temperature_relative']=alt['max_hidden_temperature_rapid_equilibrium_GeV']/h['max_hidden_temperature_rapid_equilibrium_GeV']-1
   np.savetxt(ROOT/'data/reference_cold_history.csv',raw,delimiter=',',header='N,T_GeV,Q_inverse_GeV5,Q_annihilation_GeV5,energy_per_entropy_upper_GeV,visible_collision_dTau_dN',comments='')
   np.savetxt(ROOT/'data/reference_hot_history.csv',hr,delimiter=',',header='N,T_SM_GeV,rho_hidden_GeV4,T_hidden_equilibrium_GeV,collisionless_variance_GeV2,thermal_curvature_ratio,collisionless_curvature_ratio',comments='')
  result['cases'].append(row)
  print(q,v,'k',k,'hot thermal',h['max_thermal_curvature_ratio'],'cold budget fraction',c['injected_energy_over_reference_DM_energy_BOUND'],'stress fraction',stress['injected_energy_over_reference_DM_energy_BOUND'],'collisions',stress['visible_any_collision_optical_depth_BOUND'],flush=True)
 (ROOT/'HISTORY_RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':run()
