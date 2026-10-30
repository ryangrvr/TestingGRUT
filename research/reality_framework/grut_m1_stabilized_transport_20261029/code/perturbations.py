"""Declared de Sitter vacuum to spectator charge transfer.

The inherited canonical 12x12 generator contains scalar/longitudinal mixing
but excludes metric/inflaton/radiation perturbations. k=0 uses its e=0 global
scalar limit restricted to five physical zero-gauge-charge modes. Finite k,e
controls quantify that approximation on the tested scales. The resulting
fixed-background charge power is not the full gauge-invariant CDI observable.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
import model as m
from trajectory import run
ROOT=m.ROOT
J=np.block([[np.zeros((6,6)),np.eye(6)],[-np.eye(6),np.zeros((6,6))]])

def de_sitter(m2,x0,points):
 n=len(m2);b=m2-2;w=np.sqrt(1+b/x0**2);dw=-b/(w*x0**3);u=1/np.sqrt(w)+0j;du=(1j*w-.5*dw/w)*u
 initial=np.r_[x0*u,-x0*(u+x0*du)]
 sol=solve_ivp(lambda N,y:np.r_[y[n:],-3*y[n:]-(np.exp(-2*N)+m2)*y[:n]],(-np.log(x0),max(points)),initial,method='DOP853',rtol=2e-10,atol=2e-12,t_eval=points)
 if not sol.success:raise RuntimeError(sol.message)
 return sol.y[:n],sol.y[n:]

def charge_rows(tag):
 a=np.load(ROOT/(tag+'.npz'));times=a['times'];bg=a['bg'];TT=a['trans'];rows=[];gauge=[];sym=[]
 for v,trans in zip(bg,TT):
  aa=np.exp(v[3]);H=np.sqrt((v[1]**2/2+4.5*v[0]**2+v[2])/3);chi=aa*v[4:10];chip=aa*aa*(v[10:16]+H*v[4:10]);row=np.zeros(12);row[:2]=[-chip[1],chip[0]];row[6:8]=[chi[1],-chi[0]];rows.append(row@trans)
  gg=np.r_[m.T@chi,m.T@chip];gauge.append(float(abs(row@gg)))
  sym.append([float(np.max(abs(t.T@J@t-J))) for t in trans])
 rows=np.array(rows);mask=times>=.8*times[-1];averaged=np.trapezoid(rows[mask],times[mask],axis=0)/(times[mask][-1]-times[mask][0]);nx=float(np.trapezoid(a['obs'][mask,4],times[mask])/(times[mask][-1]-times[mask][0]))
 return averaged/nx,{'gauge_charge_response_max':max(gauge),'symplectic_absolute_end':sym[-1],'late_mean_NX':nx},times,rows/nx

def main():
 z=m.initial();vals,U=np.linalg.eigh(m.hessian(z));mask=vals>1e-5;physical=U[:,mask];masses=vals[mask];row,checks,t,rr=charge_rows('holder');refrow,refchecks,_,_=charge_rows('refined');frow,fchecks,_,_=charge_rows('finite_k');row=row[0];refrow=refrow[0]
 # Both q and physical velocity projections in canonical coordinates. Initial
 # p=q+dr/d(H_i t), since a=H/H_i=1 and the zero-velocity background has g'=g.
 projection=np.block([[physical,np.zeros_like(physical)],[physical,physical]])
 long=row@projection;longref=refrow@projection;fin=frow@projection
 rows=[];points=np.array([0.,2.,4.,8.]);S,D=de_sitter(masses,1000.,points);S2,D2=de_sitter(masses,2000.,points)
 for N in [0.,2.,4.,8.,20.,50.,60.]:
  if N<=8:
   j=np.argmin(abs(points-N));s=S[:,j];v=D[:,j];ss=S2[:,j];vv=D2[:,j]
  else:
   dt=N-8;omega=np.sqrt(masses-2.25);s=np.exp(-1.5*dt)*(S[:,-1]*np.cos(omega*dt)+(D[:,-1]+1.5*S[:,-1])/omega*np.sin(omega*dt));v=-1.5*s+np.exp(-1.5*dt)*(-omega*S[:,-1]*np.sin(omega*dt)+(D[:,-1]+1.5*S[:,-1])*np.cos(omega*dt));ss=vv=None
  init=np.vstack([physical*s,physical*(s+v)]);response=row@init;power=float(np.sum(abs(response)**2)/(4*np.pi**2));refpower=float(np.sum(abs(refrow@init)**2)/(4*np.pi**2));finitepowers=np.sum(abs(frow@init)**2,axis=1)/(4*np.pi**2)
  check0=None if ss is None else float(power/(np.sum(abs(row@np.vstack([physical*ss,physical*(ss+vv)]))**2)/(4*np.pi**2))-1)
  rows.append({'assumed_efolds_before_reheating':N,'spectator_fractional_charge_power':power,'sqrt_power':float(np.sqrt(power)),'trajectory_tolerance_relative_difference':power/refpower-1,
    'finite_k_charge_power_ratio':(finitepowers/power).tolist(),'vacuum_depth_relative_difference':check0})
 # Perturb the supplied Sigma phase, keeping initial magnitudes and velocities.
 # This tangent is distinct from a fully covariant cosmological mode.
 tangent=np.zeros(12);tangent[4:6]=m.real(np.array([0j,0j,1j*z[2]]))[4:6];canonical=np.r_[tangent[:6],tangent[:6]+tangent[6:]];slope=float(row@canonical)
 fd=[]
 for step in [.001,.0005]:
  values=[]
  for sign in [-1,1]:
   zz=z.copy();zz[2]*=np.exp(1j*sign*step);tag=f'phase_{sign:+d}_{step:g}';o=run(tag,800,2e-10,ks=[],z=zz);values.append(o['late_mean_NX'])
  fd.append({'step':step,'dln_N_dphase':float((np.log(abs(values[1]))-np.log(abs(values[0])))/(2*step))})
 ref=json.loads((ROOT/'reference.json').read_text())['late_mean_NX'];nx=checks['late_mean_NX'];ratio=nx/ref;scale=(1/abs(ratio))**(2/3);mass_scenarios={'original_mass_per_charge':0.,'hypothetical_fully_cooled':.008409936530420273,'checkpoint_total':.008944364105298863}
 calibration={k:{'mass_fraction_increase':v,'H_ratio':float((abs(ref/nx)/(1+v))**(2/3)),'new_H_GeV':float(m.H*(abs(ref/nx)/(1+v))**(2/3))} for k,v in mass_scenarios.items()}
 out={'canonical_mass2_over_H2':masses.tolist(),'late_fractional_charge_response_row':row.tolist(),'physical_initial_position_velocity_response':long.tolist(),
 'checks':checks,'refined_checks':refchecks,'finite_k_checks':fchecks,'refined_projected_response_relative_norm':float(np.linalg.norm(long-longref)/np.linalg.norm(long)),
 'finite_k_projected_response_relative_norms':(np.linalg.norm(fin-long,axis=1)/np.linalg.norm(long)).tolist(),'finite_k_over_Hi':[.001,.0001],
 'phase_sensitivity_tangent':slope,'phase_finite_differences':fd,'vacuum_rows':rows,'holder_charge_relative_to_original':float(ratio),'fixed_dimensionless_family_calibration':calibration,
 'scope':'Spectator scalar/longitudinal transfer on fixed reheating slicing with an assumed de Sitter Bunch-Davies state and duration. Metric/inflaton/radiation perturbations, delta R, entropy perturbations, curvature correlation, formation response and the observational transfer are absent. These are charge-power proxies, not actual P_CDI or a Planck likelihood. Calibration co-scales all heavy masses and decay rates with H_i; light mass and dimensionless couplings stay fixed.'}
 (ROOT/'perturbations.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');np.savez_compressed(ROOT/'charge_transfer.npz',t=t,rows=rr,averaged_response=row,physical_basis=physical,eigenvalues=masses);print(json.dumps(out,indent=2))
if __name__=='__main__':main()
