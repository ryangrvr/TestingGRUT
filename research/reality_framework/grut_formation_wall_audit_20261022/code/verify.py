"""Separate full linear equations and same-author isolated numerical replay."""
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigvals
from scipy.optimize import minimize_scalar
import calculate as calc


def leaves(x,key=''):
 if isinstance(x,dict):
  for k,v in x.items():yield from leaves(v,key+'/'+k)
 elif isinstance(x,list):
  for i,v in enumerate(x):yield from leaves(v,key+'/'+str(i))
 elif isinstance(x,(int,float)) and not isinstance(x,bool):yield key,float(x)


def run(replay=False):
 root=calc.ROOT;out=json.loads((root/'RESULTS.json').read_text())
 # A denser control avoids machine cancellation of the tiny physical
 # phonon frequency when independently diagonalizing the full 6x6 system.
 tests=[]
 for f2 in [1e-5,1e-4,1e-3]:
  w=np.sqrt(calc.G**2-calc.ATTR*f2);phi2=1-calc.G**2*f2/calc.LAM
  number=w*f2
  for y in [.2,1.,1.8,2.2]:
   k2=calc.ATTR*f2*y
   a=2*calc.LX*f2;b=2*calc.LAM*phi2;c=2*calc.G**2*np.sqrt(phi2*f2)
   mat=np.zeros((6,6));mat[:3,3:]=np.eye(3)
   mat[3,0]=-(k2+a);mat[3,2]=-c;mat[3,4]=2*w
   mat[4,1]=-k2;mat[4,3]=-2*w
   mat[5,0]=-c;mat[5,2]=-(k2+b)
   actual=max(eigvals(mat).real);expected=calc.growth(number,1.,y)
   if y<2:assert abs(actual/expected-1)<2e-7,(f2,y,actual,expected)
   else:assert abs(actual)<1e-12
   tests.append({'f_squared':f2,'y':y,'full_matrix_growth':float(actual),'determinant_growth':expected})
 # Independently evolve one full linear perturbation on the same denser
 # stationary background, with both radial and phase equations retained.
 f2=1e-3;w=np.sqrt(calc.G**2-calc.ATTR*f2);phi2=1-calc.G**2*f2/calc.LAM
 n=w*f2;k2=calc.ATTR*f2;rate=calc.growth(n,1.,1.)
 a=2*calc.LX*f2;b=2*calc.LAM*phi2;c=2*calc.G**2*np.sqrt(phi2*f2)
 p=-2*w*rate/(k2+rate*rate);u=-c/(k2+b+rate*rate)
 initial=np.array([1,p,u,rate,rate*p,rate*u])
 def rhs(t,z):
  return [z[3],z[4],z[5],-(k2+a)*z[0]-c*z[2]+2*w*z[4],-k2*z[1]-2*w*z[3],-(k2+b)*z[2]-c*z[0]]
 sol=solve_ivp(rhs,(0,5/rate),initial,method='Radau',rtol=2e-9,atol=1e-12)
 ode_error=float(sol.y[0,-1]/np.exp(5)-1)
 assert sol.success and abs(ode_error)<1e-6
 ref=next(r for r in out['cases'] if r['Q']==39944 and r['v_phi_GeV']==1.)
 n=ref['linear_condensate']['supplied_charge_density_GeV3']
 dense_y=np.linspace(.95,1.05,10001)
 peak=max(calc.growth(n,1.,y) for y in dense_y)
 assert abs(peak/ref['linear_condensate']['maximum_growth_GeV']-1)<1e-10
 # Thermodynamic small-k comparison improves with decreasing k.
 convergence=[]
 fsq,w,_=calc.condensate(n,1.)
 cs=ref['linear_condensate']['cs_squared_from_charge_compressibility']
 for y in [1e-3,1e-5,1e-7]:
  value=-calc.growth(n,1.,y)**2/(calc.ATTR*fsq*y)
  convergence.append(float(abs(value/cs-1)))
 assert convergence[2]<1e-7 and convergence[2]<convergence[1]<convergence[0]
 # Conservation identity for a supplied carrier population at each case.
 for r in out['cases']:
  assert np.isclose(r['supplied_charge_yield']*r['free_carrier_mass_GeV']-calc.DM_ENTROPY,r['binding_release_GeV_per_visible_entropy'],rtol=1e-12)
  assert r['NR_Yukawa_zero_energy_Birman_Schwinger_trace']<1
  assert r['empty_cold_portal_max_DM_fraction_enlarged_response']<1.2e-7
  assert r['wall_tension_over_2026_081MeV_cubed']>1e7
 result={'scope':'Same author; full independent linear matrix and ODE organizations, not a nonlinear formation simulation.',
  'matrix_controls':tests,'full_linear_ODE_five_efolds_relative_residual':ode_error,
  'small_k_thermodynamic_convergence_relative_errors':convergence,
  'physical_peak_dense_grid_relative_residual':float(peak/ref['linear_condensate']['maximum_growth_GeV']-1)}
 if replay:
  with tempfile.TemporaryDirectory(prefix='grut-formation-replay-') as path:
   base=Path(path);packet=base/root.name;emp=base/calc.EMP.name
   (packet/'code').mkdir(parents=True);emp.mkdir()
   shutil.copy2(root/'code/calculate.py',packet/'code/calculate.py')
   for name in ['DETECTOR_RESULTS.json','HISTORY_RESULTS.json']:shutil.copy2(calc.EMP/name,emp/name)
   subprocess.run([sys.executable,str(packet/'code/calculate.py')],env=dict(os.environ,OPENBLAS_NUM_THREADS='1'),check=True,stdout=subprocess.DEVNULL)
   old=dict(leaves(out));new=dict(leaves(json.loads((packet/'RESULTS.json').read_text())))
   assert old.keys()==new.keys()
   for k in old:assert np.isclose(old[k],new[k],rtol=3e-10,atol=1e-35),(k,old[k],new[k])
   for name in ['reference_dispersion.csv','neutral_wall_profile.csv']:
    assert np.allclose(np.loadtxt(root/'data'/name,delimiter=',',skiprows=1),np.loadtxt(packet/'data'/name,delimiter=',',skiprows=1),rtol=3e-10,atol=1e-35)
   result['fresh_replay']={'numeric_leaves':len(old),'CSV_arrays':2,'same_author':True}
 (root/'VERIFICATION_RESULTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k!='matrix_controls'},indent=2))


if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--replay',action='store_true');run(p.parse_args().replay)
