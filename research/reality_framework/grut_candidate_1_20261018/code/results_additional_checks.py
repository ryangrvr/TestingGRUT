"""Supplementary diagnostics on the saved, frozen final profile; no re-fit."""
import sys,json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.sparse import diags,bmat
from scipy.sparse.linalg import eigsh
sys.path.insert(0,str(Path(__file__).parent))
import qball
out=Path(__file__).parent
p=np.load(out/'results_profile.npz')
class Saved: pass
sol=Saved();sol.x=p['rho'];sol.p=np.array([float(p['Omega'])]);sol.sol=CubicSpline(p['rho'],p['y'],axis=1)
checks=qball.checks(sol)
N=2000;L=40.;h=L/(N+1);r=h*np.arange(1,N+1);F,_,S,_=sol.sol(r);O=sol.p[0]
A=diags([-np.ones(N-1)/h**2,2/h**2+S*S+.3*F*F-O*O,-np.ones(N-1)/h**2],[-1,0,1])
B=diags([-np.ones(N-1)/h**2,2/h**2+3*S*S-1+F*F,-np.ones(N-1)/h**2],[-1,0,1])
C=diags(2*F*S)
H=bmat([[A,C],[C,B]],format='csr')
local=eigsh(H,k=3,which='SA',tol=1e-10,return_eigenvectors=False)
extra={'unconstrained_local_l0_lowest':np.sort(local),'local_l0_grid':{'N':N,'L':L,'h':h},'checks':checks,'note':'Local l0 excludes the fixed-charge rank-one term; it is diagnostic, not the physical fixed-Q Hessian.'}
(out/'results_additional_checks.json').write_text(json.dumps(extra,indent=2,default=qball.serial)+'\n')
result=json.loads((out/'results.json').read_text());result['checks']=checks;result['configuration']=qball.CONFIG;result['additional_diagnostics']=extra
(out/'results.json').write_text(json.dumps(result,indent=2,default=qball.serial)+'\n')
(out/'results_prescription.json').write_text(json.dumps({k:result[k] for k in ['configuration','versions','command','preliminary_failure_archives','preparation_search']},indent=2,default=qball.serial)+'\n')
print('phase l0 refinements:',[(a['N'],a['L'],a['lowest_eigenvalues'][0]) for a in checks['phase_sector_spectrum'] if a['ell']==0])
print('local l0 without fixed charge:',np.sort(local))
