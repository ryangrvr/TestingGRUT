"""Independent fixed-charge second variation, reduced radial self-adjoint form.

This does not evolve a spherical truncation in lieu of the 3D experiment.
The angular barrier orders ALL harmonics; numerical radial eigenvalues bound
the remaining candidate directions, modulo symmetry zero modes.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import hashlib, json, time
import numpy as np
from scipy.integrate import solve_bvp, simpson
from scipy.sparse import diags, bmat
from scipy.sparse.linalg import eigsh, LinearOperator

P=Path(__file__).resolve().parent
OMEGA=.8000000019763136
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def profile(om,initial):
    r=initial[0];y=initial[1]
    # y = F,F',S,S'; radial regularity supplied by solve_bvp's singular term.
    def rhs(r,y):
        F,Fp,S,Sp=y
        return np.array([Fp,(S*S+.1*F*F-om*om)*F,Sp,(S*S-1+F*F)*S])
    def bc(ya,yb):return np.array([ya[1],ya[3],yb[0],yb[2]-1])
    singular=np.diag([0,-2,0,-2])
    sol=solve_bvp(rhs,bc,r,y,S=singular,tol=2e-8,max_nodes=20000)
    if sol.status:raise RuntimeError(sol.message)
    rr=np.linspace(0,r[-1],10001);F,Fp,S,Sp=sol.sol(rr)
    I=4*np.pi*simpson(rr*rr*F*F,x=rr);Q=om*I
    E=4*np.pi*simpson(rr*rr*(.5*om*om*F*F+.5*(Fp*Fp+Sp*Sp)+.25*(S*S-1)**2+.5*S*S*F*F+.025*F**4),x=rr)
    return sol,{'Omega':om,'Qbar':float(Q),'Ebar':float(E),'E_over_Q':float(E/Q),'bvp_nodes':len(sol.x),'max_rms_residual':float(max(sol.rms_residuals)),'Ibar':float(I),'minimum_F':float(min(F)),'maximum_positive_Fprime':float(max(Fp)),'minimum_Sprime':float(min(Sp))}
def eig(A,k=3):
    n=A.shape[0];v0=np.cos(np.arange(n)*.013)+.1
    vals,vec=eigsh(A,k=k,which='SA',tol=2e-9,maxiter=100000,v0=v0)
    order=np.argsort(vals);vals=vals[order];vec=vec[:,order]
    residuals=[float(np.linalg.norm(A@vec[:,i]-vals[i]*vec[:,i])) for i in range(k)]
    return vals.tolist(),residuals
def spectrum(sol,R,h):
    n=int(round(R/h))-1;r=np.arange(1,n+1)*h;F,Fp,S,Sp=sol.sol(r)
    K=diags([-np.ones(n-1)/h**2,2*np.ones(n)/h**2,-np.ones(n-1)/h**2],[-1,0,1],format='csr')
    coupl=diags(2*S*F,format='csr')
    def plus(l):
        angular=l*(l+1)/r**2
        return bmat([[K+diags(S*S+.3*F*F-OMEGA**2+angular),coupl],[coupl,K+diags(3*S*S-1+F*F+angular)]],format='csr')
    minus0=K+diags(S*S+.1*F*F-OMEGA**2)
    H=plus(0);f=np.r_[r*F,np.zeros(n)];J=float(h*np.dot(r*F,r*F));alpha=4*OMEGA**2*h/J
    HQ=LinearOperator(H.shape,matvec=lambda v:H@v+alpha*f*np.dot(f,v),dtype=np.float64)
    out={'R':R,'dr':h,'dimension':2*n}
    for name,A in [('phase_l0',minus0),('unconstrained_l0',H),('fixed_charge_l0',HQ),('translation_l1',plus(1)),('shape_l2',plus(2))]:
        out[name],out[name+'_residuals']=eig(A)
    tr=np.r_[r*Fp,r*Sp];phase=r*F
    out['translation_generator_residual']=float(np.linalg.norm(plus(1)@tr)/np.linalg.norm(tr))
    out['phase_generator_residual']=float(np.linalg.norm(minus0@phase)/np.linalg.norm(phase))
    return out
if __name__=='__main__':
    protocol=P/'SPECTRAL_PROTOCOL.json'
    if not protocol.exists():
      spec={'local_freeze_unix':time.time(),'public_preregistration':False,'source_sha256':sha(Path(__file__)),'grids':[[50,.1],[50,.05],[50,.025],[60,.05]],'omega':OMEGA,'criteria':'Fixed-Q l=0 lowest eigenvalue positive and convergent; l=1 translation and l=0 phase approach zero as O(dr^2); l=2 positive. Inspect all residuals. Angular operator ordering extends l>=2 and l>=1 phase sectors. No claim of computer-assisted rigorous spectral enclosure.','mathematics':'Hminus_l=-d_r^2+l(l+1)/r^2+S^2+b F^2-Omega^2; Hplus_l diagonals add 2bF^2 and a(3S^2-1)+F^2, off-diagonal 2SF. Fixed Q adds 4Omega^2 |rF><rF| / integral(r^2F^2 dr) to l=0 amplitude block.','posthoc_scope':'Additive diagnostic proposed after coarse/medium finite runs; no empirical target fitted. Profile and exact symmetry modes predate these runs.'}
      protocol.write_text(json.dumps(spec,indent=2)+'\n')
    spec=json.loads(protocol.read_text())
    if spec['source_sha256']!=sha(Path(__file__)):raise RuntimeError('spectral source changed')
    raw=np.loadtxt(P/'profile.dat');initial=(raw[:,0],raw[:,[1,2,3,4]].T)
    sol,base=profile(OMEGA,initial)
    # An enlarged boundary is re-solved, not a polynomial extrapolation.
    r60=np.linspace(0,60,1201);guess=np.zeros((4,len(r60)));guess[2]=1;mask=r60<=50;guess[:,mask]=sol.sol(r60[mask]);sol60,large=profile(OMEGA,(r60,guess))
    spectra=[]
    for R,h in spec['grids']:
      v=spectrum(sol60 if R==60 else sol,R,h);spectra.append(v);print('SPECTRAL',R,h,'Q-l0',v['fixed_charge_l0'][0],'l1',v['translation_l1'][0],flush=True)
    neighbours=[]
    for om in [OMEGA-.002,OMEGA+.002]:neighbours.append(profile(om,(sol.x,sol.y))[1])
    result={'base':base,'larger_domain':large,'radial_spectra':spectra,'nearby_profiles':neighbours,'dQ_dOmega_centered':(neighbours[1]['Qbar']-neighbours[0]['Qbar'])/.004,'interpretation':'Numerical local constrained-energy minimum at one isolated stationary branch, modulo phase and translations, supported by analytic angular ordering; not a global-minimizer certificate or a finite-collision Lyapunov exponent.'}
    (P/'SPECTRAL_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
