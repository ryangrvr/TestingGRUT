#!/usr/bin/env python3
"""Limited classical radial checks; no full perturbation/quantum certificate."""
import json
import numpy as np
from scipy.integrate import simpson, cumulative_trapezoid
from scipy.sparse import diags, bmat
from scipy.sparse.linalg import LinearOperator, eigsh
from profiles import ROOT, solve, save


def main():
    profiles = json.loads((ROOT/'PROFILE_RESULTS.json').read_text())
    row = next(x for x in profiles['profiles'] if x['Q']==39944 and x['lambda_phi']==.021)
    a = np.loadtxt(ROOT/row['profile_file'],delimiter=',',skiprows=1)

    class ReadSeed:
        x=a[:,0]
        p=np.array([row['Omega']])
        def sol(self,r):
            j=cumulative_trapezoid(self.x**2*a[:,1]**2,self.x,initial=0)
            return np.array([np.interp(r,self.x,v) for v in [*a[:,1:].T,j]])

    s=solve(39944,.021,ReadSeed(),length=120,tolerance=2e-9)
    r=save(s,39944,.021,True)
    r['tolerance']=2e-9
    profiles['old_charge_refinement']=r
    (ROOT/'PROFILE_RESULTS.json').write_text(json.dumps(profiles,indent=2)+'\n')
    results=[]
    for q,lam in [(39944,.021),(1000000,.03)]:
        p=next(x for x in profiles['profiles'] if x['Q']==q and x['lambda_phi']==lam)
        a=np.loadtxt(ROOT/p['profile_file'],delimiter=',',skiprows=1)
        om=p['Omega']
        length=p['radial_domain']
        for n in [300,600,1000]:
            h=length/(n+1)
            r=h*np.arange(1,n+1)
            f=np.interp(r,a[:,0],a[:,1])
            ss=np.interp(r,a[:,0],a[:,3])
            off=-np.ones(n-1)/h**2
            aa=diags([off,2/h**2+ss**2+.3*f**2-om**2,off],[-1,0,1])
            bb=diags([off,2/h**2+lam/.01*(3*ss**2-1)+f**2,off],[-1,0,1])
            cc=diags(2*f*ss)
            matrix=bmat([[aa,cc],[cc,bb]],format='csr')
            norm=4*np.pi*simpson(a[:,0]**2*a[:,1]**2,x=a[:,0])
            vec=np.r_[np.sqrt(4*np.pi)*r*f,np.zeros(n)]
            op=LinearOperator(matrix.shape,matvec=lambda x:matrix@x+4*om**2/norm*h*vec*np.dot(vec,x),dtype=float)
            ev=float(eigsh(op,k=1,which='SA',tol=1e-10,maxiter=30000,v0=np.linspace(.5,1,2*n),return_eigenvectors=False)[0])
            results.append(dict(Q=q,lambda_phi=lam,N=n,lowest_fixed_charge_radial_amplitude_eigenvalue=ev))
    (ROOT/'RADIAL_RESULTS.json').write_text(json.dumps(dict(scope='Finite-grid classical fixed-charge radial amplitude Hessian only; no full dynamical or quantum certification.',rows=results),indent=2)+'\n')
    print(results)


if __name__=='__main__':
    main()
