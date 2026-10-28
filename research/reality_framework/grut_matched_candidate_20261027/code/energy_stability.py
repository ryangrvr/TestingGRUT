"""Fixed-charge radial energy and translation-mode controls on the NEW profile.

Finite grids provide numerical evidence, not a certified continuum/full-gauge spectrum.
"""
import json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.sparse import diags,bmat
from scipy.sparse.linalg import eigsh,LinearOperator
from deformation import ROOT,G,LX,LAM


def main():
    a=np.loadtxt(ROOT/'PROFILE_M1.csv',delimiter=',',skiprows=1)
    f=CubicSpline(a[:,0],a[:,1]);s=CubicSpline(a[:,0],a[:,3])
    d=json.loads((ROOT/'DEFORMATION_RESULTS.json').read_text());omega=d['primary']['Omega'];eps=.002
    rows=[]
    for radius,dr in ((80.,.2),(80.,.1),(80.,.06666666666666667),(120.,.1)):
        n=int(round(radius/dr))-1;r=np.arange(1,n+1)*dr
        ff=f(r);ss=s(r)
        lap=diags([-np.ones(n-1)/dr**2,2*np.ones(n)/dr**2,-np.ones(n-1)/dr**2],[-1,0,1])
        vff=ss*ss+3*LX/G**2*ff*ff-omega*omega
        vss=LAM/G**2*(3*ss*ss-1)-6*eps/G**2*(ss-1)+ff*ff
        mix=diags(2*ff*ss)
        base=bmat([[lap+diags(vff),mix],[mix,lap+diags(vss)]],format='csr')
        w=np.r_[r*ff,np.zeros(n)]
        rank_coefficient=4*omega*omega/(np.sum((r*ff)**2)*dr)
        op=LinearOperator((2*n,2*n),matvec=lambda x:base@x+rank_coefficient*dr*w*np.dot(w,x),dtype=float)
        radial=np.sort(eigsh(op,k=3,which='SA',tol=2e-9,return_eigenvectors=False))
        ell1=base+diags(np.tile(2/r**2,2))
        translation=np.sort(eigsh(ell1,k=3,which='SA',tol=2e-9,return_eigenvectors=False))
        phase=lap+diags(ss*ss+LX/G**2*ff*ff-omega*omega)
        phase_eigs=np.sort(eigsh(phase,k=2,which='SA',tol=2e-9,return_eigenvectors=False))
        rows.append({'radius':radius,'dr':dr,'radial_fixed_charge_lowest':radial.tolist(),
                     'ell1_lowest_translation_expected_zero':translation.tolist(),
                     'phase_lowest_U1_expected_zero':phase_eigs.tolist()})
        print('STABILITY',radius,dr,radial[0],translation[0],phase_eigs[0],flush=True)
    (ROOT/'STABILITY_RESULTS.json').write_text(json.dumps({'scope':__doc__,'rows':rows,'full_gauge_or_quantum_spectrum':False},indent=2)+'\n')


if __name__=='__main__':main()
