"""Even axisymmetric M1 pulse in Legendre coefficients on a radial grid."""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.special import roots_legendre,eval_legendre
from initial import geometry,values,OMEGA

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser();p.add_argument('--h',type=float,required=True);p.add_argument('--L',type=float,default=200.)
    p.add_argument('--lmax',type=int,default=4);p.add_argument('--nq',type=int,default=12)
    p.add_argument('--amp',type=float,default=.05);p.add_argument('--quad',type=float,default=1.);p.add_argument('--tag',required=True);a=p.parse_args()
    if a.lmax<0 or a.lmax%2:raise ValueError('lmax must be nonnegative and even')
    if a.nq<2*a.lmax+1:raise ValueError('Gauss nodes must integrate quartic modal products')
    r,vol,_,_=geometry(a.h,a.L);f,s=values(r);ls=np.arange(0,a.lmax+1,2);nm=len(ls);shape=np.exp(-(r/12)**2)
    x=np.zeros((len(r),nm));y=np.zeros_like(x);ss=np.zeros_like(x)
    x[:,0]=f*(1+a.amp*shape);ss[:,0]=s+a.amp*(1-s)*shape
    if nm>1:
        x[:,1]=a.amp*a.quad*f*shape*r*r/(r*r+4)
        ss[:,1]=a.amp*a.quad*(1-s)*shape*r*r/(r*r+4)
    px=np.zeros_like(x);py=-OMEGA*x;ps=np.zeros_like(x);weights=1/(2*ls+1)
    Q=np.sum(vol[:,None]*(-x*py)*weights[None,:]);py*=2000/Q
    path=ROOT/'.work'/(a.tag+'_initial.csv')
    np.savetxt(path,np.column_stack((r,x,y,ss,px,py,ps)),delimiter=',',comments='',header='r,'+','.join(f'{v}{l}' for v in ('x','y','s','px','py','ps') for l in ls))
    mu,w=roots_legendre(a.nq);table=np.column_stack((mu,w,*[eval_legendre(int(l),mu) for l in ls]))
    np.savetxt(ROOT/'.work'/(a.tag+'_angles.csv'),table,delimiter=',',comments='',header='mu,w,'+','.join('P'+str(l) for l in ls))
    out={'configuration':vars(a),'modes':ls.tolist(),'renormalized_charge':float(np.sum(vol[:,None]*(-x*py)*weights)/.01),
         'scope':'Supplied axisymmetric even-Legendre control; quadrupole differs from the earlier triaxial (x^2-y^2) 3D perturbation.'}
    (ROOT/'.work'/(a.tag+'_initial.json')).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
